'''
Usage:
    pipeline.py <r1> <r2> --config <config> [-o <outdir>] [--log <log>]

Options:
    -c <config>, --config <config>
    --log <log>
'''
import sh
import itertools
from toolz.dicttoolz import keymap,valfilter,keyfilter
from docopt import docopt
from path import Path
import types
import yaml
import sys
from functools import partial
import shutil
import os
# TODO: LZW!
# TODO: Log commands as doing them
# TODO: Remove log file at start if it exists.


#############
# Utilities #
#############
class Config:
    def __init__(self, entries):
        self.__dict__.update(entries)
        for k,v in self.__dict__.items():
            if type(v) == types.DictType:
                setattr(self, k, Config(v))
class Sh_(object):
    def __getattr__(self, attr):
        #cmd = getattr(sh, attr)
        def command(*args, **kwargs):
            #fixedargs = keymap("-{}".format, kwargs)
            bools = valfilter(lambda x: type(x) is bool, kwargs)
            vargs = keymap("-{}".format, keyfilter(lambda x: x not in ['_err', '_out'], valfilter(lambda x: not type(x) is bool, kwargs)))
            #bools.update(vargs)
            fixedargs = itertools.chain(vargs.items())
            getattr(sh, attr)(*(list(args) + list(fixedargs)), **bools)
        return command
sh_ = Sh_()

############
# Parts    #
############

def star(log, cfg, in1, in2):
  sh.STAR('--readFilesIn', in1, in2, #readFilesIn=unlist(in1, in2),
          genomeDir=cfg.star.starDB,
          outSAMtype="SAM",
          outReadsUnmapped="Fastx",
          _out=log, _err=log)

def pricefilter(log, cfg, in1, in2, o1, o2):
    cfg = cfg.pricefilter
    sh_.PriceSeqFilter('-fp', in1, in2,
                       '-op', o1, o2,
                       '-rqf', cfg.highQualPercent, cfg.highQualMin,
                       rnf=cfg.calledPercent)

def cdhitdup(log, cfg, r1, r2, o1, o2):
    sh_.cd_hit_dup(i=r1, i2=r2, o=o1, o2=o2, e=cfg.cdhitdup.minDifference, _err=log, _out=log)

# LZW!

def bowtie_sensitive(log, cfg, r1, r2, o1):
    args = {'1' : r1, '2' : r2,
                'very_sensitive_local' : True,
                'un_conc' : Path(o1).splitext()[0],
                'x' : cfg.bowtie2.bowtieDB,
            '_err' : log, '_out' : log}
    sh.bowtie2(**args)

def rapsearch(log, cfg, fq, out):
    sh.rapsearch(o=out, d=cfg.rapsearch.rapsearchDB, q=fq, _err=log, _out=log)

def blastn(log, cfg, fq, out):
    sh_.blastn(outfmt=8, db=cfg.ncbi.ntDB, q=fq, _err=log, _out=out)

def krona(log, cfg, blast, out):
    sh.ktImportBlast(blast, o=out, _err=log, _out=log) # probably need config for kronadb!


############
# Pipeline #
############


def run(cfg, input1, input2, log=None):
  p = partial(os.path.join, cfg.outdir)
  _star1 = p("Unmapped.out.mate1")
  _star2 = p("Unmapped.out.mate2")
  star1 = p("star.r1.fq")
  star2 = p("star.r2.fq")
  psf1 =  p( "psf.r1.fq"             )
  psf2 =  p( "psf.r2.fq"             )

  cd1 =       p( "cd.r1.fq" )
  cd2 =       p( "cd.r2.fq" )

  bowtie1 =   p( "bowtie.r1.fq" )
  bowtie2 =   p( "bowtie.r2.fq" )
  nr1     =   p( "rapsearch.r1.blast" )
  nr2     =   p( "rapsearch.r2.blast" )

  nt1 =       p( "r1.blast" )
  nt2 =       p( "r2.blast" )
  kronaNT1  = p( "r1.NT.html" )
  kronaNT2  = p( "r2.NT.html" )
  kronaNR1  = p( "r1.NR.html" )
  kronaNR2  = p( "r2.NR.html" )

  if not log:
    log = sys.stdout

  need = lambda p: not os.path.exists(p)

  if need(_star1):
    star(log, cfg, input1, input2)

  if need(star1):
    shutil.copy(_star1, star1)
    shutil.copy(_star2, star2)

  if need(psf1):
    pricefilter(log, cfg, star1, star2, psf1, psf2)

  if need(cd1):
    cdhitdup(log, cfg, psf1, psf2, cd1, cd2)

  if need(bowtie1):
    bowtie_sensitive(log, cfg, cd1, cd2, bowtie1)

  if need(nt1):
    blastn(log, cfg, bowtie1, nt1)
    blastn(log, cfg, bowtie2, nt2)

  if need(nr1):
    rapsearch(log, cfg, bowtie1, nr1)
    rapsearch(log, cfg, bowtie2, nr2)

  if need(kronaNT1):
    krona(log, cfg, nt1, kronaNT1)
    krona(log, cfg, nt2, kronaNT2)
    krona(log, cfg, nr1, kronaNR1)
    krona(log, cfg, nr2, kronaNR2)

def main():
  args = docopt(__doc__, version='Version 1.0')
  cfg = args['--config']
  cfg = yaml.load(open(cfg))
  cfg = Config(cfg)
  cfg.outdir = args['-o'] or "."
  if args['--log']:
    log = open(args['--log'], 'a')
  else:
    log = sys.stdout
  try:
    run(cfg, args['<r1>'], args['<r2>'], log)
  except Exception as e:
    log.write(str(e))
  if args['--log']: log.close()
  sys.exit(0)

if __name__ == '__main__': main()
