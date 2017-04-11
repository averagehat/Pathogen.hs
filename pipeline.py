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
from Bio import SeqIO
import subprocess
# TODO: LZW!
# TODO: Log commands as doing them
def lzw(sequence):
# https://github.com/betegonm/gen/blob/64aef21cfeefbf27b1e2bd6587c555d4df4f6913/gen.py#L294
  output = []
  table = dict(dict((chr(i), i) for i in range(256)))
  s = ''
  for ch in sequence:
    it = s + ch
    if it in table:
      s = it
    else:
      output.append(table[s])
      table[it] = len(table)
      s = ch
  output.append(table[s])
  return len(output)

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
    out = out.splitext()[0] # rapsearch adds an m8 extension
    sh.rapsearch(o=out, d=cfg.rapsearch.rapsearchDB, q=fq, _err=log, _out=log)

def blastn(log, cfg, fq, out):
    print "attempting blast with %s %s" % (fq, out)
    #sh_.blastn(outfmt=6, db=cfg.ncbi.ntDB, query=fq, _err=log, _out=out)
    sh.blastn(outfmt=6, db=cfg.ncbi.ntDB, query=fq, _err=log, _out=out, _long_prefix='-')

def krona(log, cfg, blast, out):
    sh.ktImportBLAST(blast, o=out, _err=log, _out=log) # probably need config for kronadb!


def blastx(log, cfg, fq, out):
    sh.blastx(outfmt=6, db=cfg.ncbi.nrDB, query=fq, _err=log, _out=out, _long_prefix='-')

def abyss(log, cfg, r1, r2, out):
    print 'trying :('
    dir = out.dirname()
    f1 = dir.relpathto(r1)
    f2 = dir.relpathto(r2)
    prefix=out.basename().split('-')[0]
    print f1, f2, 'name=%s' % prefix, 'k=%s' % 25
    sh.run_abyss(f1, f2, 'name=%s' % prefix, 'k=%s' % 25, C=dir, _err=log, _out=log)

def lzw_filter(log, cfg, r1, r2, out1, out2):
    un_comp_len = len(str(x.seq))
    comp_len = sum(imap(len, sh.gzip(f=True, _in=str(x.seq))))

def filter_pair_fastq(func, r1, r2, o1, o2):
    fwd = SeqIO.parse(r1, 'fastq')
    rev = SeqIO.parse(r2, 'fastq')
    filtered = ((x, y) for (x, y) in izip(fwd, rev)
                if func(x) and func(y))
    with open(o1, 'w') as f1, open(o2, 'w') as f2:  # because SeqIO is slow for writing single records
        for (x, y) in izip(fwd, rev):
            if func(x) and func(y):
                '@
                f1.write(x.format(form))
                f2.write(y.format(form))



    #comp_len = sh.wc(sh.gzip(f=True, _in=str(x.seq)), c=True)

    #sh.run_abyss(r1, r2, 'name=%s' % prefix, 'k=%s' % 25, _err=log, _out=log)
#    args = ["in=%s %s" % (r1, r2), 'name=%s' % prefix, 'k=%s' % 25]
#    subprocess.call('abyss-pe ' + ' '.join(args))
    #subprocess.call('abyss-pe', ' '.join(args))
#    print ['abyss-pe'] + args
    #subprocess.call(['abyss-pe'] + args) #, stdout=log, stderr=log, shell=True)
#    ex = "abyss-pe in=%s %s" % (r1, r2) +  ' name=%s ' % prefix + ' k=%s' % 25
#    subprocess.call(ex, stdout=log, stderr=log, shell=True)

#    sh.abyss_pe('in=\'%s %s\'' % (r1, r2), 'name=%s' % prefix, 'k=%s' % 25, # '-n'  dryrun
#                _err=log, _out=log)
#    sh.abyss_pe("in='%s %s'" % (r1, r2), name=prefix, k=25, # '-n'  dryrun
#                _err=log, _out=log, _long_prefix='', _short_prefix='')

#    abyss-pe k=25 name=test     in='test-data/reads1.fastq test-data/reads2.fastq'



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

  _bowtie1 =   p( "bowtie.1.r1" )
  _bowtie2 =   p( "bowtie.2.r1" )

  contigs   = p("abyss-contigs.fa")

  nr = p('contigs.nr.tsv')
  nt = p('contigs.nt.tsv')
  kronaNT = p('contigs.nt.html')
  kronaNR = p('contigs.nr.html')

#  bowtie1 =   p( "bowtie.r1.fa" )
#  bowtie2 =   p( "bowtie.r2.fa" )
#  nr1     =   p( "rapsearch.r1.blast.m8" ) # rapsearch automatically adds .m8 extension
#  nr2     =   p( "rapsearch.r2.blast.m8" ) # rapsearch automatically adds .m8 extension
#
#  nt1 =       p( "r1.blast" )
#  nt2 =       p( "r2.blast" )
#  kronaNT1  = p( "r1.NT.html" )
#  kronaNT2  = p( "r2.NT.html" )
#  kronaNR1  = p( "r1.NR.html" )
#  kronaNR2  = p( "r2.NR.html" )

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

  if need(_bowtie1):
    bowtie_sensitive(log, cfg, cd1, cd2, _bowtie1)

  if need(contigs):
    abyss(log, cfg, _bowtie1, _bowtie2, contigs)
  if need(nt):
    blastn(log, cfg, contigs, nt)
  if need(nr):
    blastx(log, cfg, contigs, nr)
  if need(kronaNT):
    krona(log, cfg, nt, kronaNT)
  if need(kronaNR):
    krona(log, cfg, nr, kronaNR)

#  if need(bowtie1):
#    SeqIO.convert(_bowtie1, 'fastq', bowtie1, 'fasta')
#    SeqIO.convert(_bowtie2, 'fastq', bowtie2, 'fasta')
#
#  if need(nt1):
#    blastn(log, cfg, bowtie1, nt1)
#    blastn(log, cfg, bowtie2, nt2)
#
#  if need(nr1):
#    #rapsearch(log, cfg, bowtie1, nr1)
#    #rapsearch(log, cfg, bowtie2, nr2)
#    blastx(log, cfg, bowtie1, nr1)
#    blastx(log, cfg, bowtie2, nr2)
#
#  if need(kronaNT1):
#    krona(log, cfg, nt1, kronaNT1)
#    krona(log, cfg, nt2, kronaNT2)
#    krona(log, cfg, nr1, kronaNR1)
#    krona(log, cfg, nr2, kronaNR2)

def main():
  args = docopt(__doc__, version='Version 1.0')
  cfg = args['--config']
  cfg = yaml.load(open(cfg))
  cfg = Config(cfg)
  cfg.outdir = args['-o'] or "."
  if args['--log']:
    _log = Path(args['--log'])
    if _log.exists():
      print "Removing old log file %s" % _log
      _log.remove()
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
