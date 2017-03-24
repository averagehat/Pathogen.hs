
'''
Usage:
      pipeline.py <r1> <r2> --config <config> [-o <outdir>]
'''
import sh
import itertools
from toolz.dicttoolz import keymap,valfilter
from docopt import docopt
from path import Path
import types
# TODO: LZW!


#############
# Utilities #
#############
class Config:
    def __init__(self, **entries):
        self.__dict__.update(entries)
        for k,v in self.__dict__.items():
            if type(v) == types.DictType:
                setattr(self, k, Struct(**v))
class Sh_(object):
    def __getattr__(self, attr):
        #cmd = getattr(sh, attr)
        def command(*args, **kwargs):
            #fixedargs = keymap("-{}".format, kwargs)
            bools = valfilter(lambda x: type(x) is bool, kwargs)
            vargs = keymap("-{}".format, valfilter(lambda x: not type(x) is bool, kwargs))
            #bools.update(vargs)
            print bools
            fixedargs = itertools.chain(vargs.items())
            getattr(sh, attr)(*(list(args) + list(fixedargs)), **bools)
        return command
sh_ = Sh_()
def unlist(*xs):
  return " ".join(xs)

############
# Parts    #
############

def star(cfg, in1, in2):
  sh.STAR(readFilesIn=unlist(in1, in2),
          genomeDir=cfg.star.starDB,
          outSamype="SAM",
          outReadsUnmapped="Fastx")

def pricefilter(cfg, in1, in2, o1, o2):
    cfg = cfg['pricefilter']
    sh_.PriceSeqFilter(fp=unlist(in1, in2),
                       op=unlist(o1, o2),
                       rnf=cfg.calledPercent,
                       rqf=unlist(cfg.highQualPercent, cfg.highQualMin))

def cdhitdup(cfg, r1, r2, o1, o2):
    sh_.cd_hit_dup(i=r1, i2=r2, o=o1, o2=o2, e=cfg.cdhitdup.minDifference)

# LZW!

def bowtie_sensitive(cfg, r1, r2, o1, log):
    args = {'1' : r1, '2' : r2,
                'very_sensitive_local' : True,
                'un_conc' : Path(o1).splitext()[0],
                'x' : cfg.bowtie2.bowtieDB}
    sh.bowtie2(**args)

def rapsearch(cfg, fq, out):
    sh.rapsearch(o=out, d=cfg.rapsearch.rapsearchDB, q=fq)

def blastn(cfg, fq, out):
    sh_.blastn(outfmt=8, db=cfg.ncbi.ntDB, o=out, q=fq)

def krona(cfg, blast, out):
    sh.ktImportBlast(blast, o=out) # probably need config for kronadb!


############
# Pipeline #
############


star1 = "star.r1.fq"
star2 = "star.r2.fq"

psf1 = "psf.r1.fq"
psf2 = "psf.r2.fq"

cd1 = "cd.r1.fq"
cd2 = "cd.r2.fq"

bowtie1 = "bowtie.r1.fq"
bowtie2 = "bowtie.r2.fq"

nr1     = "rapsearch.r1.blast"
nr2     = "rapsearch.r2.blast"

nt1 = "r1.blast"
nt2 = "r2.blast"

kronaNT1  = "r1.NT.html"
kronaNT2  = "r2.NT.html"
kronaNR1  = "r1.NR.html"
kronaNR2  = "r2.NR.html"
def run(cfg, input1, input2):
  star(cfg, input1, input2)

  pricefilter(cfg, star1, star2, psf1, psf2)

  cdhitdup(cfg, psf1, psf2, cd1, cd2)

  bowtie_sensitive(cfg, cd1, cd2, bowtie1)

  blastn(cfg, bowtie1, nt1)
  blastn(cfg, bowtie2, nt2)

  rapsearch(cfg, bowtie1, nr1)
  rapsearch(cfg, bowtie2, nr2)

  krona(cfg, nt1, kronaNT1)
  krona(cfg, nt2, kronaNT2)
  krona(cfg, nt1, kronaNR1)
  krona(cfg, nt2, kronaNR2)

def main():
  args = docopt(__doc__, version='Version 1.0')
  cfg = args['<config>']
  cfg = yaml.load(open('template.yaml'))
  cfg = Config(cfg)
  run(cfg, args['<r1>'], args['<r2>'])
  sys.exit(0)

if __name__ == '__main__': main()
