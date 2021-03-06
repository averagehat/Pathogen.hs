'''
Usage:
    pipeline.py <r1> <r2> --config <config> [-o <outdir>] [--log <log>]

Options:
    -c <config>, --config <config>
    --log <log>
'''
import sh
from itertools import imap, izip, tee, chain, groupby, ifilter, starmap
from toolz.dicttoolz import keymap,valfilter,keyfilter,merge
from toolz.itertoolz import mapcat
from docopt import docopt
from path import Path
import types
import yaml
import sys
from functools import partial
import shutil
import os
from Bio import SeqIO
from operator import itemgetter as get
from collections import Counter
import csv
from ete2 import NCBITaxa
import plumbum
#import numpy as np
#import matplotlib.pyplot as plt
# TODO: Log commands as doing them
# TODO: BLAST Contig results with mapped reads to get abundance:
#  - Does abyss record number of reads into contig? # no, it's just length and "kmer coverage"
#  - could duplicate contig blast entry for each read that maps to it and pass to krona
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
            fixedargs = chain(vargs.items())
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
    sh.blastn('-max_target_seqs', '1', outfmt=6, db=cfg.ncbi.ntDB, query=fq, _err=log, _out=out, _long_prefix='-')

def krona(log, cfg, blast, out):
    sh.ktImportBLAST(blast, o=out, _err=log, _out=log) # probably need config for kronadb!


def blastx(log, cfg, fq, out):
    sh.blastx('-max_target_seqs', '1', outfmt=6, db=cfg.ncbi.nrDB, query=fq, _err=log, _out=out, _long_prefix='-')

def abyss(log, cfg, r1, r2, out):
    dir = out.dirname()
    f1 = dir.relpathto(r1)
    f2 = dir.relpathto(r2)
    prefix=out.basename().split('-')[0]
    print f1, f2, 'name=%s' % prefix, 'k=%s' % 25
    sh.run_abyss(f1, f2, 'name=%s' % prefix, 'k=%s' % 25, C=dir, _err=log, _out=log)


#def lzw_filter(log, cfg, r1, r2, out1, out2):
def lzw_filter_single(min_complexity, x):
    un_comp_len = len(str(x.seq))
    comp_len = sum(imap(len, sh.gzip(f=True, _in=str(x.seq))))
    complexity =  comp_len / float(un_comp_len)
    return complexity >= min_complexity

def unzip(seq):
    t1, t2 = tee(seq)
    return imap(get(0), t1), imap(get(1), t2)

def filter_pair(func, r1, r2, o1, o2, format):
    fwd = SeqIO.parse(r1, format)
    rev = SeqIO.parse(r2, format)
    filtered = ((x, y) for (x, y) in izip(fwd, rev)
                if func(x) and func(y))
    res1, res2 = unzip(filtered)
    with open(o1, 'w') as f1:
        with open(o2, 'w') as f2:
            SeqIO.write(res1, f1, format)
            SeqIO.write(res2, f2, format)

def lzw_filter_fastq(log, cfg, r1, r2, out1, out2):
    lzw_func = partial(lzw_filter_single, cfg.lzwfilter.maxCompressionScore)
    filter_pair(lzw_func, r1, r2, out1, out2, 'fastq')

def sum_sam_by_ref(log, cfg, sam):
    res = sh.samtools.view(sam, F=260)
    refs = imap(lambda x: x.split('\t')[2], res)
    return Counter(refs)

def dup_blast(log, sam, blst, out):
    counter = sum_sam_by_ref(None, None, sam)
    log.write("Skipped Contigs:\n======\n")
    with open(out, 'w') as f:
        with open(blst, 'r') as blast:
            for k,v in groupby(blast, lambda x: x.split('\t')[0]):
                if k not in counter:
                    contig_length = next(v).split('\t')[3]
                    log.write("Contig %s of length %s had no mapped reads.\n" % (k, contig_length))
                else:
                    f.writelines(list(v) * counter[k])
    log.write("======\n")



#                # because SeqIO is slow for writing single records
#        for (x, y) in izip(fwd, rev):
#            if func(x) and func(y):
#                '@
#                f1.write(x.format(form))
#                f2.write(y.format(form))



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


# taxid = 1056490

def blastdbcmd(**opts):
    cmd_opts = keymap('-{}'.format, opts).items()
    process = plumbum.local['blastdbcmd'][cmd_opts]
    print process
    for line in process.popen().iter_lines(retcode=None):
        yield line[0]

def get_taxid(db, seqids): # (Path, str) -> dict[str,str]
   #res = sh.blastdbcmd(db=db, outfmt="'%g %T'", entry="'%s'" % seqid, _long_prefix='-')
   max_ids = 1000/80
   if len(seqids) > max_ids:
      xs = [seqids[i*max_ids:(i+1)*max_ids] \
         for i in range(len(seqids) / max_ids)]
      xs.extend(seqids[sum(map(len, xs))-1:])
   else: xs = [seqids]
#   xs = chain(*xs)
#   print map(len, xs)
   print seqids
   print xs
   res = mapcat(lambda x: blastdbcmd(db=db, outfmt="'%g %T'", entry="'%s'" % ','.join(x)), xs)
   #res = blastdbcmd(db=db, outfmt="'%g %T'", entry="'%s'" % ','.join(xs))
   #print res
   #print res
   res = ifilter(bool, res)
   res = imap(lambda s: s.strip("'"), res)
   return dict(imap(unicode.split, res))

def taxonomy(ncbi, taxid):
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    def make_d(lineage, ranks, names):
        for lin in lineage:
            if ranks[lin] == 'no rank':
                continue
            yield (ranks[lin], names[lin])
    return dict(make_d(lineage, ranks, names))
def dictmap(f, d): return starmap(f, d.items())
def blast2summary_dict(db, blastpath): # (Path, Path) -> list[dict]

  """Reading in a blast output file, lookup all seqids to get taxids with a single blastdbcmd.
  Then, lookup the taxonomy using ETE2 via the taxid, and add that info to the blast info."""
  rows = csv.DictReader(open(blastpath), delimiter='\t',fieldnames=['qseqid', 'sseqid','pid', 'alnlen','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
  rows = list(rows)
  seqids = map(get('sseqid'), rows)
  taxids = get_taxid(db, seqids)
  gis = (s.split('|')[1] for s in seqids)
  matches = dict((taxids[gi], row) for gi, row in zip(gis,rows) if gi in taxids)
  ncbi = NCBITaxa() # downloads database and creates SQLite database if needed
  return dictmap(lambda tid,row: merge(row, taxonomy(ncbi, tid)), matches)

def blast2summary(db, blastpath, outpath): # (Path,Path,Path) -> None
    with_taxonomies = list(blast2summary_dict(db, blastpath))
    head = with_taxonomies[0]
    print with_taxonomies
    with open(outpath, 'w') as out:
       writer = csv.DictWriter(out, head.keys(), delimiter='\t')
       #writer = csv.DictWriter(out, fieldnames=head.keys(), delimiter='\t')
       writer.writeheader()
       for row in with_taxonomies:
         writer.writerow(row)

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

  lzw1 = "lzw.r1"
  lzw2 = "lzw.r2"

  _bowtie1 =   p( "bowtie.1.r1" )
  _bowtie2 =   p( "bowtie.2.r1" )

  contigs   = p("abyss-contigs.fa")
  contigs_sam = 'contigs.sam'

  contig_nr = p('contigs.nr.blast')
  contig_nt = p('contigs.nt.blast')

  dup_nt = p('contigs.nt.blast.dup')
  dup_nr = p('conrigs.nr.blast.dup')
  contig_kronaNT = p('contigs.nt.html')
  contig_kronaNR = p('contigs.nr.html')

  contig_nt_tsv = p("contigs.nt.tsv")
  conrig_nr_tsv = p("conrigs.nr.tsv")
  nt_tsv = p('nt.tsv')
  nr_tsv = p('nr.tsv')
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

  if need(lzw1):
    lzw_filter_fastq(log, cfg, cd1, cd2, lzw1, lzw2)

  if need(_bowtie1):
    bowtie_sensitive(log, cfg, lzw1, lzw2, _bowtie1)

  if need(contigs):
    abyss(log, cfg, _bowtie1, _bowtie2, contigs)
    contigs_index = 'contigs-b2'
    sh.bowtie2_build(contigs, contigs_index)
    sh.bowtie2(**{'1' : _bowtie1, '2' : _bowtie2, 'x' : contigs_index,
                  '_err' : log, '_out' : contigs_sam})

  if need(contig_nt):
    blastn(log, cfg, contigs, contig_nt)
    # TODO: fix below
    dup_blast(log, contigs_sam, contig_nt, dup_nt)
  if need(contig_nr):
    blastx(log, cfg, contigs, contig_nr)
    dup_blast(log, contigs_sam, contig_nr, dup_nr)
  if need(contig_kronaNT):
    krona(log, cfg, contig_nt, contig_kronaNT)
  if need(contig_kronaNR):
    krona(log, cfg, contig_nr, contig_kronaNR)
  if need(contig_nt_tsv):
    blast2summary(cfg.ncbi.ntDB, contig_nt, contig_nt_tsv)

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
  run(cfg, args['<r1>'], args['<r2>'], log)
#  try:
#    run(cfg, args['<r1>'], args['<r2>'], log)
#  except Exception as e:
#    log.write(str(e))
  if args['--log']: log.close()
  sys.exit(0)

if __name__ == '__main__': main()



