module Pipeline where
import Data.Char(ord,chr,toUpper)
import qualified Data.List as L
import Data.Word (Word8)
import System.Directory(doesFileExist)
import Control.Applicative
import Bio.Sequence.FastQ (readSangerQ, writeSangerQ)
import Bio.Core.Sequence (BioSeqQual, BioSeq, unSD, seqheader, seqdata, seqqual, unQD, unQual)
import qualified Data.ByteString.Lazy.Char8 as C
import qualified Data.ByteString.Lazy as BB
import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util
import qualified Data.Compressed.LZ78 as LZW
import Config
import qualified Data.Yaml as Y

data ReadPair = RP FilePath FilePath
type Fastq = FilePath
type Sam = FilePath
data SeqType = NT | NR

--block :: Config -> Fastq -> Fastq -> Action ()
block config input1 input2 = shakeArgs shakeOptions{shakeFiles="_build"} $ do

  let ignore = "bar"
   -- paste 'cat's files "vertically"
  let paste src1 src2 out = cmd Shell "paste -d \t" src1 src2 ">" out

  want ["rapsearch.tsv", "gsnap.tsv"] -- final result of build

  "rapsearch.tsv" %> \out -> do
    taxa <- needExt "tax" out
    m8   <- needExt "m8" out
    paste m8 taxa out

  "gsnap.tsv" %> \out -> do
    taxa <- needExt "tax" out
    sam   <- needExt "sam" out
    paste sam taxa out

  "rapsearch.tax" %> \out -> do
    acc2tax NR out =<< needExt "gi" out

  "gsnap.tax" %> \out -> do
    acc2tax NT out =<< needExt "gi" out

  "rapsearch.gi" %> \out -> do
    cmd Shell "cut -f 2 | cut -d \\| -f 2 >" out-- need to extract the GI
--    let rows' = dropWhile (\xs -> (head xs) == '@') rows
--    let gids  = map (\x -> head $ splitOn "|" $ head $ splitOn "\t" x) rows'
--    writeFile' out $ L.intersperse "\n" gids

  "gsnap.gi" %> \out -> do
    cmd "gibber!"
--  see above

  "rapsearch.m8" %> \out -> do
    (r1, r2) <- need2 "r1.bowtie" "r2.bowtie"
    cmd "rapsearch" "-e" r1 r2 "-o" out "-d" $ rapsearchDB $ rapsearch config

  "gsnap.samsv" %> \out -> do
    sam <- needExt "sam" out
    writeFile' out "query\tsubject\tcigar"
    cmd Shell "grep -v ^@" sam "| awk -v OFS='\t' '{print $1, $3 $6}' >" out

  "gsnap.sam" %> \out -> do
    (r1, r2) <- need2 "r1.bowtie" "r2.bowtie"
    cmd "gsnapl"r1 r2 "-A" out "-d" $ gsnapDB $ gsnap config 

  ["r1.bowtie", "r2.bowtie"] &%> \[o1, o2] -> do
    sam2fastq (RP o1 o2) =<< need' "bowtie.sam"

  "bowtie.sam" %> \out -> do
    (r1, r2) <- need2 "R1.deduped" "R2.deduped"
    cmd "bowtie2" "-1" r1 "-2" r2 "-S" out "--very-sensitive-local" "-x" $ bowtieDB $ bowtie2 config

  ["R1.deduped", "R2.deduped"] &%> \[out1, out2] -> do
    (r1, r2) <- needExt2 "lzw" out1 out2
    runCdHitDup (cdhitdup config) (RP r1 r2) (RP out1 out2)

  ["R1.lzw", "R2.lzw"] &%> \[out1, out2] -> do
    (src1, src2) <- needExt2 "deduped" input1 input2
    liftIO $ filterPair (not . (lowComplex 0.45)) src1 src2 out1 out2

  ["r1.psf", "r2.psf"] &%> \[o1, o2] -> do
    (src1, src2) <- needExt2 "star" o1 o2
    runPriceFilter (pricefilter config) (RP src1 src2) (RP o1 o2)

  ["r1.star", "r2.star"]  &%> \[r1, r2] -> do
    let src = need' "star.sam"
    sam2fastq (RP r1 r2) =<< need' "star.sam"

  "star.sam" %> \out -> do
    need [input1, input2]
    runSTAR (star config) (RP input1 input2) out

  where
    sam2fastq :: ReadPair -> Sam -> Action () 
    sam2fastq (RP r1 r2) sam = cmd "java -jar picard.jar SamToFASTQ" 
      ("I" <=> sam) ("FASTQ" <=> r1) ("SECOND_END_FASTQ" <=> r2)
      where x <=> y = x ++ "=" ++ y

    --(input1, input2) = ("R1.fq", "R2.fq") 
    runBowtie (BowtieOpts {bowtieDB=db}) (RP r1 r2) out = cmd' where
      cmd' = cmd "bowtie2" "-1" r1 "-2" r2 "-x" db "-S" out "--very-sensitive-local"

    runCdHitDup :: CDHitOpts -> ReadPair -> ReadPair -> Action ()
    runCdHitDup CDHitOpts {minDifference=mindiff} (RP r1 r2) (RP o1 o2) = cmd' where
      cmd' = cmd "cd-hit-dup"
        "-i" r1  "-i2" r2
        "-o" o1 "-o2" o2
        "-e" $ show mindiff

    runPriceFilter :: PriceSeqFilterOpts -> ReadPair -> ReadPair -> Action ()
    runPriceFilter PriceSeqFilterOpts {calledPercent=cp, highQualPercent=hqp, highQualMin=hqm} (RP r1 r2) (RP o1 o2) = cmd' where
      cmd' = cmd "PriceSeqFilter"
        "-fp" r1 r2
        "-op" o1 o2
        "-rnf" (show cp)
        "-rqf" (show hqp) (show hqm)
    runSTAR :: STAROpts -> ReadPair -> FilePath -> Action ()
    runSTAR STAROpts{starDB=db'} (RP r1 r2) out = cmd'
      where
        cmd' = cmd "--readFilesIn" r1 r2 "--genomeDir" db' "--runThreadN" (show $ threads config) "--outFileNamePrefix" (out -<.> "") "--outSAMtype" "SAM"

    acc2tax :: SeqType -> FilePath -> FilePath -> Action ()
    acc2tax st fpout fpin = cmd "acc2tax" "--gi" "-i" fpin "-o" fpout "--database" db "--entries" ents type'
      where
        ents = show 1114054124
        (db, type') = case st of -- also switch to get NR/NT db
          NT -> (ntDB $ ncbi config, "--nucleotide")
          NR -> (nrDB $ ncbi config, "--protein")

main = printConfig

need' x = need [x] >> return x
needExt  ext fp = need' (fp -<.> ext)

needExt2 :: String -> FilePath -> FilePath -> Action (FilePath, FilePath)
needExt2 ext a b = need2 (a -<.> ext) (b -<.> ext)

need2 :: FilePath -> FilePath -> Action (FilePath, FilePath)
need2        a b = need [a, b] >> return (a, b)

lowComplex n s = compScore < n  where
  compScore    = (ucSize - cSize) / ucSize
  cSize        = length' $ LZW.encode uncompressed -- encode expects [a]
  ucSize       = length' uncompressed
  uncompressed = BB.unpack $ unSD $ seqdata s -- so need to unpack
  length' xs = foldr (\_ z -> z + 1) 0 xs

filterPair f in1 in2 o1 o2 = do
  r1 <- readSangerQ in1
  r2 <- readSangerQ in2
  let (r1', r2') = unzip $ filter pred $ zip r1 r2
  writeSangerQ o1 r1'
  writeSangerQ o2 r2'
  where
    pred (fwd, rev) = (f fwd) && (f rev)
