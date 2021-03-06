module Pipeline where
import Data.Char(ord,chr,toUpper)
import qualified Data.List as L
import Data.Word (Word8)
import System.Directory(doesFileExist)
import Control.Applicative
import Bio.Sequence.FastQ (readSangerQ, writeSangerQ)
import Bio.Sequence.Fasta (writeFasta, Sequence(..)) --, mkSeq, bleh)
import Bio.Core.Sequence (BioSeqQual, BioSeq, unSD, seqheader, seqdata, seqqual, unQD, unQual, seqid)
import qualified Data.ByteString.Lazy.Char8 as C
import qualified Data.ByteString.Lazy as BB
import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util
import Data.Traversable (for)
import qualified Data.Compressed.LZ78 as LZW
import Config
import qualified Data.Yaml as Y
import System.Time.Extra (offsetTime, showDuration)
import System.IO.Unsafe (unsafePerformIO)
import Options.Generic  (getRecord, unHelpful)
import qualified Data.Text as T
import Index (buildIndex)
data ReadPair = RP FilePath FilePath
type Fastq = FilePath
type Sam = FilePath
data SeqType = NT | NR


-- TODO: get some toy databases (not hg38)

--block :: Config -> Fastq -> Fastq -> Action ()
-- replace "_build" with outdir
--block config input1 input2 = shakeArgs shakeOptions{shakeFiles="_build"} $ do
block config input1 input2 = shake shakeOptions $ do
  --(input1, input2) = ("R1.fq", "R2.fq") 
  let r1_deduped = input1 -<.> "deduped"
  let r2_deduped = input2 -<.> "deduped"

  let ignore = "bar"

  --want ["rapsearch.tsv", "gsnap.tsv"] -- final result of build
  want ["NR.tax.R1", "NR.tax.R2", "NT.tax.R1", "NT.tax.R2"]

  let r1_star = "Unmapped.out.mate1.fq"
  let r2_star = "Unmapped.out.mate2.fq"

  let r1_bowtie = "bowtie.1.unmapped"
  let r2_bowtie = "bowtie.2.unmapped"

  let r1_psf  = "R1.psf.fq"
  let r2_psf  = "R2.psf.fq"

  pairRule "R1" r1_bowtie

  pairRule "R2" r2_bowtie

  [r1_bowtie, r2_bowtie, "bowtie.mapped"] &%> \[o1, _, log] -> do
    (r1, r2) <- need2 "R1.deduped" "R2.deduped"
    cmd Shell "bowtie2" "-1" r1 "-2" r2 "--very-sensitive-local" 
      "--un-conc" (dropExtension o1)
      "-x" (bowtieDB $ bowtie2 config)
      ">" log

  ["R1.lzw", "R2.lzw"] &%> \[out1, out2] -> do
    (src1, src2) <- needExt2 "deduped" out1 out2
    liftIO $ filterPair (not . lowComplex 0.45) src1 src2 out1 out2 

  ["R1.deduped", "R2.deduped"] &%> \[out1, out2] -> do
    --(r1, r2) <- needExt2 "psf" out1 out2
    (r1, r2) <- need2 r1_psf r2_psf
    runCdHitDup (cdhitdup config) (RP r1 r2) (RP out1 out2)

  [r1_psf, r2_psf] &%> \[o1, o2] -> do
    need [r1_star, r2_star]
    runPriceFilter (pricefilter config) (RP r1_star r2_star) (RP o1 o2)

  [r1_star, r2_star] &%> \[o1, _] -> do
    need [input1, input2]
    () <- cmd "STAR" "--readFilesIn" input1 input2 
        "--genomeDir" (starDB $ star config) 
        "--runThreadN" (show $ threads config) 
        "--outSAMtype SAM" "--outReadsUnmapped Fastx" 
    -- need this business because Price only accepts files with .fq/.fastq extension.
    () <- cmd "mv" (dropExtension r1_star) r1_star
    cmd "mv" (dropExtension r2_star) r2_star

  where
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
        cmd' = cmd "STAR" "--readFilesIn" r1 r2 "--genomeDir" db' "--runThreadN" (show $ threads config) "--outFileNamePrefix" (out -<.> "") "--outSAMtype SAM"

    acc2tax :: SeqType -> FilePath -> FilePath -> Action ()
    acc2tax st fpout fpin = cmd "acc2tax" "--gi" "-i" fpin "-o" fpout "--database" db "--entries" ents type'
      where
        ents = show 1114054124
        (db, type') = case st of
          NT -> (ntDB $ ncbi config, "--nucleotide")
          NR -> (nrDB $ ncbi config, "--protein")

    pairRule :: String -> FilePath -> Rules ()
    pairRule s fq = do

     for ["NT", "NR"] $ \t -> do

       t <.> "tsv" <.> s %> \out -> do
         taxa <- need' $ t <.> "tax" <.> s
         m8   <- need' $ t <.> "m8" <.> s
         paste m8 taxa out

       t <.> "tax" <.> s %> \out ->
         acc2tax NR out =<< need' (t <.> "gi" <.> s)

       t <.> "gi" <.> s %> \out -> 
          writeFile' out . getGI . T.pack =<< readFile' =<< need' (t <.> "m8" <.> s)

     "NR.m8" <.> s %> \out -> 
       cmd "rapsearch -o" out "-d" (rapsearchDB $ rapsearch config) "-q" =<< need' fq
  
     "NT.m8" <.> s %> \out ->
       cmd "blastn -outfmt 8 -db" (ntDB $ ncbi config) "-o" out "-q" =<< need' fq 
      -- paste 'cat's files "vertically"
      where paste src1 src2 out = cmd Shell "paste -d \t" src1 src2 ">" out

need' x = need [x] >> return x
needExt  ext fp = need' (fp -<.> ext)

needExt2 :: String -> FilePath -> FilePath -> Action (FilePath, FilePath)
needExt2 ext a b = need2 (a -<.> ext) (b -<.> ext)

need2 :: FilePath -> FilePath -> Action (FilePath, FilePath)
need2        a b = need [a, b] >> return (a, b)

getGI :: T.Text -> String
getGI s = let rows = dropWhile (\xs -> T.head xs == '#') $ T.lines s
              gis  = map (\x -> flip (!!) 1 $ T.splitOn (T.pack "|") $ T.splitOn (T.pack "\t") x !! 1) rows  in
    L.intercalate "\n" $ map T.unpack gis

lowComplex n s = compScore < n  where
  compScore    = (ucSize - cSize) / ucSize
  cSize        = length' $ LZW.encode uncompressed -- encode expects [a]
  ucSize       = length' uncompressed
  uncompressed = BB.unpack $ unSD $ seqdata s -- so need to unpack
  length' xs = foldr (\_ z -> z + 1) 0 xs

fastq2fasta fq fa = do
 fq' <- readSangerQ fq
 let fas = map (\x -> Seq (seqid x) (seqdata x) Nothing) fq'
 writeFasta fa fas

filterPair f in1 in2 o1 o2 = do
  r1 <- readSangerQ in1
  r2 <- readSangerQ in2
  let (r1', r2') = unzip $ filter pred $ zip r1 r2
  writeSangerQ o1 r1'
  writeSangerQ o2 r2'
  where
    pred (fwd, rev) = f fwd && f rev

-- https://github.com/ndmitchell/shake/issues/483

{-# NOINLINE logTime #-}
--logTime :: IO Seconds
logTime = unsafePerformIO offsetTime

time :: IO (CmdLine, CmdTime, a) -> IO a
time act = do
    (CmdLine msg, CmdTime tim, res) <- act
    tot <- logTime
    putStrLn $ "[BAKE-TIME] " ++ showDuration tim ++ " (total of " ++ showDuration tot ++ "): " ++ msg
    return res

time_ :: IO (CmdLine, CmdTime) -> IO ()
time_ act = time $ do (a,b) <- act; return (a,b,()) 

main = do
  opts <- getRecord $ T.pack "Running Pathogen.hs" :: IO CommandArgs
  case opts of 
    Index{} -> buildIndex opts
    Run{}   -> runPipeline opts

runPipeline opts = do 
  cfg <- readYaml $ unHelpful $ config opts
  case cfg of
    Right cfg' -> block cfg' (unHelpful $ r1 opts) (unHelpful $ r2 opts) 
    Left err -> putStrLn err 


--  "NT.gi" %> \out -> do
--    writeFile' out . getGI . T.pack =<< readFile' =<< need' "NT.tsv"
--  
--  "NR-R1.m8" %> \out ->
--    cmd "rapsearch -o" out "-d" (rapsearchDB $ rapsearch config) "-q" =<< need' r1_bowtie
--
--  "NR-R2.m8" %> \out ->
--    cmd "rapsearch -o" out "-d" (rapsearchDB $ rapsearch config) "-q" =<< need' r2_bowtie
--
--  "NT-R1.tsv" %> \out -> do --- should make e-value etc. configurable
--     cmd "blastn -m 8 -db" (ntDB $ ncbi config) "-o" out "-q" =<< need' r1_bowtie
--    
--  "NT-R2.tsv" %> \out -> do --- should make e-value etc. configurable
--     cmd "blastn -m 8 -db" (ntDB $ ncbi config) "-o" out "-q" =<< need' r2_bowtie

--  "gsnap.gi" %> \out -> do
--    src <- need' "gsnap.samsv"
--    cmd "gibber!"
----  see above
----
--  "gsnap.samsv" %> \out -> do
--    sam <- needExt "sam" out
--    writeFile' out "query\tsubject\tcigar"
--    cmd Shell "grep -v ^@" sam "| awk -v OFS='\t' '{print $1, $3 $6}' >" out
--
--  "gsnap.sam" %> \out -> do
--    (r1, r2) <- need2 r1_bowtie r2_bowtie
--    cmd "gsnapl"r1 r2 "-A" out "-d" $ gsnapDB $ gsnap config 

--  "rapsearch.tsv" %> \out -> do
--    taxa <- needExt "tax" out
--    m8   <- needExt "m8" out
--    paste m8 taxa out
--
--  "gsnap.tsv" %> \out -> do
--    taxa <- needExt "tax" out
--    sam   <- needExt "sam" out
--    paste sam taxa out
--
--  "rapsearch.tax" %> \out -> do
--    acc2tax NR out =<< needExt "gi" out
--
--  "gsnap.tax" %> \out -> do
--    acc2tax NT out =<< needExt "gi" out
--
--  "rapsearch.gi" %> \out -> do
--    src <- need' "rapsearch.m8"
--    --cmd Shell "cut -f 2 | cut -d \\| -f 2 >" out-- need to extract the GI
--    s <- T.pack <$> readFile' src
--    writeFile' out $ getGI s
--
--  "rapsearch.m8" %> \out -> do
--    (r1, r2) <- need2 "r1.bowtie" "r2.bowtie"
--    cmd "rapsearch" "-e" r1 r2 "-o" out "-d" $ rapsearchDB $ rapsearch config
--  
