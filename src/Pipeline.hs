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

block cfg = shakeArgs shakeOptions{shakeFiles="_build"} $ do

  want ["rapsearch.tsv", "gsnap.tsv"] -- final result of build
  let ignore = "bar"
  ["R1.lzw", "R2.lzw"] &%> \[out1, out2] -> do
    let [src1, src2] = [input1 -<.> "deduped", input2 -<.> "deduped"]
    need [src1, src2]
    liftIO $ filterPair (not . (lowComplex 0.45)) src1 src2 out1 out2

  "rapsearch.gi" %> \out -> do
    src  <- need' "rapsearch.m8"
    rows <- lines <$> readFile' src
    liftIO $ putStrLn "foo"
--    let rows' = dropWhile (\xs -> (head xs) == '@') rows
--    let gids  = map (\x -> head $ splitOn "|" $ head $ splitOn "\t" x) rows'
--    writeFile' out $ L.intersperse "\n" gids

  "rapsearch.tsv" %> \out -> do
    src <- need' $ out -<.> "gi"
    acc2tax NT src out

  ["r1.star", "r2.star"]  &%> \[r1, r2] -> do
    let src = need' "star.bam"
    -- somehow extract fwd and reverse out of the bam file
    liftIO $ putStrLn "hey"

  "star.bam" %> \out -> do
    need [input1, input2]
    runSTAR (STAROpts {starDB = "/baz"}) (RP input1 input2) out

  where
    (input1, input2) = ("R1.fq", "R2.fq")
    threads = 16
    runSTAR :: STAROpts -> ReadPair -> FilePath -> Action ()
    runSTAR STAROpts{starDB=db'} (RP r1 r2) out = cmd'
      where
        cmd' = cmd "--readFilesIn" r1 r2 "--genomeDir" db' "--runThreadN" (show threads) "--outFileNamePrefix" (out -<.> "") "--outSAMtype" "BAM"

    acc2tax :: SeqType -> FilePath -> FilePath -> Action ()
    acc2tax st fpin fpout = cmd "acc2tax" "--gi" "-i" fpin "-o" fpout "--database" db "--entries" ents type'
      where
        ents = show 1114054124
        (db, type') = case st of -- also switch to get NR/NT db
          NT -> ("NT", "--nucleotide")
          NR -> ("NR", "--protein")

main = printConfig
--  cfg <- (Y.decodeEither <$> readFile "config.yaml" :: Maybe Config)
--  either block error' $ (bowtieDB . bowtie2) <$> cfg where
--    error' = putStrLn

need' x = do
  need [x]
  return x

data SeqType = NT | NR



--lowComplex :: BioSeqQual s => Double -> s -> Bool
--lowComplex :: Double -> Sequence -> Bool
--lowComplex :: BioSeq s => Double -> s -> Bool
lowComplex n s = compScore < n  where
  compScore    = (ucSize - cSize) / ucSize
  cSize        = length' $ LZW.encode uncompressed -- encode expects [a]
  ucSize       = length' uncompressed
  uncompressed = BB.unpack $ unSD $ seqdata s -- so need to unpack
  --length' :: Foldable f => (f a) -> Int
  length' xs = foldr (\_ z -> z + 1) 0 xs

--filterPair :: BioSeq s => (s -> Bool) -> FilePath -> FilePath -> FilePath -> FilePath -> IO ()
-- filterPair :: (Sequence -> Bool) -> FilePath -> FilePath -> FilePath -> FilePath -> IO ()
filterPair f in1 in2 o1 o2 = do
  r1 <- readSangerQ in1
  r2 <- readSangerQ in2
  let (r1', r2') = unzip $ filter pred $ zip r1 r2
  writeSangerQ o1 r1'
  writeSangerQ o2 r2'
  where
    pred (fwd, rev) = (f fwd) && (f rev)


-- | Split a string (note that | is needed for doctest
-- >>> splitOn 'x' "AAxAA"
-- ["AA", "AA"]
splitOn delimiter = foldr f [[]]
            where f c l@(x:xs) | c == delimiter = []:l
                                | otherwise = (c:x):xs


-- cmd' :: Show a => [(String, a)] -> IO ()
-- cmd' opts = cmd $ concat $  map (\(opt,val) -> ["--" ++ opt, show val]) opts
