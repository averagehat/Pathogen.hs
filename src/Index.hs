{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds, TypeOperators, DeriveGeneric #-}
module Index where
import GHC.Generics
import Options.Generic  -- ((<?>), ParseRecord)
import qualified Data.Yaml
import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util

import qualified Data.ByteString.Char8 as B

import Config 

data Inputs = Inputs {
  nrFasta :: Fasta
, panTro4 :: Fasta
, hg38    :: Fasta
, ntFasta :: Fasta
} deriving (Generic, Show)
instance ParseRecord Inputs

block cfg = shakeArgs shakeOptions{shakeFiles="_build"} $ do
  let gsdb = gsnapDB $ gsnap cfg 
  let starDir = starDB $ star cfg
  let bowtiedb = bowtieDB $ bowtie2 cfg
  bowtiedb <.> "1.bt21" %> \out -> do
    cmd "bowtie2-build" (intercalate "," panTro4 hg38) bowtiedb "--threads" (show $ threads cfg) "--large-index"
  
  (gsdb </> (dropDirectory gsdb) <.> "version") %> \out -> do
    let gscfg = gsnapCfg $ gsnap cfg
    () <- cmd "mkdir -p" gscfg  -- could be its own rule 
    cmd "gmap_build" "-D" gscfg "-d" gsdb "-k" "16" ntFasta
  
  (rapsearchDB $ rapsearch cfg) %> \out -> do
    cmd "prerapsearch" "-d" nrFasta "-n" out
  
  starDir </> "Genome" %> \_ -> do
    () <- cmd "mkdir" starDir
    cmd "STAR --runMode genomeGenerate"
        "--sjdbOverhang 249"
        "--genomeFastaFiles" hg38 panTro4
        "--runThreadN" (show $ threads cfg)
        "--genomeDir" starDir
  
  
  -- hgmap-db.version
  
  
  
  
  
  --bowtie-index.1.bt2l      bowtie-index.2.bt2l      bowtie-index.3.bt2l      bowtie-index.4.bt2l      bowtie-index.rev.1.bt2l  bowtie-index.rev.2.bt2l
  
  
config' = Config { threads = 4
                   , star = STAROpts { starDB = p </> "stardb" }
                   , pricefilter = PriceSeqFilterOpts { calledPercent = 95 , highQualPercent = 80 , highQualMin = 0.01 }
                   , lzwfilter = LZWOpts { maxCompressionScore = 55 }
                   , cdhitdup  = CDHitOpts { minDifference = 15 }
                   , bowtie2   = BowtieOpts { bowtieDB = p </> "bowtiedb" }
                   , gsnap   = GSNAPOpts { gsnapDB = p </> "gsnapdb", gsnapCfg = p </> "gsnapcfg" }
                   , rapsearch   = RapSearchOpts { rapsearchDB = p </> "rapsearchdb" }
                   , ncbi = Acc2TaxOpts { ntDB = nt, nrDB = nr}
                   , assembler = AssemblerOpts { foo = "bar" }
                   } where 
                      p = "/fakepath/"
                      (nr, nt) = ("/", "/")
