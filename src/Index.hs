{-# LANGUAGE DataKinds, TypeOperators, DeriveGeneric #-}
module Index where
import GHC.Generics
import Options.Generic  -- ((<?>), ParseRecord)
import qualified Data.Yaml
import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util
import Data.List (intercalate)
import qualified Data.ByteString.Char8 as B

import Config 

buildIndex' args = shake shakeOptions $ do
--  let gsdb = gsnapDB $ gsnap cfg 
--  let starDir = starDB $ star cfg
--  let bowtiedb = bowtieDB $ bowtie2 cfg
  let dir  = indexDir args

  let gsdb = dir </> "gsdb"
  let starDir = dir </> "stardb"
  let bowtiedb = dir </> "bowitiedb"
  let rsdb = dir </> "rsdb"
  let gscfg = dir </> "gscfg" --gsnapCfg $ gsnap cfg

  let panTro4' = panTro4 args
  let hg38' = hg38 args

  bowtiedb <.> "1.bt21" %> \out -> do
    cmd "bowtie2-build" (intercalate "," [panTro4', hg38']) bowtiedb "--threads" (show $ cpu args) "--large-index"
  
  (gsdb </> (dropDirectory1 gsdb) <.> "version") %> \out -> do
    () <- cmd "mkdir -p" gscfg  -- could be its own rule 
    cmd "gmap_build" "-D" gscfg "-d" gsdb "-k" "16" $ ntFasta args
  
  --(rapsearchDB $ rapsearch cfg) %> \out -> do
  rsdb %> \out -> do
    cmd "prerapsearch" "-d" (nrFasta args) "-n" out
  
  starDir </> "Genome" %> \_ -> do
    () <- cmd "mkdir" starDir
    cmd "STAR --runMode genomeGenerate"
        "--sjdbOverhang 249"
        "--genomeFastaFiles" hg38' panTro4'
        "--runThreadN" (show $ cpu args)
        "--genomeDir" starDir

  
  
  -- hgmap-db.version
  
  --bowtie-index.1.bt2l      bowtie-index.2.bt2l      bowtie-index.3.bt2l      bowtie-index.4.bt2l      bowtie-index.rev.1.bt2l  bowtie-index.rev.2.bt2l
  
  
config' p = Config { threads = 4
                   , star = STAROpts { starDB = p </> "stardb" }
                   , pricefilter = PriceSeqFilterOpts { calledPercent = 95 , highQualPercent = 80 , highQualMin = 0.01 }
                   , lzwfilter = LZWOpts { maxCompressionScore = 55 }
                   , cdhitdup  = CDHitOpts { minDifference = 15 }
                   , bowtie2   = BowtieOpts { bowtieDB = p </> "bowtiedb" }
                   , gsnap   = GSNAPOpts { gsnapDB = p </> "gsdb"} -- , gsnapCfg = p </> "gscfg" }
                   , rapsearch   = RapSearchOpts { rapsearchDB = p </> "rsdb" }
                   , ncbi = Acc2TaxOpts { ntDB = nt, nrDB = nr}
                   , assembler = AssemblerOpts { foo = "bar" }
                   } where 
                      (nr, nt) = ("/", "/")

buildIndex opts = do
  --opts <- getRecord $ T.pack "Running Pathogen.hs" :: IO CommandArgs
  buildIndex' opts
  let cfg = config' $ indexDir opts
  B.putStrLn $ Data.Yaml.encode cfg
