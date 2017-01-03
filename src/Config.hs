{-# LANGUAGE OverloadedStrings, DeriveGeneric #-}
{-# LANGUAGE DataKinds, TypeOperators #-}  -- these two needed for help output
module Config where
import GHC.Generics
import Options.Generic  -- ((<?>), ParseRecord)
import qualified Data.Yaml
import Data.Aeson hiding (encode)
import qualified Data.ByteString.Char8 as B
import Data.Maybe 
-- TODO: give everything lenses via (Aeson's?) makeLens?
type Fasta = FilePath
data CommandArgs = Run {
    r1 :: FilePath <?> "Forward"
  , r2 :: FilePath <?> "Reverse"
  , i1 :: Maybe FilePath <?> "Forward Index"
  , i2 :: Maybe FilePath <?> "Reverse Index"
  , outdir :: Maybe FilePath <?> "Output Directory"
  , config :: FilePath <?> "Path to config.yaml" }
  | Index {
  cpu :: Int
, nrFasta :: Fasta
, panTro4 :: Fasta
, hg38    :: Fasta
, ntFasta :: Fasta
, indexDir  :: FilePath } 
  deriving (Generic, Show)
instance ParseRecord CommandArgs

type Percent = Int
data Config = Config {threads :: Int
                     , star :: STAROpts
                     , pricefilter :: PriceSeqFilterOpts
                     , cdhitdup :: CDHitOpts
                     , lzwfilter :: LZWOpts
                     , bowtie2    :: BowtieOpts
                     , gsnap     :: GSNAPOpts
                     , rapsearch :: RapSearchOpts
                     , ncbi :: Acc2TaxOpts
                     , assembler :: AssemblerOpts }
  deriving (Generic, Show)
instance FromJSON Config
instance ToJSON Config -- for generateing config.yaml.example. but requires all others to have toJSON as well.

data STAROpts = STAROpts { starDB :: FilePath }
  deriving (Generic, Show)
instance FromJSON STAROpts
instance ToJSON STAROpts

data PriceSeqFilterOpts = PriceSeqFilterOpts { calledPercent :: Int -- <= 100
                                  , highQualPercent :: Int -- <=100
                                  , highQualMin :: Float } -- <= 1
  deriving (Generic, Show)
instance FromJSON PriceSeqFilterOpts
instance ToJSON PriceSeqFilterOpts

data CDHitOpts = CDHitOpts { minDifference :: Percent } -- <= 1
  deriving (Generic, Show)
instance FromJSON CDHitOpts
instance ToJSON CDHitOpts

data LZWOpts = LZWOpts     { maxCompressionScore :: Percent }
  deriving (Generic, Show)
instance FromJSON LZWOpts
instance ToJSON LZWOpts

-- note multiple db accessors not allowed
data BowtieOpts = BowtieOpts { bowtieDB :: FilePath } deriving (Generic, Show)
instance FromJSON BowtieOpts
instance ToJSON BowtieOpts

data GSNAPOpts = GSNAPOpts   { gsnapDB :: FilePath }
  deriving (Generic, Show)
instance FromJSON GSNAPOpts
instance ToJSON GSNAPOpts

data RapSearchOpts = RapSearchOpts { rapsearchDB :: FilePath }
  deriving (Generic, Show)
instance FromJSON RapSearchOpts
instance ToJSON RapSearchOpts

data Acc2TaxOpts = Acc2TaxOpts { nrDB :: FilePath, ntDB :: FilePath }
  deriving (Generic, Show)
instance FromJSON Acc2TaxOpts
instance ToJSON Acc2TaxOpts

data AssemblerOpts = AssemblerOpts { foo :: String }
  deriving (Generic, Show)
instance FromJSON AssemblerOpts
instance ToJSON AssemblerOpts
-- could have the index-building script output a yaml file??
type Error = String
readYaml :: FilePath -> IO (Either Error Config)
readYaml fp = do
  cfgFile <- B.readFile fp
  let cfg = (Data.Yaml.decodeEither cfgFile :: Either Error Config)
  return cfg

--printConfig = do
--  cfg <- readYaml "template.yaml"
--  --cfg <- ((decodeEither <$> readFile "config.yaml") :: IO (Either String Config))
----  bimap  putStrLn putStrLn $ (bowtieDB . bowtie2) <$> cfg
----  either putStrLn putStrLn $ (bowtieDB . bowtie2) <$> cfg
--  -- is equivalent to the below
--  case cfg of
--    Right cfg' -> putStrLn $ show $ bowtieDB $ bowtie2 $ cfg'
--    Left err -> putStrLn err
--  putStrLn $ show $ Data.Yaml.encode config'
--  Data.Yaml.encodeFile "template.yaml" config'

-- have to get input1 and input2 from command args
-- have to get threads and gd from config file
