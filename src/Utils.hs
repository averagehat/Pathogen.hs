{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DeriveGeneric #-} 
{-# LANGUAGE RecordWildCards #-}
module Utils where
import GHC.Generics (Generic())
import Data.Csv
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Data.Text (Text())
import qualified Data.ByteString.Lazy.Char8 as C8
import Data.Vector (Vector())
import qualified Data.Vector as V
import Data.Char (ord)
import Control.Monad (mzero)
import Data.List (nub,groupBy)
data M8 = M8 {
   query       :: !Text
 , subject     :: !Text
 , identity    :: !Double
 , alnLen      :: !Double
 , mismatch    :: !Int
 , gapOpenings :: !Int
 , qstart      :: !Int -- 0-based
 , qend        :: !Int
 , sstart      :: !Int
 , send        :: !Int
 , eValue      :: !Double
 , bitScore    :: !Double
 } deriving (Show, Eq)
-- Fields: Query	Subject	identity	aln-len	mismatch	gap-openings	q.start	q.end	s.start	s.end	log(e-value)	bit-score
instance FromNamedRecord M8 where
    parseNamedRecord m = M8 <$>
                         m .: "Query" <*>
                         m .: "Subject" <*>
                         m .: "identity" <*>
                         m .: "aln-len" <*>
                         m .: "mismatch" <*>
                         m .: "gap-openings" <*>
                         m .: "q.start" <*>
                         m .: "q.end" <*>
                         m .: "s.start" <*>
                         m .: "s.end" <*>
                         m .: "log(e-value)" <*>
                         m .: "bit-score"
type Error = String
tabDec = defaultDecodeOptions {      decDelimiter = fromIntegral (ord '\t') }
tabEnc = defaultEncodeOptions { encDelimiter = fromIntegral $ ord '\t' }
showM8 :: M8 -> String
showM8 m = show m

readM8 :: FilePath -> IO (Either Error (Vector M8))
readM8 p = decodeWith tabDec HasHeader . C8.pack . T.unpack . last . T.splitOn "# Fields: " <$> T.readFile p

readTax :: FilePath -> IO (Either Error (Vector Taxonomy))
readTax p = decodeWith tabDec NoHeader <$> C8.readFile p

combine' :: [M8] -> [Taxonomy] -> [WithTax]
combine' m8 tax = zipWith (combine'' hasOdd) m8 tax where
  hasOdd = (length $ nub $ fmap genus tax) > 1

combine :: Vector M8 -> Vector Taxonomy -> [WithTax]
combine m8 tax = let m8s = groupBy (\x y -> query x == query y) $ V.toList m8
                     taxs = takeLengths (V.toList tax) m8s in
                 concatMap (uncurry combine') $ zip m8s taxs
                 where takeLengths xs ys = fst $ foldl f ([], xs) ys
                       f (x, y) a = let (x', y') = splitAt (length a) y in (x ++ [x'], y')

  
combine'' :: Bool -> M8 -> Taxonomy -> WithTax
combine'' odd M8{..} Taxonomy{..} = WithTax {
    _taxid=taxid
    , _phyla=phyla
    , _kingdom=kingdom
    , _order=order
    , _genus=genus
    , _species=species
    , _query=query
    , _subject=subject
    , _identity=identity
    , _alnLen=alnLen
    , _mismatch=mismatch
    , _gapOpenings=gapOpenings
    , _qstart=qstart
    , _qend=qend
    , _sstart=sstart
    , _send=send
    , _eValue=eValue
    , _bitScore=bitScore
    , oddGenus=odd}

readWithTax :: FilePath -> IO (Either Error (Vector WithTax))
readWithTax p =   do
    s <- C8.readFile p
    return $ fmap snd $ decodeByNameWith tabDec s 

-- e.g.
--  "rapsearch.tsv.top" %> \out -> do
--    src <- needExt "" out
--    topBlast src out

topBlast :: FilePath -> FilePath -> IO ()
topBlast fp out = do
  xs <- readWithTax fp
  let grps = groupBy (\x y -> _query x == _query y) . V.toList <$> xs
  let tops = fmap head <$> grps
  either error (C8.writeFile out . encodeDefaultOrderedByNameWith tabEnc) tops 

addTax :: FilePath -> FilePath -> FilePath -> IO ()                              
addTax m8 tax out = do
  m8' <- readM8 m8
  tax' <- readTax tax 
  let xs = combine <$> m8' <*> tax'
  case xs of
    (Right xs') -> C8.writeFile out $ encodeDefaultOrderedByNameWith tabEnc xs'
    (Left err) -> error err

instance FromRecord M8 where
    parseRecord m = M8 <$>
                         m .! 0 <*>
                         m .! 1 <*>
                         m .! 2 <*>
                         m .! 3 <*>
                         m .! 4 <*>
                         m .! 5 <*>
                         m .! 6 <*>
                         m .! 7 <*>
                         m .! 8 <*>
                         m .! 9 <*>
                         m .! 10 <*>
                         m .! 11
                         
data Taxonomy = Taxonomy {
    taxid   :: !Int
  , phyla   :: !Text
  , kingdom :: !Text
  , order   :: !Text
  , genus   :: !Text
  , species :: !Text
  } deriving (Show, Eq)
instance FromRecord Taxonomy where
  parseRecord m = do
    m' <- fmap toField . V.fromList . T.splitOn "," <$> (m .! 1)
    Taxonomy <$>
     m  .! 0 <*>
     m' .! 0 <*>
     m' .! 1 <*>
     m' .! 2 <*>
     m' .! 3 <*>
     m' .! 4 
     
data WithTax = WithTax {
   _query      :: !Text
 , _subject     :: !Text
 ,  _taxid      :: !Int
 , _phyla      :: !Text
 , _kingdom    :: !Text
 , _order      :: !Text
 , _genus      :: !Text
 , _species    :: !Text
 , _identity    :: !Double
 , _alnLen      :: !Double
 , _mismatch    :: !Int
 , _gapOpenings :: !Int
 , _qstart      :: !Int -- 0-based
 , _qend        :: !Int
 , _sstart      :: !Int
 , _send        :: !Int
 , _eValue      :: !Double
 , _bitScore    :: !Double
 , oddGenus     :: Bool
 } deriving (Show, Eq, Generic)

instance ToField Bool where
  toField = toField . show
instance FromField Bool where
  parseField "True" = pure True
  parseField "False" = pure False
  parseField _       = mzero
  
  
instance ToNamedRecord WithTax
instance DefaultOrdered WithTax
instance FromNamedRecord WithTax
