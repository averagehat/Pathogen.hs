module Main where

--import System.FilePath.Glob (glob)
import Test.DocTest (doctest)

main = do
  doctest ["src/Config.hs", "src/Pipeline.hs"]
--main = glob "src/**/*.hs" >>= doctest
