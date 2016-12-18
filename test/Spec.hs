{-# LANGUAGE OverloadedStrings #-}
module Main where
import Test.Framework (defaultMain, testGroup)
import Config
import Test.Framework.Providers.HUnit
import Test.HUnit

--main :: IO ()
--main = putStrLn "No Tests"

main = defaultMain [tests]

tests = testGroup "HUnit tests" [ testCase "a config test!"  testTemplateParse
                                , testCase "another passing test!" $ 6 @?= 6 ]

testTemplateParse = do
  cfg <- readYaml "template.yaml"
  let fakepath = bowtieDB <$> bowtie2 <$> cfg
  fakepath @?= Right "/fakepath/"


--tests = [testGroup "Specifics" [ testProperty "config examples" prop_example_expand ] ]
