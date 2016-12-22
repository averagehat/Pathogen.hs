module Make where

import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util

-- Makefile for installing pipeline and setting it up 
--
--
--bin = "/media/People/michael.panciera/myconda/bin"
block bin = shakeArgs shakeOptions{shakeFiles="_build"} $ do
  
  let condaExecs = map (bin </>) ["gsnapl", "STAR", "RepeatMasker", "bowtie2"]
  let psf = bin </> "PriceSeqFilter"
  let rs = [bin </> "prerapsearch", bin </> "rapsearch"]
  let cdd = bin </> "cd-hit-dup"

  want $ condaExecs ++ rs ++ [psf, cdd]

  condaExecs &%> \execs -> do
    cmd "conda install" execs
  
  psf %> \exec -> do
    let source = "PriceSource140408"
    () <- wgetAndUnTar $ Url $ "http://derisilab.ucsf.edu/software/price/" ++ source ++ ".tar.gz" 
    cmd (Cwd source) "ln -s PriceSeqFilter" exec
  
  
  rs &%> \[exec1, exec2] -> do
    () <- cmd "git clone https://github.com/zhaoyanswill/RAPSearch2"
    () <- cmd "sh RAPSearch2/install"
    () <- cmd "ln -s RAPSearch2/bin/rapsearch" exec1
    cmd "ln -s RAPSearch2/bin/prerapsearch" exec2
  
  cdd %> \exec -> do
    () <- cmd "git clone https://github.com/weizhongli/cdhit"
    () <- cmd "make -C cdhit/cd-hit-auxtools"
    cmd "ln -s cdhit/cd-hit-auxtools/cd-hit-dup"  exec

newtype Url = Url String
wgetAndUnTar :: Url -> Action () 
wgetAndUnTar (Url url) = do 
  () <- cmd "wget" url
  let tar = reverse $ takeWhile (/= '/') $ reverse url
  cmd "tar" "-xf" tar
