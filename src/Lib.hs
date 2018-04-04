{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Lib
  where
import           Data.Array.IArray  (IArray)
import           Data.Array.Unboxed (UArray, amap, array, bounds, ixmap, range,
                                     (!))
import qualified Data.Array.Unboxed as A

v' :: UArray (Int,Int,Int) Double
v' = array ((1,1,1),(3,3,3))
            [((i,j,k), x*x+y*y+z*z) | i <- [1..3], j <- [1..3], k <- [1..3],
                                      let x = 1 + fromIntegral i,
                                      let y = 1 + fromIntegral j,
                                      let z = 1 + fromIntegral k]

toMatrix :: UArray (Int,Int,Int) Double -> Int  -> UArray (Int,Int) Double
toMatrix a k = ixmap ((1,1), (m,n)) (\(i,j) -> (i,j,k)) a
  where
  (m,n,_) = snd (bounds a)

minorMatrix :: IArray UArray a => UArray (Int,Int) a -> Int -> Int -> UArray (Int,Int) a
minorMatrix a r c = ixmap ((1,1), (m-1,n-1)) f a
  where
  (m,n) = snd (bounds a)
  f (i,j) = (g r i, g c j)
    where
    g h k = if k<h then k else k+1

matricialSum :: (Num a, IArray UArray a) => UArray (Int,Int) a -> UArray (Int,Int) a -> UArray (Int,Int) a
matricialSum m1 m2 = array bds [(ij, m1!ij + m2!ij) | ij <- range bds]
  where
  bds = bounds m1

scaledMatrix :: (Num a, IArray UArray a) => a -> UArray (Int,Int) a -> UArray (Int,Int) a
scaledMatrix k = amap (* k)

faceType :: UArray (Int,Int) Double -> Int -> Int -> Double -> Double -> UArray (Int,Int) Int
faceType v nx ny level maxvol = foldr matricialSum v1 [v2,v3,v4]
  where
  compare = if level == maxvol then fromEnum . (>= level) else fromEnum . (> level)
  v0 = amap compare v
  v1 = minorMatrix v0 nx ny
  v2 = scaledMatrix 2 (minorMatrix v0 1 ny)
  v3 = scaledMatrix 4 (minorMatrix v0 1 1)
  v4 = scaledMatrix 8 (minorMatrix v0 nx 1)
