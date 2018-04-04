{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Matrices
  where
import           Data.Array.IArray  (IArray)
import           Data.Array.Unboxed (UArray, amap, array, bounds, indices,
                                     ixmap, range, (!))
import           Data.Tuple.Extra   (swap)
import           Utils              (findIndicesAndItems)

toMatrix :: UArray (Int,Int,Int) Double -> Int  -> UArray (Int,Int) Double
toMatrix a k = ixmap ((1,1), (m,n)) (\(i,j) -> (i,j,k)) a
  where
  (m,n,_) = snd (bounds a)

minorMatrix :: IArray UArray a => UArray (Int,Int) a -> Int -> Int -> UArray (Int,Int) a
minorMatrix mat r c = ixmap ((1,1), (m-1,n-1)) f mat
  where
  (m,n) = snd (bounds mat)
  f (i,j) = (g r i, g c j)
    where
    g h k = if k<h then k else k+1

matricialSum :: (Num a, IArray UArray a) => UArray (Int,Int) a -> UArray (Int,Int) a -> UArray (Int,Int) a
matricialSum m1 m2 = array bds [(ij, m1!ij + m2!ij) | ij <- range bds]
  where
  bds = bounds m1

scaledMatrix :: (Num a, IArray UArray a) => a -> UArray (Int,Int) a -> UArray (Int,Int) a
scaledMatrix k = amap (* k)

findIndicesAsIntegers :: IArray UArray a => (a -> Bool) -> UArray (Int,Int) a
                      -> ([Int], [(Int,Int)])
findIndicesAsIntegers predicate m =
    unzip $ findIndicesAndItems (\ij -> predicate (m ! swap ij)) (indices m)
{-    (integers, indices')
  where
  allIndices = indices m
  integers = findIndices (\ij -> predicate (m ! swap ij)) allIndices
  indices' = [allIndices!!i | i <- [0 .. length allIndices - 1], i `elem` integers] -}
