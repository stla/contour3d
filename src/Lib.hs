{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Lib
  where
import           Data.Array.IArray  (IArray)
import           Data.Array.Unboxed (UArray, amap, array, bounds, indices,
                                     ixmap, range, (!))
import qualified Data.Array.Unboxed as A
import           Data.List
import           Data.Tuple.Extra   (swap)

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

faceType :: UArray (Int,Int) Double -> Int -> Int -> Double -> Double -> UArray (Int,Int) Int
faceType v nx ny level maxvol = foldr matricialSum v1 [v2,v3,v4]
  where
  compare = if level == maxvol then fromEnum . (>= level) else fromEnum . (> level)
  v0 = amap compare v
  v1 = minorMatrix v0 nx ny
  v2 = scaledMatrix 2 (minorMatrix v0 1 ny)
  v3 = scaledMatrix 4 (minorMatrix v0 1 1)
  v4 = scaledMatrix 8 (minorMatrix v0 nx 1)

findIndicesAndItems :: (a -> Bool) -> [a] -> [(Int, a)]
findIndicesAndItems predicate = filter (predicate . snd) . zip [0..]

findIndicesAsIntegers :: IArray UArray a => (a -> Bool) -> UArray (Int,Int) a
                      -> ([Int], [(Int,Int)])
findIndicesAsIntegers predicate m =
    unzip $ findIndicesAndItems (\ij -> predicate (m ! swap ij)) (indices m)
{-    (integers, indices')
  where
  allIndices = indices m
  integers = findIndices (\ij -> predicate (m ! swap ij)) allIndices
  indices' = [allIndices!!i | i <- [0 .. length allIndices - 1], i `elem` integers] -}

levCells :: UArray (Int,Int,Int) Double -> Double -> Double
         -> ([Int],[Int],[Int],[Int])
levCells a level maxvol = -- (concatMap (fst.f) [1 .. (nz-1)], concatMap (snd.f) [1 .. (nz-1)])
  (i, j, k, concatMap snd cellsAndTypes)
  where
  ((_,_,_),(nx,ny,nz)) = bounds a
  types = map (\k -> faceType (toMatrix a k) nx ny level maxvol) [1 .. nz]
  f k = (map (+ (nx-1)*(ny-1)*(k-1)) contourCells, -- inutile de mettre +1
         [cellTypes ! swap ij | ij <- snd intind])
    where
    cellTypes = matricialSum (types!!(k-1)) (scaledMatrix 16 (types!!k))
    intind = findIndicesAsIntegers (\x -> x>0 && x<255) cellTypes
    contourCells = fst intind
  cellsAndTypes = map f [1 .. (nz-1)]
  cells = concatMap fst cellsAndTypes
  i = map ((+1) . (`mod` (nx-1))) cells
  j = map ((+1) . (`mod` (ny-1)) . (`div` (nx-1))) cells
  k = map ((+1) . (`div` ((nx-1)*(ny-1)))) cells
