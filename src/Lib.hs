module Lib
  where
import           Data.Array.Unboxed (UArray, amap, array, bounds, (!))
-- import           Data.Array.Unboxed (UArray, amap, array, bounds, indices, ixmap, range, (!))
-- import qualified Data.Array.Unboxed as A
import           Data.List          (transpose)
import           Data.Tuple.Extra   (swap)
import           Matrices
import           Tables
import           Utils              (toTriplet)

v' :: UArray (Int,Int,Int) Double
v' = array ((1,1,1),(3,3,3))
            [((i,j,k), x*x+y*y+z*z) | i <- [1..3], j <- [1..3], k <- [1..3],
                                      let x = 1 + fromIntegral i,
                                      let y = 1 + fromIntegral j,
                                      let z = 1 + fromIntegral k]


faceType :: UArray (Int,Int) Double -> Int -> Int -> Double -> Double -> UArray (Int,Int) Int
faceType v nx ny level maxvol = foldr matricialSum v1 [v2,v3,v4]
  where
  comparison = if level == maxvol then fromEnum . (>= level) else fromEnum . (> level)
  v0 = amap comparison v
  v1 = minorMatrix v0 nx ny
  v2 = scaledMatrix 2 (minorMatrix v0 1 ny)
  v3 = scaledMatrix 4 (minorMatrix v0 1 1)
  v4 = scaledMatrix 8 (minorMatrix v0 nx 1)

levCells :: UArray (Int,Int,Int) Double -> Double -> Double
         -> (([Int],[Int],[Int]),[Int])
levCells a level maxvol = -- (concatMap (fst.f) [1 .. (nz-1)], concatMap (snd.f) [1 .. (nz-1)])
  ((v_i, v_j, v_k), concatMap snd cellsAndTypes)
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
  v_i = map ((+1) . (`mod` (nx-1))) cells
  v_j = map ((+1) . (`mod` (ny-1)) . (`div` (nx-1))) cells
  v_k = map ((+1) . (`div` ((nx-1)*(ny-1)))) cells

getBasic :: [Int] -> UArray (Int,Int,Int) Double -> Double -> (([Int],[Int],[Int]),[Int])
         -> ([[Double]], [Int], [Int])
getBasic r vol level ((v_i,v_j,v_k),v_t) = (information, p1, cases)
  where
  cube_1 = transpose [v_i,v_j,v_k]
  index :: [[Int]]
  index = [ [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1] ]
  k1 = concat $ replicate (length v_i) index
  k2 = concatMap (replicate 8) cube_1
  cube_co = zipWith (zipWith (+)) k1 k2
  value = map (subtract level) [vol ! toTriplet ijk | ijk <- cube_co] ++ [0]
  information = transpose $ transpose (map (map fromIntegral) (cube_co ++ [[0,0,0]])) ++ [value]
  -- on verra si c'est bien de concaténer cube_co et value
  p1 = map ((+1) . (*8)) [0 .. length r -1]
  cases = [v_t !! (i-1) | i <- r]

test_levCells = levCells v' 22 48
test_getBasic = getBasic [1,2,3,4,5,6,7] v' 22 test_levCells
