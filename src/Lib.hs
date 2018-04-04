module Lib
  where
import           Data.Array.Unboxed (UArray, amap, array, bounds, (!))
-- import           Data.Array.Unboxed (UArray, amap, array, bounds, indices, ixmap, range, (!))
-- import qualified Data.Array.Unboxed as A
import           Data.List          (transpose)
import           Data.Tuple.Extra   (fst3, snd3, swap, thd3)
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

edges_p1rep_1 :: [Int] -> [Int] -> ([Int],[Int])
edges_p1rep_1 cases p1 =
  (concat edges, concatMap (uncurry replicate) (zip counts p1))
  where
  edges = map head [edgesTable!!(i-1) | i <- cases]
  counts = map length edges
{-
#       count <- sapply(edges, function(x) length(x)) # probablement qu'ici on n'a que des edges vecteurs
#       edges <- cbind(unlist(edges), rep(p1, count))
#                 # ce cbind est inutile, on resépare après
-}

getPoints :: [Int] -> [Int] -> [[Double]] -> [[Double]]
getPoints edges p1 info = out -- correspond à matrix(info, ncol = 8) : CalPoint appliqué là-dessus
  where
  x1 = [edgePoints!!(i-1)!!1 | i <- edges]
  x2 = [edgePoints!!(i-1)!!2 | i <- edges]
  lambda = map (realToFrac . floor . (/9) . fromIntegral) x1
  mu = map (\x -> 1-x) lambda
  average w w' = zipWith (+) (zipWith (*) mu w) (zipWith (*) lambda w')
  v1357 = [info!!(i-2) | i <- zipWith (+) p1 x1]
  v2468 = [info!!(i-2) | i <- zipWith (+) p1 x2]
  v35' = [info!!i | i <- p1]
  v1  = map (!!0) v1357
  v1' = [info!!(i-1)!!0 | i <- p1]
  v2  = map (!!0) v2468
  v2' = [info!!i!!0 | i <- p1]
  v3  = map (!!1) v1357
  v3' = map (!!1) v35'
  v4  = map (!!1) v2468
  v4' = [info!!(i+1)!!1 | i <- p1]
  v5  = map (!!2) v1357
  v5' = map (!!2) v35'
  v6  = map (!!2) v2468
  v6' = [info!!(i+4)!!2 | i <- p1]
  v7  = map (!!3) v1357
  v7' = lambda
  v8  = map (!!3) v2468
  v8' = map negate lambda
  out = [ average v1 v1'
        , average v2 v2'
        , average v3 v3'
        , average v4 v4'
        , average v5 v5'
        , average v6 v6'
        , average v7 v7'
        , average v8 v8'
        ]
{-
  c((1 - floor(x1 / 9)) * info[p1 + x1 - 1, 1] + # v1
      floor(x1 / 9) * info[p1, 1], # v1'
    (1 - floor(x1 / 9)) * info[p1 + x2 - 1, 1] + # v2
      floor(x1 / 9) * info[p1 + 1, 1], # v2'
    (1 - floor(x1 / 9)) * info[p1 + x1 - 1, 2] + # v3
      floor(x1 / 9) * info[p1 + 1, 2], # v3'
    (1 - floor(x1 / 9)) * info[p1 + x2 - 1, 2] + # v4
      floor(x1 / 9) * info[p1 + 2, 2], # v4'
    (1 - floor(x1 / 9)) * info[p1 + x1 - 1, 3] + # v5
      floor(x1 / 9) * info[p1 + 1, 3], # v5'
    (1 - floor(x1 / 9)) * info[p1 + x2 - 1, 3] + # v6
      floor(x1 / 9) * info[p1 + 5, 3], # v6'
    (1 - floor(x1 / 9)) * info[p1 + x1 - 1, 4] + # v7
      floor(x1 / 9) * (0 * info[p1 + 1, 3] + 1),
    (1 - floor(x1/9)) * info[p1 + x2 - 1, 4] + # v8
      floor(x1 / 9) * (0 * info[p1 + 1, 3] - 1))
-}

calPoint :: [[Double]] -> [[Double]]
calPoint info = [scale (info!!0) (info!!1), scale (info!!2) (info!!3), scale (info!!4) (info!!5)]
  where
  s = zipWith (/) (info!!6) (zipWith (-) (info!!6) (info!!7))
  scale u v = zipWith (+) u (zipWith (*) s (zipWith (-) v u))

-- preRender1 :: [[Int]] -> [Int] -> [[Double]] -> [[Double]]

test_levCells = levCells v' 22 48
test_getBasic = getBasic [1,2,3,4,5,6,7] v' 22 test_levCells
test_edges_p1rep = edges_p1rep_1 (thd3 test_getBasic) (snd3 test_getBasic)
test_getPoints = getPoints (fst test_edges_p1rep) (snd test_edges_p1rep) (fst3 test_getBasic)
