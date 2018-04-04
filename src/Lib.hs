module Lib
  where
import           Data.Array.Unboxed (UArray, amap, array, bounds, (!))
-- import           Data.Array.Unboxed (UArray, amap, array, bounds, indices, ixmap, range, (!))
-- import qualified Data.Array.Unboxed as A
import           Data.List          (transpose, zipWith4)
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
         -> (([Double], [[Double]]), [Int], [Int])
getBasic r vol level ((v_i,v_j,v_k),v_t) = ((values, information), p1, cases)
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
  values = map (subtract level) [vol ! toTriplet ijk | ijk <- cube_co] ++ [0]
  information = transpose $ transpose (map (map fromIntegral) (cube_co ++ [[0,0,0]])) ++ [values]
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

calPoint :: [[Double]] -> [[Double]]
calPoint info = [scale (info!!0) (info!!1), scale (info!!2) (info!!3), scale (info!!4) (info!!5)]
  where
  s = zipWith (/) (info!!6) (zipWith (-) (info!!6) (info!!7))
  scale u v = zipWith (+) u (zipWith (*) s (zipWith (-) v u))

preRender1 :: [Int] -> [Int] -> [[Double]] -> [[Double]]
preRender1 cases p1 information = transpose $ calPoint info
  where
  (edges, p1rep) = edges_p1rep_1 cases p1
  info = getPoints edges p1rep information

faceNo7 :: [Int] -> [Int] -> [Double] -> [Int]
faceNo7 faces p1 values = map (\i -> if i==1 then 1 else 0) index
  -- info n'intervient que par sa colonne 3, i.e. values
  where
  sign :: (Num a, Ord a) => a -> Int
  sign x = if x>0 then 1 else -1
  faces' = map abs faces
  e  = [facePoints!!(i-1) | i <- faces']
  e1 = map (!!1) e
  e2 = map (!!2) e
  e3 = map (!!3) e
  e4 = map (!!4) e
  a = [values!!(i-2) | i <- zipWith (+) p1 e1]
  b = [values!!(i-2) | i <- zipWith (+) p1 e2]
  c = [values!!(i-2) | i <- zipWith (+) p1 e3]
  d = [values!!(i-2) | i <- zipWith (+) p1 e4]
  abMINUScd = zipWith (-) (zipWith (*) a b) (zipWith (*) c d)
  index = zipWith (*) (map sign faces) (map sign abMINUScd)

face7 :: [Int] -> [Int] -> [Double] -> [Int]
face7 faces p1 values = map (\i -> if i==1 then 1 else 0) index
  where
  a0 = [values!!(i-1) | i <- p1]
  b0 = [values!!(i+2) | i <- p1]
  c0 = [values!!(i+1) | i <- p1]
  d0 = [values!!i | i <- p1]
  a1 = [values!!(i+3) | i <- p1]
  b1 = [values!!(i+6) | i <- p1]
  c1 = [values!!(i+5) | i <- p1]
  d1 = [values!!(i+4) | i <- p1]
  a1a0 = zipWith (-) a1 a0
  b1b0 = zipWith (-) b1 b0
  c1c0 = zipWith (-) c1 c0
  d1d0 = zipWith (-) d1 d0
  a = zipWith (-) (zipWith (*) a1a0 c1c0) (zipWith (*) b1b0 d1d0)
  c = zipWith (-) (zipWith (*) a0 c0) (zipWith (*) b0 d0)
  b' =  zipWith (-) (zipWith (*) a1a0 c0) (zipWith (*) c1c0 a0)
  b'' = zipWith (-) (zipWith (*) b1b0 d0) (zipWith (*) d1d0 b0)
  b = zipWith (-) b' b''
  tmax = zipWith (/) b (map (*(-2)) a)
  tmax2 = zipWith (*) tmax tmax
  maxi = zipWith (+) (zipWith (+) (zipWith (*) a tmax2) (zipWith (*) b tmax)) c
  maxim = map (\x -> if isNaN x then -1 else x) maxi
{-  cond1 = map (\x -> if x<0 then 1 else 0) a
  cond2 = map (\x -> if x>0 then 1 else 0) tmax
  cond3 = map (\x -> if x<1 then 1 else 0) tmax
  cond4 = map (\x -> if x>0 then 1 else 0) maxim -}
  cond1 = map (<0) a
  cond2 = map (>0) tmax
  cond3 = map (<1) tmax
  cond4 = map (>0) maxim
  totalcond = zipWith4 (\s t u v -> and [s,t,u,v]) cond1 cond2 cond3 cond4
  index = zipWith (*) (map (\x -> if x>0 then 1 else -1) faces)
                      (map (\x -> if x then 1 else -1) totalcond)

getInfo = snd . fst3
test_levCells = levCells v' 22 48
test_getBasic = getBasic [1,2,3,4,5,6,7] v' 22 test_levCells
test_edges_p1rep = edges_p1rep_1 (thd3 test_getBasic) (snd3 test_getBasic)
test_getPoints = getPoints (fst test_edges_p1rep) (snd test_edges_p1rep) (getInfo test_getBasic)
test_preRender1 = preRender1 (thd3 test_getBasic) (snd3 test_getBasic) (getInfo test_getBasic)
