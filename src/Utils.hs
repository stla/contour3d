module Utils
  where

findIndicesAndItems :: (a -> Bool) -> [a] -> [(Int, a)]
findIndicesAndItems predicate = filter (predicate . snd) . zip [0..]

toTriplet :: [a] -> (a,a,a)
toTriplet [x,y,z] = (x,y,z)
toTriplet _       = undefined
