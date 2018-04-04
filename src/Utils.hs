module Utils
  where

findIndicesAndItems :: (a -> Bool) -> [a] -> [(Int, a)]
findIndicesAndItems predicate = filter (predicate . snd) . zip [0..]
