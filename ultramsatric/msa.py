import copy
import numpy as np

class MSA:
    alns: Dict[str, str] # store fasta as mapping of ID to sequence

    def subset(self, subs) -> MSA:
        if not subs.issubset(self.alns.keys()):
            raise ValueError(f"{subs} is not a subset of {self.alns.keys}!")
        return MSA(alns)

   # higher-order function that takes a function from the distance module to turn a MSA into an induced distance matrix
   def to_distance(distfun) -> np.ndarray:
