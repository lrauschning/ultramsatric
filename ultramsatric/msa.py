import copy
import numpy as np
from typing import Dict, List
import os

class MSA:
    alns: Dict[str, List[chr]] # store fasta as mapping of ID to sequence

    def subset(self, subs) -> MSA:
        if not subs.issubset(self.alns.keys()):
            raise ValueError(f"{subs} is not a subset of {self.alns.keys}!")
        return MSA(alns)

    def from_file(path: os.PathLike) -> MSA:
        """Parses a FASTA file into a MSA object.
        """
        alns = dict()
        with open(path, 'rt') as f:
            curid = ''
            seq = list()
            for l in f:
                if l[0] == '>':
                    alns[curid] = seq # store the sequence we have so far
                    curid = l.split(' ')[0][1:] # extract new ID
                    seq = list() # reset sequence buffer
                    continue
                seq += list(l.strip())
            alns[curid] = seq

        return MSA(alns)



   # higher-order function that takes a function from the distance module to turn a MSA into an induced distance matrix
   def to_distance(distfun) -> np.ndarray:
