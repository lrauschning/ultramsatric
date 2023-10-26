import copy
import numpy as np
from typing import Dict, List
import os

class MSA:
    def __init__(self, alns:Dict[str, List[chr]]):
        self.alns = alns # store fasta as mapping of ID to sequence

    @classmethod
    def subset(cls, subs):
        if not subs.issubset(self.alns.keys()):
            raise ValueError(f"{subs} is not a subset of {self.alns.keys}!")
        return cls(alns[subs])

    @classmethod
    def from_file(cls, path: os.PathLike):
        """Parses a FASTA file into a MSA object.
        """
        alns = dict()
        with open(path, 'rt') as f:
            curid = ''
            seq = list()
            for l in f:
                l = l.strip()
                if l[0] == '>':
                    if len(seq) > 0: # avoid adding empty line
                        alns[curid] = seq # store the sequence we have so far
                    curid = l.split(' ')[0][1:].strip() # extract new ID
                    seq = list() # reset sequence buffer
                    continue
                seq += list(l.strip())

            if len(seq) > 0: # avoid adding empty line
                alns[curid] = seq
        return cls(alns)

    def __repr__(self) -> str:
        return '\n'.join([f">{id}\n{seq}" for id, seq in self.alns.items()])
