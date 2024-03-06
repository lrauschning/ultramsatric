from typing import Dict, List
import os

import numpy as np

class MSA:
    def __init__(self, alns:Dict[str, List[chr]]):
        self.alns = alns # store fasta as mapping of ID to sequence

    def subset(self, subs):
        if not subs.issubset(self.alns.keys()):
            raise ValueError(f"{subs} is not a subset of {self.alns.keys}!")
        return MSA({x:self.alns[x] for x in subs})

    @classmethod
    def from_file(cls, path: os.PathLike):
        """Parses a FASTA file into a MSA object.
        """
        with open(path, 'rt') as f:
            return cls.from_inputstream(f)

    @classmethod
    def from_inputstream(cls, stream):
        """Parses a FASTA stream into a MSA object.
        """
        alns = dict()
        curid = ''
        seq = list()
        for l in stream:
            l = l.strip()
            if len(l) == 0: # skip empty lines
                continue
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

    def sumofpairs(self, distfun) -> float:
        return sum([sum([
                distfun(self.alns[x], self.alns[y])
                # avoid double-counting by requiring this bound
                for y in self.alns if y > x ])
            for x in self.alns])

    def sumofpairs_avg(self, distfun) -> float:
        return self.sumofpairs(distfun)*2/(len(self.alns)**2 - len(self.alns))

    def totalcol(self, subs, gapcost, use_gaplen=True) -> float:
        """
        Calculates the total column score of a MSA.
        This implementation is able to use non-linear gapcosts by looking at the neighbouring columns.
        If disabled by passing `use_gaplen=False`, all `gapcost` calls are made as `gapcost(1)`.
        Gap lengths do not include lengths of gaps shared between the two columns that are being compared.
        This score is symmetric, i.e. for any comparison `subs(a, b)`, there will also be another comparison `subs(b, a)`. This is not the case for gaps -- `gapcost` is called only once per sequence that gap is in.

        :args: substitution model, gapcost function
        :returns: total column score
        :rtype: float
        """
        tc = 0.0
        alns = list(self.alns.values())
        for i in range(len(alns[0])): # iterate through columns
            for a in alns:
                for b in alns:
                    if a == b or (a[i] == '-' and b[i] == '-'): # do not compare a to a, or gaps to each other
                        continue

                    if a[i] == '-':
                        # add gapcost if there is a gap in the a strand;
                        # the b strand will be caught when a is b and b is a
                        gaplen = 1

                        if use_gaplen:
                            if i > 0:
                                if (a[i] == '-' and a[i-1] == '-' ):
                                    continue # this gap is not new, and has been added before

                            # find the gaplen by looking for an uninterrupted gap in a; ignore cases that are also gaps in b
                            for j in range(i, len(alns[0])):
                                if a[j] != '-':
                                    break
                                if b[j] == '-':
                                   continue
                                # the ordering is important here,
                                # in the case below, the first gap would be counted with length 4 otherwise
                                # a AC--TG--
                                # b ACAA--TC
                                gaplen += 1

                        # add gapcost, indepent of if using gaplen or not
                        tc += gapcost(gaplen)

                    else:
                        tc += subs(a[i], b[i])
        return tc
