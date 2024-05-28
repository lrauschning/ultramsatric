from ultramsatric.distance import *

def test_indexing():
    for n in [10, 15, 20]:
        for x in [0, 5, 7, 10, 15, 20, 30]:
            if x < (n*(n-1))//2:
                assert(DistMat.index(*DistMat.revindex(x, n), n) == x)

def test_len():
    pass


