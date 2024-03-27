#!/bin/env python3

import sys
import random

"""
A script that takes a score file in the format used by MSA (Carillo & Lipman), and randomly mutates the scores to produce perturbed alignments.
The way the alignments are mutated can be configured.


"""

def mutate_scores(fin, fout, mutator, mutate_gapcost=False, mutate_gapmatchcosts=True):
    """
    Mutate a scores file.
    `fin` and `fout` describe input and output channels, respectively.
    `mutator` is a lambda that performs the mutation operation on a single score.
    A suitable closure can be constructed using some of the other functions in this file.
    `mutate_gapcost` and `mutate_gapmatchcosts` are flags specifying whether to perform the mutation operation on the gapcost specification or matches involving gaps, respectively.
    They default to False and True, respectively.
    """

    if not mutate_gapcost: # the gapcost is always first; skip mutation if instructed
        print(next(fin), end='', file=fout)
    else:
        print(mutator(
                int(next(fin).strip())
                ), file=fout)

    for line in fin:
        if line.strip() == '' or line[0] == '#': # ignore and re-emit empty and comment lines
            print(line, end='', file=fout)
            continue

        fields = line.strip().split(' ')
        assert(len(fields) == 3) # check that the format is correct

        # check if this field is a match to a gap and skip if specified
        if (not mutate_gapmatchcosts) and any(f == '-' for f in fields):
            print(line, end='', file=fout)
            continue

        print(fields[0], fields[1],
              mutator(int(fields[2])),
              sep=' ', file=fout)
    #end


def add_mutator(vals, chance):
    """
    Returns a mutator function that adds a randomly chosen value of vals to the input `chance` % of the time.
    """
    return lambda x: x + random.choice(vals) if random.randrange(100) < chance else x

def mult_mutator(vals, chance):
    """
    Returns a mutator function that multiplies a randomly chosen value of vals with the input `chance` % of the time.
    """
    return lambda x: round(x * random.choice(vals)) if random.randrange(100) < chance else x


SHUFFLE_STORE = list()
def shuffle_mutator(chance):
    """
    Stores all values passed to it in a global set.
    `chance` % of the time, returns a random value from the set instead of x, sampling from the background with replacement.
    Not thread-safe.
    """
    # somewhat hacky way to call multiple statements in a lambda
    # taken from https://stackoverflow.com/a/69790955
    return lambda x: (SHUFFLE_STORE.append(x),
                      random.choice(SHUFFLE_STORE) if random.randrange(100) < chance
                      else x
                      )[-1]

def clear_shufflers():
    """
    Clears the set used by shuffle_mutators.
    """
    SHUFFLE_STORE = list()

def chain_mutators(lst):
    """
    Helper function to chain different mutator lambdas together into a single lambda.
    Precedence is strictly right-to-left.
    Chaining multiple shuffler lambdas may lead to strange behaviour!
    """
    if len(lst) == 1: # base case for the recursion
        return lst[0]
    return lambda x: lst[0](chain_mutators(lst[1:])(x))


def make_blosum_scores(fout, blosum=82):
    #TODO implement, use python blosum package
    pass


if __name__ == '__main__':
    #import argparse #TODO do some argparsing to configure the perturbation by a script

    #mutate_scores(sys.stdin, sys.stdout, add_mutator([-1, 1], 40))
    with open(sys.argv[1], 'r') as f:
        mutate_scores(f, sys.stdout,
                      chain_mutators([
                          add_mutator([-1, 1], 50),
                          mult_mutator([0.8, 1.2], 30),
                          shuffle_mutator(10)
                          ]))
        #mutate_scores(f, sys.stdout, add_mutator([-1, 1], 40))

