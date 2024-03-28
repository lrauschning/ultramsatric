#!/bin/env python3

import random

"""
A script that takes a score file in the format used by MSA (Carillo & Lipman), and randomly mutates the scores to produce perturbed alignments.
The way the alignments are mutated can be configured; lines starting with '#' in the input scores are considered comments and re-emitted.


"""

def mutate_scores(fin, fout, mutator, mutate_gapcost=False, mutate_gapmatchcosts=False):
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
    import argparse as ap

    parser = ap.ArgumentParser(description="""
    MutateScores – mutates your favourite alignment scores beyond recognition
    """)
    parser.add_argument("-i", dest='infile', default='-', type=ap.FileType('r'),
                        help="Input Scores file. Uses the format used by MSA (Carillo & Lipmann) with each line corresponding to one pairing and its score, with the first line containing the gapcost. Default stdin.")
    parser.add_argument("-o", dest='outfile', default='-', type=ap.FileType('wt'),
                        help="File to write output scores to. Default stdout.")
    parser.add_argument("-m", "--mutators", dest='mutators', default="add_mutator([-1, 1], 50);mult_mutator([0.8, 1.2], 30);shuffle_mutator(10)",
                        help="Mutators to apply to the scores. Takes python expressions separated by ';'. The mutators are applied right-to-left. The constructors defined in this file may be useful to look at. By default, the first mutator shuffles the values with a 10% chance, the second multiplies by 1.2 or 0.8 with a chance of 30%, and the third adds 1 or -1 with a chance of 50%.")
    parser.add_argument("--mutate-gapcost", dest='mutate_gapcost', action='store_true', default=False,
                        help="Whether to apply the mutation to the gapcost specified in the first line of the file. Default False.")
    parser.add_argument("--mutate-gapmatches", dest='mutate_gapmatches', action='store_true', default=False,
                        help="Whether to apply the mutation to matches between a residue and a gap. Default False.")

    
    args = parser.parse_args()

    mutator = chain_mutators([eval(mut) for mut in args.mutators.split(';')])

    mutate_scores(args.infile, args.outfile, mutator, mutate_gapcost=args.mutate_gapcost, mutate_gapmatchcosts=args.mutate_gapmatches)



    #import sys
    #mutate_scores(sys.stdin, sys.stdout, add_mutator([-1, 1], 40))
    #with open(sys.argv[1], 'r') as f:
    #    mutate_scores(f, sys.stdout,
    #                  chain_mutators([
    #                      add_mutator([-1, 1], 50),
    #                      mult_mutator([0.8, 1.2], 30),
    #                      shuffle_mutator(10)
    #                      ]))
    #    #mutate_scores(f, sys.stdout, add_mutator([-1, 1], 40))

