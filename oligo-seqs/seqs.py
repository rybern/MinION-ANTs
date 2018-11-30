import itertools
import sys
from functools import reduce
from barcode_filtering import *

NNTs = "ACGT"
ANTs = "X"
NTs = NNTs + ANTs

# from Page 20 of
# http://science.sciencemag.org.ezproxy.cul.columbia.edu/content/sci/suppl/2017/03/01/355.6328.950.DC1/Erlich.SM.pdf
PCR_prefix = "GTTCAGAGTTCTACAGTCCGACGATC"
PCR_suffix = "TGGAATTCTCGGGTGCCAAGG"

infix = "CAAGCAGAAGACGGCATACGAGAT"

suffix_count_adjustment = 21 - cycle_count(PCR_suffix)

barcodes_file = "barcodes27-2-sorted2.txt"
barcodes = []
with open(barcodes_file) as f:
    barcodes = f.read().splitlines()

def strings(alphabet, n):
    if n is 0:
        return [""]
    else:
        ls = strings(alphabet, n-1)
        return [ l + c
                 for l in ls
                 for c in alphabet]

def allKMers(k):
    return strings(NTs, k)

def longest_run(string):
    return max(len(list(l)) for n, l in itertools.groupby(string))

def validKMer(kmer):
    return all([
        1 <= len([c for c in kmer if c in ANTs]) <= 3,
        # no homopolymer runs >3, including with PCR sites
        longest_run(infix[-2:] + kmer + PCR_suffix[:2]) <= 3
    ])

k = 8
validKMers = list(filter(validKMer, allKMers(k)))
validKMers.sort(key = cycle_count, reverse = True)

def generateValidSeqs(barcodes, kmers):
    kmix = 0
    for barcode in barcodes:
        if kmix >= len(kmers):
            break
        kmer = kmers[kmix]

        # make sure barcode works with prefix
        if longest_run(PCR_prefix[-2:] + barcode + infix[:2]) > 3:
            continue

        # make sure barcode works with this kmer
        #if longest_run(barcode[-2:] + kmer[:2]) <= 3:
        kmix += 1
        yield (PCR_prefix + barcode + infix + kmer + PCR_suffix, barcode)

    if kmix < len(kmers) - 1:
        print("DID NOT FIND ENOUGH BARCODES TO PAIR WITH KMERS")

def writeOut(output_file = None):
    #cycle_score_statistics(barcodes, n=4)
    validSeqs = list(generateValidSeqs(barcodes, validKMers))

    seqs = [p[0] for p in validSeqs]
    codes = [p[1] for p in validSeqs]

    # Stop after payload and barcode: 197
    # Stop after payload: 235
    # Stop after barcode: 153
    # No stops: 204


    if len(sys.argv) > 1:
        output_file = sys.argv[1]

    if output_file:
        regions = [ # :: [(Length, CycleCount)]
            # prefix
            ( len(PCR_prefix), 1 ),
            # barcode
            ( 27, 4 ),
            # infix
            ( len(infix), 1 ),
            # payload
            ( k, 5 ),
            # suffix
            ( len(PCR_suffix), 1 ),
        ]
        cycle_score_statistics_stops(regions, seqs)
        #cycle_score_statistics(codes, n=4)

        with open(output_file, 'w') as f:
            for seq in seqs:
                f.write(seq + '\n')
    else:
        for seq in seqs:
            print(seq + '\n')

if __name__ == "__main__":
    writeOut()
