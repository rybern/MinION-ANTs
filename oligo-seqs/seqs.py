import itertools
import sys
from functools import reduce

NNTs = "GATC"
ANTs = "X"
NTs = NNTs + ANTs

# from Page 20 of
# http://science.sciencemag.org.ezproxy.cul.columbia.edu/content/sci/suppl/2017/03/01/355.6328.950.DC1/Erlich.SM.pdf
PCR_prefix = "GTTCAGAGTTCTACAGTCCGACGATC"
PCR_suffix = "TGGAATTCTCGGGTGCCAAGG"

barcodes_file = "barcodes.txt"
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
        longest_run(kmer + PCR_suffix[:2]) <= 3
    ])

validKMers = list(filter(validKMer, allKMers(8)))

def generateValidSeqs(barcodes, kmers):
    kmix = 0
    for barcode in barcodes:
        if kmix >= len(kmers):
            break
        kmer = kmers[kmix]

        # make sure barcode works with prefix
        if longest_run(PCR_prefix[:2] + barcode[:2]) > 3:
            continue

        # make sure barcode works with this kmer
        if longest_run(barcode[-2:] + kmer[:2]) <= 3:
            kmix += 1
            yield PCR_prefix + barcode + kmer + PCR_suffix

    if kmix < len(kmers) - 1:
        print("DID NOT FIND ENOUGH BARCODES TO PAIR WITH KMERS")

def writeOut():
    validSeqs = generateValidSeqs(barcodes, validKMers)

    if len(sys.argv) > 1:
        output_file = sys.argv[1]
        with open(output_file, 'w') as f:
            for seq in validSeqs:
                f.write(seq + '\n')
    else:
        for seq in validSeqs:
            print(seq + '\n')

if __name__ == "__main__":
    writeOut()
