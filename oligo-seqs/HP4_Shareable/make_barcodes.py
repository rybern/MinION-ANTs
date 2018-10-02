"""
Poorly coded by Eli Lyons
Copyright Tupac Bio, Inc. 2016

This script outputs a barcode library of user defined size N.
The batch code file and target code file provided should be in the same directory as this script.

Example usage:
        If python 3 is your default pyton and N = 1000:
        $ python make_barcodes.py 1000
        else:
        $ python3 make_barcodes.py 1000

 where 1000 is the desired number of barcodes
"""
import sys
import csv


no_barcodes = int(sys.argv[1])

batch_codes = []
target_codes = []
no_of_gen_codes = 0

with open('batch_codes.txt', newline="") as csvfile:
    seq_reader = csv.reader(csvfile, delimiter=",", skipinitialspace=True, quotechar="|")
    for line in seq_reader:
        batch_codes.append(line[0])

with open('target_codes.txt', newline="") as csvfile:
    seq_reader = csv.reader(csvfile, delimiter=",", skipinitialspace=True, quotechar="|")
    for line in seq_reader:
        target_codes.append(line[0])



output_file = open('barcodes.txt', 'w')
# can check here if it possible to generate that many barcodes

for batch_code in batch_codes:
    for target_code in target_codes:
        barcode = batch_code + target_code
        output_file.write(barcode)
        output_file.write("\n")
        no_of_gen_codes += 1
        if no_of_gen_codes == no_barcodes:
            break
    if no_of_gen_codes == no_barcodes:
            break

output_file.close()






