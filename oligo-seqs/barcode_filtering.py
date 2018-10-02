nntIxs = {nt:ix for (ix, nt) in enumerate("ACGT")}

def cycle_score(a, b):
    return ((nntIxs[b] - nntIxs[a] - 1) % 4) + 1

def cycle_count(seq):
    total = 0
    for ix in range(0, len(seq)-1):
        total += cycle_score(seq[ix], seq[ix + 1])
    return total

test = cycle_count("ACGTAAAA") == 16

def sort_barcode_file(barcode_out, barcode_in = "barcodes_10m.txt"):
    lines = []
    with open(barcode_in) as f:
        lines = f.read().splitlines()

    lines.sort(key = cycle_count)

    with open(barcode_out, "w") as f:
        f.writelines(map(lambda s: s + "\n", lines))
