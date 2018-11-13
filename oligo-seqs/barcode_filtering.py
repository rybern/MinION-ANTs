NNTs = "ACGTX"
nntIxs = {nt:ix for (ix, nt) in enumerate(NNTs)}

def cycle_score(a, b):
    return ((nntIxs[b] - nntIxs[a] - 1) % len(NNTs)) + 1

def cycle_count(seq_):
    seq = list(seq_)[::-1]
    total = 0
    prev = 'X'
    for nt in seq:
        score = cycle_score(prev, nt)
        #print(score, end=' ')
        total += score
        prev = nt
    #print()
    return total

def cycle_count_ANT(seq):
    return cycle_count(filter(lambda c: c in NNTs, seq))

test = cycle_count("ACGTAAAA") == 16

def sort_barcode_file(barcode_out, barcode_in = "barcodes_10m.txt"):
    lines = []
    with open(barcode_in) as f:
        lines = f.read().splitlines()

    lines.sort(key = cycle_count)

    with open(barcode_out, "w") as f:
        f.writelines(map(lambda s: s + "\n", lines))

def cycle_score_statistics(seqs, adj = 0):
    scores = [cycle_count_ANT(seq) + adj for seq in seqs]
    scores.sort()

    print("Oligo cycle score statistics:")
    print("\tMin:", scores[0])
    print("\tMedian:", scores[int(len(scores) / 2)])
    print("\tMax:", scores[-1])
