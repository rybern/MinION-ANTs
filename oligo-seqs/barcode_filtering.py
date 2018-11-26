NNTs = "ACGTX"
nntIxs = {nt:ix for (ix, nt) in enumerate(NNTs)}

def cycle_score(a, b, n=5):
    return ((nntIxs[b] - nntIxs[a] - 1) % n) + 1

def cycle_count(seq_, n=5):
    seq = list(seq_)[::-1]
    total = 0
    prev = NNTs[n-1]
    for nt in seq:
        score = cycle_score(prev, nt, n=n)
        #print(score, end=' ')
        total += score
        prev = nt
    #print()
    return total

test = cycle_count("AAAATGCA", 4) == 17

def sort_barcode_file(barcode_out, barcode_in = "barcodes27-2.txt"):
    lines = []
    with open(barcode_in) as f:
        lines = f.read().splitlines()

    lines.sort(key = lambda s: cycle_count(s, n=4))

    with open(barcode_out, "w") as f:
        f.writelines(map(lambda s: s + "\n", lines))

def cycle_score_statistics_stops(seqs):
    prefix_len = 26
    barcode_len = 26
    payload_len = 8
    suffix_len = 21

    payload_barcode_scores = []
    payload_barcode_prefix_scores = []
    barcode_prefix_scores = []
    barcode_scores = []
    payload_scores = []

    print("total seqs: " + str(len(seqs)))
    i = 0

    for seq in seqs:
        i += 1
        if(i % 10000 == 0):
            print("progress " + str(i))

        barcode_prefix = seq[:prefix_len+barcode_len]
        barcode_prefix_scores.append(cycle_count(barcode_prefix, 4))

        payload_barcode_prefix = seq[:prefix_len+barcode_len+payload_len]
        payload_barcode_prefix_scores.append(cycle_count(payload_barcode_prefix, 5))

        payload_barcode = seq[prefix_len:prefix_len+barcode_len+payload_len]
        payload_barcode_scores.append(cycle_count(payload_barcode, 5))

        barcode = seq[prefix_len:prefix_len+barcode_len]
        barcode_scores.append(cycle_count(barcode, 4))

        payload = seq[prefix_len+barcode_len:prefix_len+barcode_len+payload_len]
        payload_scores.append(cycle_count(payload, 5))

    # stop after payload and barcode
    score_all_stop = suffix_len + max(payload_scores) + max(barcode_scores) + prefix_len
    print(suffix_len, max(payload_scores), max(barcode_scores), prefix_len)

    # don't stop, everything at 5
    score_no_stop = [suffix_len + s for s in payload_barcode_prefix_scores]

    # stop after barcode
    score_prefix_stop = suffix_len + max(payload_barcode_scores) + prefix_len

    # stop after payload, then do barcode and prefix at 4
    constant = suffix_len + max(payload_scores)
    print(str(len(barcode_prefix_scores)))
    print(str(barcode_prefix_scores[0]))
    score_payload_stop = [constant + s for s in barcode_prefix_scores]

    print("Stop after payload and barcode: " + str(score_all_stop))
    print("Stop after payload: " + str(max(score_payload_stop)))
    print("Stop after barcode: " + str(score_prefix_stop))
    print("No stops: " + str(max(score_no_stop)))

def cycle_score_statistics(seqs, n=5):
    scores = [cycle_count(seq, n=n) for seq in seqs]
    scores.sort()

    print("Oligo cycle score statistics:")
    stats(scores)

def stats(xs):
    xs.sort()
    print("\tMin:", xs[0])
    print("\tMedian:", xs[int(len(xs) / 2)])
    print("\tMax:", xs[-1])
