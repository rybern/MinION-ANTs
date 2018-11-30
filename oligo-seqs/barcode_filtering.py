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

def join_regions(regions, ixs):
    using_regions = [regions[ix] for ix in ixs]

    start = sum(region[0] for region in regions[:ixs[0]])
    end =   start + sum(region[0] for region in using_regions)
    cycle = max(region[1] for region in using_regions)

    return (start, end, cycle)

def cycle_region(seq, region):
    (start, end, cycle) = region
    if cycle == 1:
        return end - start
    else:
        return cycle_count(seq[start : end], cycle)

def cycles_partitions(seq, regions, partition):
    return (cycle_region(seq, join_regions(regions, ixs)) for ixs in partition)

def cycle_score_statistics_stops(seqs):
    regions = [ # :: [(Length, CycleCount)]
        # prefix
        ( 26, 1 ),
        # barcode
        ( 27, 4 ),
        # infix
        ( 24, 1 ),
        # payload
        ( 8, 5 ),
        # suffix
        ( 21, 1 ),
        ]

    strategies = [
        [ [0], [1], [2], [3], [4] ],
        [ [0], [1], [2 ,  3], [4] ],
        [ [0], [1 ,  2], [3], [4] ],
        [ [0], [1 ,  2 ,  3], [4] ],
        [ [0 ,  1], [2], [3], [4] ],
        [ [0 ,  1], [2 ,  3], [4] ],
        [ [0 ,  1 ,  2], [3], [4] ],
        [ [0 ,  1 ,  2 ,  3], [4] ]
    ]

    scores = [[] for _ in strategies]

    i = 0
    for seq in seqs[:100]:
        i += 1
        if(i % 10000 == 0):
            print("done " + str(i) + "\t/" + str(len(seqs)))

        for score, partition in zip(scores, strategies):
            score.append(cycles_partitions(seq, regions, partition))

    maxScores = [list(map(max, zip(*score))) for score in scores]
    totalScores = [sum(score) for score in maxScores]

    header = "  PCR_PREFIX | BARCODE  |  INFIX  | PAYLOAD | PCR_SUFFIX"
    lines = [" [{: 6d}    ] [{: 6d}  ] [{: 6d} ] [{: 6d} ] [ {: 6d}  ]:  ",
             " [{: 6d}    ] [{: 6d}  ] [       {: 6d}    ] [ {: 6d}  ]:  ",
             " [{: 6d}    ] [     {: 6d}       ] [{: 6d} ] [ {: 6d}  ]:  ",
             " [{: 6d}    ] [            {: 6d}          ] [ {: 6d}  ]:  ",
             " [        {: 6d}       ] [{: 6d} ] [{: 6d} ] [ {: 6d}  ]:  ",
             " [        {: 6d}       ] [       {: 6d}    ] [ {: 6d}  ]:  ",
             " [              {: 6d}        ]    [{: 6d} ] [ {: 6d}  ]:  ",
             " [                    {: 6d}               ] [ {: 6d}  ]:  ",
    ]

    print(header)
    for line, total, score in zip(lines, totalScores, maxScores):
        print (line.format(*score) + str(total))

    return

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
