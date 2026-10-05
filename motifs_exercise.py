
#______________________TASK_1______________________

lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]

red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]
best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]


def count_matrix(motifs):
    l = len(motifs[0])
    counts = {base: [0] * l for base in "ACGT"}
    for motif in motifs:
        for i, base in enumerate(motif):
            counts[base][i] += 1
    return counts


def score(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    return sum(max(counts[base][i] for base in "ACGT") for i in range(l))


def consensus(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    res = []
    for i in range(l):
        best_base = "A"
        best_count = -1
        for base in "ACGT":
            if counts[base][i] > best_count:
                best_count = counts[base][i]
                best_base = base
        res.append(best_base)
    return "".join(res)


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def total_distance(pattern, sequences):
    l = len(pattern)
    tot = 0
    for seq in sequences:
        min_dist = min(
            hamming_distance(pattern, seq[i : i + l])
            for i in range(len(seq) - l + 1)
        )
        tot += min_dist
    return tot


print(score(red))
print(consensus(best), score(best))
print(hamming_distance("TAAGTT", "TGAATT"))
print(total_distance("TGCGTT", lecture_dna))




#______________________TASK_2______________________