
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


#print(score(red))
#print(consensus(best), score(best))
#print(hamming_distance("TAAGTT", "TGAATT"))
#print(total_distance("TGCGTT", lecture_dna))




#______________________TASK_2______________________

def count_matrix(motifs):
    l = len(motifs[0])
    counts = {base: [0] * l for base in "ACGT"}
    for motif in motifs:
        for i, base in enumerate(motif):
            counts[base][i] += 1
    return counts


class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.l = len(motifs[0])
        t = len(motifs)
        counts = count_matrix(motifs)
        denom = t + 4 * pseudocount
        self.ppm = {
            base: [(counts[base][i] + pseudocount) / denom for i in range(self.l)]
            for base in "ACGT"
        }

    def lmer_probability(self, lmer):
        prob = 1.0
        for i, base in enumerate(lmer):
            prob *= self.ppm[base][i]
        return prob

    def most_probable_lmer(self, sequence):
        best_lmer = sequence[: self.l]
        max_p = -1.0
        for i in range(len(sequence) - self.l + 1):
            lmer = sequence[i : i + self.l]
            p = self.lmer_probability(lmer)
            if p > max_p:
                max_p = p
                best_lmer = lmer
        return best_lmer

    def consensus(self):
        res = []
        for i in range(self.l):
            best_base = "A"
            best_p = -1.0
            for base in "ACGT":
                if self.ppm[base][i] > best_p:
                    best_p = self.ppm[base][i]
                    best_base = base
            res.append(best_base)
        return "".join(res)


profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
print(profile.consensus())
print(round(profile.lmer_probability("ATGCGTA"), 4))

two = MotifProfile(["GTAC", "TTAA"])
print(two.most_probable_lmer("ACTGGATGACCC"))
print(round(two.lmer_probability("TGAC"), 4))

#________________________________________________________________________
profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
print(profile.l)  # 7
print(profile.ppm["A"])

#________________________________________________________________________

from Bio import motifs
from Bio.Seq import Seq

bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
bio.pseudocounts = 1
print(bio.consensus)     # ATGCGTA
print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]


#______________________TASK_3______________________