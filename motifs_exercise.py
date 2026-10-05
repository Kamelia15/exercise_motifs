from Bio import motifs,SeqIO
from Bio.Seq import Seq
import itertools
import random

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
#print(profile.consensus())
#print(round(profile.lmer_probability("ATGCGTA"), 4))

two = MotifProfile(["GTAC", "TTAA"])
#print(two.most_probable_lmer("ACTGGATGACCC"))
#print(round(two.lmer_probability("TGAC"), 4))

#________________________________________________________________________
profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
#print(profile.l)  # 7
#print(profile.ppm["A"])

#________________________________________________________________________

bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
bio.pseudocounts = 1
#print(bio.consensus)     # ATGCGTA
#print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]


#______________________TASK_3______________________

class MotifFinder:
    def __init__(self, sequences, l, seed=None):
        self.sequences = sequences
        self.l = l
        self.rng = random.Random(seed)
        self.windows = [
            [seq[i : i + l] for i in range(len(seq) - l + 1)]
            for seq in sequences
        ]

    def total_distance(self, pattern):
        tot = 0
        for window_list in self.windows:
            min_dist = min(hamming_distance(pattern, w) for w in window_list)
            tot += min_dist
        return tot

    def median_string(self):
        best_pattern = None
        min_dist = float("inf")
        for p_tuple in itertools.product("ACGT", repeat=self.l):
            pattern = "".join(p_tuple)
            dist = self.total_distance(pattern)
            if dist < min_dist:
                min_dist = dist
                best_pattern = pattern
        return best_pattern, min_dist

    def randomized_search(self):
        current_motifs = [self.rng.choice(w_list) for w_list in self.windows]
        best_motifs = current_motifs
        best_score = score(best_motifs)

        while True:
            profile = MotifProfile(current_motifs, pseudocount=1)
            next_motifs = [
                profile.most_probable_lmer(seq) for seq in self.sequences
            ]
            next_score = score(next_motifs)
            if next_score > best_score:
                best_motifs = next_motifs
                best_score = next_score
                current_motifs = next_motifs
            else:
                return best_motifs, best_score

    def best_of(self, runs):
        best_motifs = None
        best_score = -1
        for _ in range(runs):
            m, s = self.randomized_search()
            if s > best_score:
                best_score = s
                best_motifs = m
        return best_motifs, best_score


if __name__ == "__main__":
    finder_lecture = MotifFinder(lecture_dna, 6, seed=42)
    med_pat, med_dist = finder_lecture.median_string()
    best_m, best_s = finder_lecture.best_of(100)
    print(f"Median string: {med_pat} (dist: {med_dist})")
    print(f"best_of(100) consensus: {consensus(best_m)} (score: {best_s})\n")

    planted_seqs = [
        str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")
    ]
    finder_planted = MotifFinder(planted_seqs, 7, seed=42)

    med_pat_p, med_dist_p = finder_planted.median_string()
    best_m_p, best_s_p = finder_planted.best_of(100)
    best_consensus_p = consensus(best_m_p)

    print(f"Planted median string: {med_pat_p} (dist: {med_dist_p})")
    print(f"Planted best_of(100) consensus: {best_consensus_p} (score: {best_s_p})")
    print("Motifs selected per sequence:")
    for i, m in enumerate(best_m_p):
        print(f"  Seq {i+1}: {m}")
    print(f"Do consensus and median string agree? {best_consensus_p == med_pat_p}\n")

    single_runs = 200
    successes = 0
    eval_finder = MotifFinder(planted_seqs, 7, seed=123)
    for _ in range(single_runs):
        m, _ = eval_finder.randomized_search()
        if consensus(m) == med_pat_p:
            successes += 1
    p = successes / single_runs
    print(f"Single run success rate (p): {p:.4f}\n")

    print(f"{'R':<6} | {'Measured Success Rate':<22} | {'Predicted 1 - (1-p)^R':<22}")
    print("-" * 56)
    for R in [5, 10, 20, 50]:
        r_successes = 0
        trials = 50
        for _ in range(trials):
            bm, _ = eval_finder.best_of(R)
            if consensus(bm) == med_pat_p:
                r_successes += 1
        measured = r_successes / trials
        predicted = 1 - (1 - p) ** R
        print(f"{R:<6} | {measured:<22.2f} | {predicted:<22.4f}")
    print()

    found = motifs.create([Seq(m) for m in best_m_p])
    found.pseudocounts = 1
    pssm = found.pssm

    for fpr in [0.01, 0.001]:
        threshold = pssm.distribution().threshold_fpr(fpr)
        total_hits = 0
        found_sites_count = 0
        reverse_hits = 0

        for idx, sequence in enumerate(planted_seqs):
            target_lmer = best_m_p[idx]
            for position, hit_score in pssm.search(
                Seq(sequence), threshold=threshold
            ):
                total_hits += 1
                if position < 0:
                    reverse_hits += 1
                else:
                    lmer_hit = sequence[position : position + 7]
                    if lmer_hit == target_lmer:
                        found_sites_count += 1

        expected_chance = 10 * 54 * 2 * fpr
        print(f"Threshold FPR = {fpr} (Score >= {threshold:.2f}):")
        print(f"  Total hits:             {total_hits}")
        print(f"  Found sites among them: {found_sites_count}")
        print(f"  Reverse strand hits:    {reverse_hits}")
        print(f"  Expected by chance:     ~{expected_chance:.1f}\n")



# 1) The median string is AGATAG with a distance of 2 and best_of(100) returns a consensus of AGATAG with a score of 34.
#
# 2) Find the planted motif. --> YES
# The median string is GCTAAAG and best_of(100) finds the exact same consensus.
#
# 3) How many restarts do you need?
# A single run succeeds only 3 to 4% of the time, but success grows to 18% for R=5, 38% for R=10, 68% for R=20,
# and 88% for R=50, matching the predicted curve. Single runs fail often because random initial picks trap the search
# in local maxima instead of the global optimum.
#
# 4) From motif to scanner
# With threshold_fpr(0.01), you get 19 total hits, including all 10 found sites, 5 hits on the reverse strand,
# and ~11 expected by pure chance, meaning the extra hits are just random background noise. Switching to
# threshold_fpr(0.001) drops total hits down to 10.
#  ----> threshold_fpr(0.001)
