# TODO(Zijin): off-target scoring: affine gap penalty: match=3, mismatch=-2, gap opening=-10 and gap extension=-1 (Hiranniramol et al. 2020)
# TODO(Zijin): the match and mismatch scoring can be input as parameters in __init__
# TODO(Zijin): However, I still donot understand why the path_matrix's coordinates should be transposed....
import numpy as np


class Smi_Wat_Ali:

    def __init__(self, CFD_mismatch, CFD_deletion, CFD_insertion, ali_score_dict):
        self.CFD_mismatch = CFD_mismatch
        self.CFD_deletion = CFD_deletion
        self.CFD_insertion = CFD_insertion
        self.score_dict = ali_score_dict
        self.score = 1
        self.DNA_RNA = {b'A'[0]: b'A', b'C'[0]: b'C', b'G'[0]: b'G', b'U'[0]: b'U'}

    def compare(self, a, b):
        if a == b:
            return "match"
        else:
            return "mismatch"

    def max_p(self, l, m, n, c, r):
        max_val = max(l, m, n, 0)
        if l == max_val:
            return [l, [r - 1, c - 1]]
        elif m == max_val:
            return [m, [r - 1, c]]
        elif n == max_val:
            return [n, [r, c - 1]]
        else:
            return [0, [0, 0]]

    # build scoring matrix
    def scoring_matrix(self, seq1, seq2):

        s1 = b" " + seq1[:]
        s2 = b" " + seq2[:]
        self.score = 1
        s_matrix = np.zeros((len(s1), len(s2)))
        path_matrix = np.zeros((len(s1), len(s2), 2))
        exten_countr = 0
        exten_countc = 0
        for r in range(1, len(s1)):
            for c in range(1, len(s2)):
                p1 = s_matrix[r - 1][c - 1] + self.score_dict[self.compare(s1[r], s2[c])]
                p2 = s_matrix[r - 1][c] + self.score_dict['gap_opening'] + self.score_dict[
                    'gap_extension'] * exten_countc
                p3 = s_matrix[r][c - 1] + self.score_dict['gap_opening'] + self.score_dict[
                    'gap_extension'] * exten_countr
                p4 = r
                p5 = c
                s_matrix[r][c] = self.max_p(p1, p2, p3, p4, p5)[0]
                path_matrix[c][r] = self.max_p(p1, p2, p3, p4, p5)[1]
                if s_matrix[r][c] != 0:
                    if (r == path_matrix[r][c][0]) & (c != path_matrix[r][c][1]):
                        exten_countr += 1
                    elif (r != path_matrix[r][c][0]) & (c == path_matrix[r][c][1]):
                        exten_countc += 1

        end = np.argwhere(s_matrix == s_matrix.max())[0]
        self.traceBack(s1, s2, s_matrix, path_matrix, end)
        return self.score

    # traceback for the logest way
    def traceBack(self, s1, s2, s_matrix, path_matrix, end):
        r = int(end[0])
        c = int(end[1])
        while s_matrix[r][c] != 0:

            if (r != path_matrix[r][c][0]) & (c != path_matrix[r][c][1]):

                if s1[r] == s2[c]:
                    self.score *= 1
                else:
                    p = self.CFD_mismatch[bytes((s1[r],s2[c]))][r - 1]
                    self.score *= p

            elif (r != path_matrix[r][c][0]) & (c == path_matrix[r][c][1]):
                p = self.CFD_insertion[bytes((s1[r],))][r - 1]
                self.score *= p

            else:
                p = self.CFD_deletion[bytes((s2[c],))][r - 1]
                self.score *= p

            end = path_matrix[r][c]
            r = int(end[0])
            c = int(end[1])


# t0 = CFD_input(r'C:\Users\ZijinDesktop2\Desktop\Pyproject\res\STable 19 FractionActive_dlfc_lookup.xlsx')
# t01 = t0.mismatch_build()
# t02 = t0.deletion_build()
# t03 = t0.insertion_build()
# t1 = Smi_Wat_Ali(t01,t02,t03)

# t1.scoring_matrix("UAUAUAUGGG","UAUAUAAUGG")
