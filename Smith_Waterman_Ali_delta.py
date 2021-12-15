# TODO(Zijin): off-target scoring: affine gap penalty: match=3, mismatch=-2, gap opening=-10 and gap extension=-1 (Hiranniramol et al. 2020)
# TODO(Zijin): the match and mismatch scoring can be input as parameters in __init__
# TODO(Zijin): However, I still donot understand why the path_matrix's coordinates should be transposed....
import numpy as np
from CFD_input import CFD_input

class Smi_Wat_Ali:

    def __init__(self, CFD_mismatch, CFD_deletion, CFD_insertion, ali_score_dict):
        self.CFD_mismatch = CFD_mismatch
        self.CFD_deletion = CFD_deletion
        self.CFD_insertion = CFD_insertion
        self.score_dict = ali_score_dict
        self.score = 1
        self.max_score = 0
        self.max_position = (0,0)
        self.r1 = b""
        self.r2 = b""
    def compare(self, a, b):
        if a == b:
            return "match"
        else:
            return "mismatch"

    def max_p(self, l, m, n, c, r):
        max_val = max(l, m, n, 0)
        if 0 == max_val:
            return [0, 0]
        elif l == max_val:
            return [l, 1]
            #return [l, [r - 1, c - 1]]
        elif m == max_val:
            return [m, 2]
            #return [m, [r - 1, c]]
        elif n == max_val:
            return [n, 3]
            #return [n, [r, c - 1]]



    # build scoring matrix
    def scoring_matrix(self, seq1, seq2):

        s1 = b" " + seq1[:]
        s2 = b" " + seq2[:]
        self.score = 1

        path_matrix = np.zeros((len(s1), len(s2), 2))
        exten_countr = 0
        exten_countc = 0
        for r in range(1, len(s1)):
            for c in range(1, len(s2)):
                p1 = path_matrix[r - 1][c - 1][0] + self.score_dict[self.compare(s1[r], s2[c])]
                p2 = path_matrix[r - 1][c][0] + self.score_dict['gap_opening'] + self.score_dict[
                    'gap_extension'] * exten_countc
                p3 = path_matrix[r][c - 1][0] + self.score_dict['gap_opening'] + self.score_dict[
                    'gap_extension'] * exten_countr
                p4 = r
                p5 = c


                path_matrix[c][r] = self.max_p(p1, p2, p3, p4, p5)
                path_trans = path_matrix[r][c][1]
                path_score = path_matrix[r][c][0]
                if path_score > self.max_score:
                    self.max_score = path_score
                    self.max_position = (r,c)
                if path_score != 0:
                    if path_trans == 3:
                        exten_countr += 1
                    elif path_trans == 2:
                        exten_countc += 1

        end = self.max_position
        self.traceBack(s1, s2, path_matrix, end)
        return self.score

    # traceback for the logest way
    def traceBack(self, s1, s2, path_matrix, end):
        r = int(end[0])
        c = int(end[1])
        path_score = path_matrix[c][r][0]
        path_trans = path_matrix[c][r][1]
        while path_score != 0:

            if path_trans == 1:
                #print((r,c))

                #print(s1[r])
                if s1[r] == s2[c]:
                    self.score *= 1
                    #self.r1 = bytes((s1[r],))+self.r1
                    #self.r2 = bytes((s2[c],))+ self.r2
                else:
                    p = self.CFD_mismatch[bytes((s1[r],s2[c]))][r - 1]
                    self.score *= p
                r = r-1
                c = c-1

            elif path_score == 2:
                p = self.CFD_insertion[bytes((s1[r],))][r - 1]
                self.score *= p
                r = r-1
                c = c
                #self.r1 = bytes((s1[r],))+ self.r1
                #self.r2 = bytes((s2[c], ))+b'-'
            else:
                p = self.CFD_deletion[bytes((s2[c],))][r - 1]
                self.score *= p
                r = r
                c = c-1
                #self.r1 = bytes((s1[r], ))+b'-'
                #self.r2 = bytes((s2[c], ))+self.r2
            path_score = path_matrix[c][r][0]


ali_score_dict = {"match": 3, "mismatch": -2, 'gap_opening': -10,'gap_extension': -1}
t0 = CFD_input(r'STable 19 FractionActive_dlfc_lookup(1)(1).xlsx')
t01 = t0.mismatch_build()
t02 = t0.deletion_build()
t03 = t0.insertion_build()
t1 = Smi_Wat_Ali(t01,t02,t03,ali_score_dict)

'''t1.scoring_matrix(b"UAUAUAUGGGAGAGAGAGAGAGA",b"UAUAUAUGGAGAGAGAGAGAGAG")
print(t1.r1)
print(t1.r2)'''
