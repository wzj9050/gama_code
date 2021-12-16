# TODO(Zijin): transform base letters into four-bite codes
# TODO(Zijin): return leaf path(Done); return sequence(Done); return suffix; dynamic size of TST
# TODO(Zijin): apply prefix DAWG to this algorithm to collapse space
import numpy as np
class TSTbuilder:

    def __init__(self, k, n):

        self.nex = np.zeros((4**int(3*k/5),4))

        self.cnt = 0
        self.exist = [False] * 4
        self.exist_f = [False] * 4
        self.subleaf = []
        self.seqs = {}
        self.only_seqs = []
        for i in range(50):
            exec('self.path_l{} = []'.format(i))

    num_to_base = {0.0: b'A', 1.0: b'U', 2.0: b'G', 3.0: b'C'}
    base_to_num = {0: 0, 19: 1, 6: 2, 2: 3, 13:0}# 13 is to transform 'N' into 'A'

    # insert seq meanwhile removing all duplicate seqs
    def insert(self, s, id):
        s = s.upper()
        p = 0
        for i in range(0, len(s)):
            c = self.base_to_num[s[i] - b'A'[0]]

            if self.nex[p][c] == 0:
                self.cnt += 1
                self.nex[p][c] = self.cnt

                self.exist_f.append(False)
                self.exist.append(False)
            p = self.nex[p][c]
            p = int(p)
        self.exist[p] = 'True' + " " + str(id)
        self.exist_f[p] = True
        
    # find if seq exist in tree and mark duplicate seqs and seqs with excessive overlap
    def find(self, s, k_f):
        s = s.upper()
        p = 0
        f_count = 0
        for i in range(0, len(s)):
            c = self.base_to_num[s[i] - b'A'[0]]
            if (self.nex[p][c] == 0) & (f_count < k_f):
                return False
            elif (self.nex[p][c] == 0) & (f_count >= k_f):
                return True
            p = self.nex[p][c]
            p = int(p)
            f_count += 1
        self.exist[p] = 'True' + " " + str(-1)
        return self.exist[p].split(sep=" ")[0]

    # return the built TNT
    def tst_return(self):
        return self.nex

    # return exist_table
    def exist_return(self):
        return self.exist

    # judege the seqs with excessive overlap
    def tst_filter(self):
        exist_tem = self.exist_f[:]
        a = 0
        b = 0
        while (a < len(self.exist_f) - 1) & (b < len(self.exist_f) - 1):
            if b == 0:
                a = self.exist_f.index(True)
            elif 1 in self.exist_f[b:]:
                a = b
            else:
                break
            if 1 in self.exist_f[a + 1:]:
                b = self.exist_f.index(True, a + 1, len(self.exist_f))
            else:
                break
                # print()
            if (b - a) < 3:
                self.exist_f[a] = False
                self.exist_f[b] = False
        for i in range(len(self.exist_f)):
            if self.exist_f[i] == False:
                self.exist[i] == False

    # traverse TST and return their original seqs and marks
    def seq_return(self, layer, n):
        if sum(self.nex[layer]) == 0:
            if layer == 0:
                self.seqs[int(self.exist[layer].split(sep=" ")[1])] = ""
            else:
                exec('self.seqs[int(self.exist[layer].split(sep=" ")[1])] = self.path_l{}'.format(n - 1))

        else:
            for node in range(4):
                if self.nex[layer][node] != 0:
                    if n == 0:
                        self.path_l0 = self.num_to_base[node]
                    else:
                        exec('self.path_l{} = self.path_l{}[:]'.format(n, n - 1))
                        exec('self.path_l{} += self.num_to_base[node]'.format(n))
                    sublayer = int(self.nex[layer][node])
                    TSTbuilder.seq_return(self, sublayer, n + 1)
        self.seqs[-1] = ""
        del self.seqs[-1]
        return self.seqs





'''test1 = TSTbuilder(5,20)
test2 = TSTbuilder(5,20)
test1.insert(b"TTTGC",1)
test1.insert(b"AAAGC",2)
if not test1.find(b"AAAGC",3):
    test1.insert(b"AAAGC", 2)
test1.insert(b"ATTGG",0)
if not test1.find(b'AAAGG',-3):
    test1.insert(b'AAAGG',4)
test1.insert(b"AAAGC",-2)

test1.insert(b"AATGC",5)
if not test1.find(b"AATGC",-3):
    test1.insert(b"AATGC", -2)

print(test1.find(b"AAAGC",3))
test1.tst_filter()
print(test1.seq_return(0,0))
print(test1.exist)'''
