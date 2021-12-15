import pandas as pd


class CFD_input:
    def __init__(self, resp):
        self.resp = resp
        self.mismatch = pd.read_excel(self.resp, sheet_name=0).values
        self.deletion = pd.read_excel(self.resp, sheet_name=1).values
        self.insertion = pd.read_excel(self.resp, sheet_name=2).values
        self.mismatch_table = {}
        self.deletion_table = {}
        self.insertion_table = {}

    # build mismatch table
    def mismatch_build(self):
        for i in range(len(self.mismatch)):
            transform = (self.mismatch[i][0]).encode()
            # p = self.mismatch_table[i][1]
            score = self.mismatch[i][2]
            if not transform in self.mismatch_table:
                self.mismatch_table[transform] = [score]
            else:
                self.mismatch_table[transform].append(score)
        return self.mismatch_table

    # build deletion table
    def deletion_build(self):
        for i in range(len(self.deletion)):
            delet = (self.deletion[i][0]).encode()
            score = self.deletion[i][2]
            if not delet in self.deletion_table.keys():
                self.deletion_table[delet] = [1]
            self.deletion_table[delet].append(score)
        return self.deletion_table

    # build insertion table
    def insertion_build(self):
        for i in range(len(self.insertion)):
            insert = (self.insertion[i][0]).encode()
            score = self.insertion[i][2]
            if not insert in self.insertion_table.keys():
                self.insertion_table[insert] = [1]
            self.insertion_table[insert].append(score)
        return self.insertion_table

#t1 = CFD_input(r'STable 19 FractionActive_dlfc_lookup(1)(1).xlsx')
#print(t1.insertion_build())
# print(len(t1.deletion_build()['U']))
# print(len(t1.insertion_build()['A']))
