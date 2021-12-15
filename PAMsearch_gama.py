#TODO(Zijin): PAM search refers to (Aho and Corasick, 1975)
#TODO(Zijin): the positive strand(five prime to three prime. Therefore, sgRNA sequence
# is on the upstream of PAM)
#TODO(Zijin): How to deal with N in PAN?
class PAMsearch:
    def __init__(self):
        self.sg_seq = {}
        self.cri_cleave = []
        self.ca_seq = set()

    def revcom(self,s):
        basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        letters = list(s[::-1])
        letters = [basecomp[base] for base in letters]
        return ''.join(letters)


    #build the search event tree
    def F_table(self,str):
        f_return = 0
        goto_graph = list(str)
        f_table = []
        for i in range(len(str)):
            if str[i] in goto_graph[0:i]:
                f_return = str[0:i].rfind(str[i])
                f_table.append(f_return)
            else:
                f_table.append(f_return)
        return f_table

    # get candidate sequences-structure:{cleavage_position:seq}
    def goto_function(self,str,pam,k,gap,start):
        state = 0
        str_num = str.find(pam[state])
        f_table = PAMsearch.F_table(self,pam)
        str_num = str_num + 1
        self.sg_seq = {}
        self.cri_cleave = []
        while(str_num < len(str)):
            if (str[str_num] == pam[state+1]):
                state = state + 1
            else:
                state = f_table[state]
            str_num = str_num + 1
            l = len(pam)
            if (state == l-1):
                if str_num >= l + k:
                    self.sg_seq[str_num-l-gap+start] = str[str_num-l-k:str_num-l].upper()
                    self.cri_cleave.append(str_num - l - gap+start)
                state = 0
        return self.sg_seq

    #reverse and complementary
    def rc_goto_function(self,str,pam,k,gap,start):
        pam = self.revcom(pam)
        state = 0
        str_num = str.find(pam[state])
        f_table = PAMsearch.F_table(self,pam)
        str_num = str_num + 1

        while(str_num < len(str)):
            if (str[str_num] == pam[state+1]):
                state = state + 1
            else:
                state = f_table[state]
            str_num = str_num + 1
            l = len(pam)
            if (state == l-1):
                if str_num+k <= len(str):

                    self.sg_seq[-(str_num-l-gap+start)] = self.revcom(str[str_num:(str_num+k)].upper())
                    self.cri_cleave.append(-(str_num - l - gap+start))
                state = 0
        return self.sg_seq

    #candidate seqs
    def candidate_seq(self,str,pam,k):
        state = 0
        str_num = str.find(pam[state])
        f_table = PAMsearch.F_table(self,pam)
        str_num = str_num + 1
        self.sg_seq = {}
        self.cri_cleave = []
        while(str_num < len(str)):
            if (str[str_num] == pam[state+1]):
                state = state + 1
            else:
                state = f_table[state]
            str_num = str_num + 1
            l = len(pam)
            if (state == l-1):
                if str_num >= l + k:
                    self.ca_seq.add(str[str_num-l-k:str_num-l].upper())


                state = 0
        return self.ca_seq

    def rc_candidate_seq(self,str,pam,k):
        pam = self.revcom(pam)
        state = 0
        str_num = str.find(pam[state])
        f_table = PAMsearch.F_table(self,pam)
        str_num = str_num + 1
        self.sg_seq = {}
        self.cri_cleave = []
        while(str_num < len(str)):
            if (str[str_num] == pam[state+1]):
                state = state + 1
            else:
                state = f_table[state]
            str_num = str_num + 1
            l = len(pam)
            if (state == l-1):

                if str_num+k <= len(str):
                    self.ca_seq.add(self.revcom(str[str_num:(str_num+k)].upper()))

                state = 0
        return self.ca_seq


'''t1 = PAMsearch()
t1.goto_function('ATGCATGAG','ATG',2,1,1)
t1.rc_goto_function('ATGCATGAG','CAT',2,1,1)
print(t1.sg_seq)
t1.candidate_seq('ATGCATGAG','ATG',2)
t1.rc_candidate_seq('ATGCATGAG','ATG',2)
print(t1.ca_seq)
print(t1.revcom('ATGCATGC'))'''