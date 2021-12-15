'''chrV	20924180
chrX	17718942
chrIV	17493829
chrII	15279421
chrI	15072434
chrIII	13783801
chrM	13794'''

import pandas as pd

# TODO(Zijin): write codes to make sure install additional packages


class Fasta_reader:

    def __init__(self, resp):
        self.seqinfo = ''
        self.seq = ''
        self.resp = resp
        self.fa_dict = {}
        self.num = 0
        self.chroms = []
        self.switch = 0
        self.start = 0
        self.end = 0

    def read_generator(self):
        fasta_f = open(self.resp, 'rb')

        while 1:
            line = fasta_f.readline()
            line = line.rstrip(b'\n')
            # print(line)

            if (line.startswith(b'>') or not line) and self.seqinfo:
                self.chroms.append(self.seqinfo)
                # print(self.switch)
                if self.switch == 1:
                    self.end = self.num - 1
                    self.fa_dict[self.seqinfo][(self.start,self.end)] = bytes(self.seq)

                # self.fa_dict[self.seqinfo] = self.seq
                self.num = 0
                self.switch = 0
            if not line:
                break
            if line.startswith(b'>'):
                self.seqinfo = str(line[1:], encoding="utf-8")
                self.fa_dict[self.seqinfo] = {}
                self.seq = bytearray()


            elif line.startswith(b'N') & (bytes((line[-1],)) == b'N'):
                if self.switch == 1:
                    self.end = self.num - 1
                    self.switch = 0
                    self.fa_dict[self.seqinfo][(self.start,self.end)] = bytes(self.seq)
                self.num += len(line)
                self.seq = bytearray()
            else:
                if self.switch == 0:
                    self.start = self.num
                    self.switch = 1
                self.seq.extend(line)
                self.num += len(line)

    def fasta_gene_extract(self, start, end):
        return self.seq[start:end + 1]



#t1 = Fasta_reader(r'.mm10.fa')
#t1.read_generator()
#print(len(t1.fa_dict.keys()))
#print(t1.fa_dict.keys())
