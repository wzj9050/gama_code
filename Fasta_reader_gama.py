'''chrV	20924180
chrX	17718942
chrIV	17493829
chrII	15279421
chrI	15072434
chrIII	13783801
chrM	13794'''
#TODO(Zijin): write codes to make sure install additional packages


class Fasta_reader:

    def __init__(self,resp):
        self.seqinfo= ''
        self.seq = ''
        self.resp = resp
        self.fa_dict = {}
        self.num = 0
        self.chroms = []
        self.switch = 0
        self.start = 0
        self.end = 0
    def read_generator(self):
        fasta_f = open(self.resp)

        while 1:
            line = fasta_f.readline()
            line = line.strip('\n')
            #print(line)

            if (line.startswith('>') or not line) and self.seqinfo:
                self.chroms.append(self.seqinfo)
               # print(self.switch)
                if self.switch == 1:
                    self.end = self.num - 1
                    exec('self.{}[(self.start,self.end)]=self.seq'.format(self.seqinfo))
                #self.fa_dict[self.seqinfo] = self.seq
                self.num = 0
                self.switch = 0
            if not line:
                break
            if line.startswith('>'):
                self.seqinfo = line[1:]
                self.seq = ""
                exec('self.{} = dict()'.format(self.seqinfo))

            elif line.startswith('N')& (line[-1]=='N'):
                if self.switch ==1:
                    self.end = self.num - 1
                    self.switch = 0
                    exec('self.{}[(self.start,self.end)]=self.seq'.format(self.seqinfo))
                self.num += len(line)
                self.seq = ''
            else:
                if self.switch ==0:
                    self.start = self.num
                    self.switch = 1
                    print(self.start)
                self.seq = u'{}{}'.format(self.seq,line)
                self.num += len(line)



    def fasta_gene_extract(self,start,end):
        return self.seq[start:end+1]

if __name__ == "__main__":
    t1 = Fasta_reader(r'..\chrM.fa')
    t1.read_generator()
    pass
    print(t1.chrM)

