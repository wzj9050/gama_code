# space complexity: k*len(GTF)*2
# TODO(Zijin): it is noted that "start" in FASTA starts from 1
# TODO(Zijin): increase the cleavage efficiency: temperaturre(Biopython Tm_staluc function)


# transcript = UTR+
import re


class GTF_reader:
    def __init__(self, resp):  # `resp` is the storage path of GTF
        self.resp = resp
        self.gtf_dict = {}  # self.gtf_dict store the information read in, its structure is shown below:
        # {gene:[[(transcript start, transcript end)], [start_codon],[(exon start, exon end)],[(cds start, cds end)],+or-,chrom]}
        self.gtf_chrom_strand = {}
        self.gene_extract_gtf = {}

    def GTF_return(self):
        gtf_in = open(self.resp)
        while 1:
            line = gtf_in.readline()  # `readline()` is a generator
            # print(line)
            line = line.strip("\n")
            if not line:
                break
            elif not line.startswith("#"):
                line_spl = line.split(sep="\t")
                pattern = re.compile("\".*\"")
                attribute = line_spl[8].split(sep=';')
                # gene_id = pattern.search(attribute[0]).group()[1:-1]
                gene_name = pattern.search(attribute[-2]).group()[1:-1]
                # gene_name = gene_name + '_'+line_spl[2]#add type to the gene_id

                # if the gene haven't been input, initiate its list
                if not gene_name in self.gtf_dict.keys():
                    self.gtf_dict[gene_name] = [[], [], [], [], line_spl[6].encode(), line_spl[0]]

                # add transcripts start and end
                if line_spl[2] == 'transcript':
                    self.gtf_dict[gene_name][0].append((line_spl[3], line_spl[4]))


                # add start codons' start and end
                elif line_spl[2] == 'start_codon':
                    self.gtf_dict[gene_name][1].append(line_spl[4])


                # add exons' start and end
                elif line_spl[2] == 'exon':
                    self.gtf_dict[gene_name][2].append((line_spl[3], line_spl[4]))


                # add CDSs' start and end
                elif line_spl[2] == 'CDS':
                    self.gtf_dict[gene_name][3].append((line_spl[3], line_spl[4]))

        return self.gtf_dict#all strings are bytes

    # return GTF of specified gene
    def gene_extract(self, genes):#`genes` is string list
        self.gene_extract_gtf = {}
        for key in self.gtf_dict.keys():
            # print(key)
            if key in genes:
                self.gene_extract_gtf[key] = self.gtf_dict[key]
        return self.gene_extract_gtf

    # return GTF of specified chromosome
    def chrom_strand_extract(self, chrom, p_or_n):#`chrom` and `p_or_n` are bytes
        for key in self.gtf_dict.keys():
            if (self.gtf_dict[key][5] == chrom) & (self.gtf_dict[key][4] == p_or_n):
                self.gtf_chrom_strand[key] = self.gtf_dict[key]
        return self.gtf_chrom_strand
# Partial test
#gtf = GTF_reader(r'..\mm10.refGene.gtf')
#gtf.GTF_return()
#print(gtf.gene_extract(['Gm20931']))
# print(gtf.gene_extract(['WBGene00010966.1']))
# print((len(gtf.gene),gtf.transcript,gtf.start,gtf.exon,gtf.cds))
# (46904, 61451, 33575, 273640, 225595)
