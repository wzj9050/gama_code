# TODO(Zijin): add exceptions when error occur
# import needed class
from CFD_input import CFD_input
from Fasta_reader_gama import Fasta_reader
from GTF_read_gama import GTF_reader
from MergeOn_gama import MergeOn
from PAMsearch_gama import PAMsearch
from Smith_Waterman_Ali_delta import Smi_Wat_Ali
from TSTbuilder import TSTbuilder

# variable parameters
k = 20  # the length of sgRNA(without pam)
k_f = 4  # front matching threhold to filter candidate sequences during building TST
resp_gtf = r'mm10.refGene.gtf'#r'ce11.ensGene.gtf'##r'mm10.refGene.gtf'  # user can use their own gtf
chrom_num = 'chrY'# the chromosome the user interested in
gene_name = ['Gm20931']#['WBGene00014472.1']#

resp_fa = r'mm10.fa'  # user can use their own fasta file, the first column should be the name of chromosomes
pam = b'GG'  # the PAM sequence
gap = 1  # the gap length between slice position
resp_cfd = r'STable 19 FractionActive_dlfc_lookup(1)(1).xlsx'  # user can use their own cds file and extend the limination of 20 input, the example structure is appended
ali_score_dict = {"match": 3, "mismatch": -2, 'gap_opening': -10,'gap_extension': -1}  # affine gap penalty for the Smith-Waterman algorithm,the user can adapt the value in this dictionary
starts = []
ends = []

# read the gtf of annotations(e.x. chrM)
gtf = GTF_reader(resp_gtf)
gtf.GTF_return()
# chro_m_p = gtf.chrom_strand_extract(chrom, '+')
# chro_m_n = gtf.chrom_strand_extract(chrom, '-')
chrom_22_g = gtf.gene_extract(gene_name)
del gtf
for key in chrom_22_g.keys():
    for i in chrom_22_g[key][0]:
        starts.append(int(i[0]))
        ends.append(int(i[1]))
start = min(starts)
end = max(ends)
#print(type(start))


print('gtf read Done!')
# {gene:[[(transcript start, transcript end)], [start_codon],[(exon start, exon end)],[(cds start, cds end)],chrom]}
#print(chrom_22_g)

# read the fasta of mm39
genome = Fasta_reader(resp_fa)
genome.read_generator()
monochr_dict = genome.fa_dict[chrom_num]
seq_baseline = genome.chroms
print('Fasta_read done!')
# search the PAM and return the cleavage point and candidate sg sequences

pam_search = PAMsearch()
for key in monochr_dict.keys():

    if (start >= key[0])&(start < key[1]): 
        seq_gene = monochr_dict[key][start:]
        pam_search.goto_function(seq_gene, pam, k, gap, start)
        pam_search.rc_goto_function(seq_gene, pam, k, gap, start)
        if end <= key[0]: 
            break
print('Target PAM search done!')



candidate_sg_seq = pam_search.sg_seq
candidate_sg = pam_search.cri_cleave
#print(candidate_sg_seq)
i_counter = 0
for chrom in seq_baseline:

    monochr_dict = genome.fa_dict[chrom]
    for seq in monochr_dict.values():
        pam_search.candidate_seq(seq, pam, k)
        pam_search.rc_candidate_seq(seq, pam, k)
candidate_sg_seq_baseline = pam_search.ca_seq
del pam_search#release PAMSEARCH
del genome
print('baseline PAM search done!')
# [cleavage_position]
#print(candidate_sg)
# print(candidate_sg_seq)


# build the TST of sg

tst_builder = TSTbuilder(k, len(candidate_sg))
for key in candidate_sg:
    if not tst_builder.find(candidate_sg_seq[key], k_f):
        tst_builder.insert(candidate_sg_seq[key], key)
candidate_sg_seq = tst_builder.seq_return(0, 0)
for baseline in candidate_sg_seq_baseline:
    if not tst_builder.find(baseline,k_f):
        tst_builder.insert(baseline,0)
candidate_sg_seq_baseline = tst_builder.seq_return(0, 0).values()
del tst_builder
candidate_sg = list(candidate_sg_seq.keys())
candidate_sg.sort()
print('TST filtration done!')
# {cleavage_position: sgRNA_seq}
#print(candidate_sg_seq)


# process alignment
cfd_input = CFD_input(resp_cfd)
cfd_mismatch = cfd_input.mismatch_build()
cfd_deletion = cfd_input.deletion_build()
cfd_insertion = cfd_input.insertion_build()
del cfd_input
candidate_sg_seq_score = {}
for key in candidate_sg:
    print(i_counter)
    i_counter+=1
    smith_waterman = Smi_Wat_Ali(cfd_mismatch, cfd_deletion, cfd_insertion, ali_score_dict)
    score = 0
    for seq in candidate_sg_seq_baseline:
        score += smith_waterman.scoring_matrix(candidate_sg_seq[key], seq)
    candidate_sg_seq_score[key] = [candidate_sg_seq[key], 100 / (100 + score)]
    del smith_waterman
candidate_sg = list(candidate_sg_seq_score.keys())
candidate_sg.sort()
print('Alignment done!')
#print(candidate_sg_seq_score)
# {cleavage:[seq,cfd_score]}


# merge the sequence and gene name and score
merge_on = MergeOn(chrom_22_g, candidate_sg, candidate_sg_seq, candidate_sg_seq_score, chrom_num)
candidate_sg_ongene = merge_on.return_result()
del merge_on
print('Scoring done!')
# {gene:[[cleavage_point, seq, min_distance to start_codon, on_exon?, on_cds?, CFD_score, CG_content, chrom]]}
print(candidate_sg_ongene)
