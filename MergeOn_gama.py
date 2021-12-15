# {gene:[[(transcript start, transcript end)], [start_codon],[(exon start, exon end)],[(cds start, cds end)],chrom]}
# {gene:[[cleavage_point, seq, min_distance to start_codon, on_exon?, on_cds?, CFD_score, CG_content, chrom]]}

class MergeOn:
    def __init__(self, gtf, slice, slice_sq, slice_sq_score, chromo):
        self.gtf = gtf

        self.slice = slice
        self.slice_sq = slice_sq
        self.result = {}
        self.slice_sq_score = slice_sq_score
        self.chromo = chromo

    # merge all information from fasta, gtf, cdf, and alignment
    def return_result(self):
        break_start = 0
        break_end = 0
        for i in range(len(self.slice)):
            # du_count = 0
            bp = int(self.slice[i])
            if bp < 0:
                p_or_n = b'-'

            else:
                p_or_n = b'+'

            for key in self.gtf.keys():
                on_gene = False
                g_info = self.gtf[key]
                for s in g_info[0]:
                    start = int(s[0])
                    end = int(s[1])
                    if (bp >= start) & (bp <= end):
                        on_gene = True
                if on_gene == True:

                    if not key in self.result.keys():
                        self.result[key] = []

                    if len(g_info[1]) > 0:
                        du = []
                        for j in g_info[1]:
                            if (bp - int(j)) > -2:
                                du.append(bp - int(j))
                        if len(du) > 0:
                            du = min(du)
                        else:
                            du = b'-'.decode()
                    else:
                        du = b'-'.decode()

                    on_exon = b'-'.decode()
                    if len(g_info[2]) > 0:
                        on_exon = False
                        for k in g_info[2]:
                            if (bp >= int(k[0])) & (bp <= int(k[1])):
                                on_exon = True

                    on_cds = b'-'.decode()
                    if len(g_info[3]) > 0:
                        on_cds = False
                        for l in g_info[3]:
                            if (bp >= int(l[0])) & (bp <= int(l[1])):
                                on_cds = True
                    seq = self.slice_sq[bp]
                    cg_content = (seq.count(b'C') + seq.count(b'G')) / len(seq)
                    self.result[key].append(
                        [bp, seq.decode(), du, on_exon, on_cds, self.slice_sq_score[bp][1], cg_content, p_or_n.decode(), self.chromo])

                elif bp < start:
                    break

        return self.result
