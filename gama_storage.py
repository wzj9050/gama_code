from Fasta_reader_gama import Fasta_reader

mm39 = Fasta_reader(r'C:\Users\ZijinDesktop2\Desktop\Pyproject\res\mice\mm39.fa')
mm39.read_generator()
data = mm39.fa_dict

with open(r'C:\Users\ZijinDesktop2\Desktop\Pyproject\res\mice\mm39.txt') as f:
    f.write(data)