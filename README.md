# codon_optimizer
Codon optimizer is a simple tool to exchange codons that have a low frequency to codons with the highest frequency in the organism of reference like Homo sapiens.

# Parameters
1.  -h, --help            show this help message and exit
2.  -input INPUT          Input path to the fasta file containing sequences to
                        optimize (default: local)
3.  -reffreqtab REFERENCE_FREQUENCY_TABLE
                        Input path to the reference frequency table (default:
                        Homo sapiens)
4.  -output OUTPUT        Output fasta file with codon optimized sequences
                        (default: local)
5.  -nolog                Detailed output (default: False)
6.  --version             show program's version number and exit
