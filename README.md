# codon_optimizer
A simple tool to produce DNA sequences with the most abundant codons in a specific organism.

# Parameters
*  -h, --help            show this help message and exit
*  -input FILENAME          Input path to the fasta file containing sequences to
                        optimize (default: local)
*  -reffreqtab FILENAME
                        Input path to the reference frequency table should be represented in a specific format (see reference_seq folder) (default:
                        Homo sapiens)
*  -output FILENAME        Output fasta file with codon optimized sequences
                        (default: local)
*  -nolog                Detailed output (default: False)
*  --version             show program's version number and exit

 # Usage
 
 ```python
# Most simple execution
 python codon_optimizer.py -input myfile.fa
# disable log
python codon_optimizer.py -input myfile.fa -nolog
# use a specific codon usage table
python codon_optimizer.py -input myfile.fa -reffreqtab codonusage.txt
  ```

