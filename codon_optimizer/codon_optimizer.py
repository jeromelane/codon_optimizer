import argparse
import os
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from numpy import where, arange
from time import asctime, localtime, time
from pandas import read_csv

def program_setting():
    parser = argparse.ArgumentParser(description="Codon optimizer is a simple tool to exchange codons that have a low frequency to codons with the highest frequency in the organism of reference like Homo sapiens ")

    parser.add_argument('-input', action='store', dest='input', required=True,
                        default="",
                        help='input path to the fasta file containing sequences to optimize (default: local)')

    parser.add_argument('-reffreqtab', action='store', dest='reference_frequency_table',
                        default="",
                        help='input path to the reference frequency table, can be found: http://www.kazusa.or.jp/codon/ (default: Homo sapiens)')

    parser.add_argument('-output', action='store', dest='output',
                        default="",
                        help='output fasta file with codon optimized sequences (default: local)')

    parser.add_argument('-nolog', action='store_true', dest='nolog',
                        default=False,
                        help='disable log (default: False)')

    parser.add_argument('--version', action='version', version='codon optimizer 0.2')
    return parser

class Options(object):
    def __init__(self, args):
        self.args = args

    def __str__(self):
        mystring = "Options: \n"
        index = 1
        for k in self.args.__dict__:
            if self.args.__dict__[k] is not None:
                mystring += str(index) + "- " + k + "=" + str(self.args.__dict__[k]) + "\n"
                index += 1
        return mystring

    def __repr__(self):
        return self.format()

class Sequence(Seq):
    def __init__(self, sequence):
        Seq.__init__(self, sequence.upper())
        self.sequence = str(self)
        self.qc_msg = []

    def get_optimized_codon(self,codon, codon_frequency_table):
        codon_indexes = where(codon_frequency_table['Codon'] == codon.upper())[0]
        codon_row = codon_frequency_table.iloc[codon_indexes]

        amino_acid = codon_row['AmAcid'].iloc[0]

        amino_acid_indexes = where(codon_frequency_table['AmAcid'] == amino_acid)[0]
        amino_acid_rows = codon_frequency_table.ix[amino_acid_indexes]
        optimized_codon_index = where(amino_acid_rows["/1000"] == max(amino_acid_rows["/1000"]))[0]

        optimized_codon = amino_acid_rows["Codon"].iloc[optimized_codon_index].iloc[0]
        return optimized_codon

    def get_codon(self, start):
        return Sequence(self.sequence[start:start + 3])

    def get_codon_optimized_sequence(self, codon_frequency_table):
        optimized_codons = []
        for index in arange(0, len(self.sequence), 3):
            codon = self.get_codon(index)
            if codon.isx3():
                optimized_codons.append(self.get_optimized_codon(codon, codon_frequency_table))
            else:
                self.qc_msg.append("Sequence " + codon + " not multiple of 3")
                optimized_codons.append(codon)
        return Sequence("".join(optimized_codons))

    def isLongerThan(self, ln):
        return len(self.sequence) > ln

    def isStopCodonAtIts3prime(self):
        lngth = len(self.sequence)
        last_codon = self.sequence[(lngth-3):lngth]
        if last_codon not in ["TGA", "TAG", "TAA"]:
            self.qc_msg.append("No Stop codon at its 3 prime: " + last_codon)
            return False
        return True

    def isMetCodonAtIts5prime(self):
        start_codon = self.sequence[0:3]
        if start_codon not in ["ATG"]:
            self.qc_msg.append("No Methionine codon at its 5 prime: " + start_codon)
            return False
        return True

    def isNucleotide(self):
        if self.sequence == "":
            return False
        for nt in self.sequence:
            if nt not in ["A", "T", "G", "C", "U"]:
                self.qc_msg.append("Not nucleotide: " + nt)
                return False
        return True

    def isDNA(self):
        for nt in self.sequence:
            if not nt in ["A","T","G","C"]:
                if nt == "U":
                    return False
                else:
                    self.qc_msg.append("Not nucleotide: " + nt)
        return True

    def isRNA(self):
        for nt in self.sequence:
            if nt == "U":
                return True
        return False

    def convertRNAtoDNA(self, inplace=True):
        converted_seq = self.sequence.replace("U", "T")
        if not inplace:
            return converted_seq
        self.sequence = converted_seq

    def isStopCodonInCodingFrame(self,coding_frame=0):
        sequence = self.sequence[coding_frame:]
        for pos in arange(0, len(sequence)-3, 3):
            mySeqToCompare = sequence[pos:(pos + 3)]
            if Sequence(mySeqToCompare).isLongerThan(2):
                if mySeqToCompare in ["TGA","TAG","TAA"]:
                    return True
        return False

    def general_info(self):
        self.isStopCodonAtIts3prime()
        self.isMetCodonAtIts5prime()

    def isValid(self):
        if not self.isNucleotide():
            self.qc_msg.append("Sequence has unrecognized character.")
            return False
        if not self.isLongerThan(2):
            self.qc_msg.append("Sequence is too short.")
            return False
        if self.isStopCodonInCodingFrame():
            self.qc_msg.append("Sequence contains stop codon.")
            return False
        return True

    def isx3(self):
        return (len(self.sequence) % 3) == 0

def codon_optimizer():
    parser = program_setting()
    args = parser.parse_args()

    execution_path = os.path.dirname(__file__)
    os.chdir(execution_path)

    output = args.output
    if not os.path.exists(args.input):
        raise IOError("Please provide a fasta file containing sequences to analyze: file path=" + args.input)

    if not args.nolog:
        logging.basicConfig(filename=os.path.dirname(output) + '/codon_optimizer.log', level=logging.INFO)

    logging.info("--- Program starts at: %s ---" % asctime(localtime(time())))
    logging.info("Codon optimizer")
    logging.info(Options(args))

    not_optimized = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

    if args.reference_frequency_table != "":
        if not os.path.exists(args.reference_frequency_table):
            raise IOError("Please provide a correct path to the reference file.")
        codon_frequency_table = read_csv(args.reference_frequency_table)
    else:
        codon_frequency_table = read_csv(execution_path + "/reference_seq/Homo_sapiens_codon_frequency.csv")

    valid_sequence_ids=[]
    optimized_sequences = []
    for seq_id in not_optimized.keys():
        seqToTest = Sequence(str(not_optimized[seq_id].seq))
        if seqToTest.isRNA():
            seqToTest.convertRNAtoDNA()
        if seqToTest.isValid():
            seqToTest.general_info()
            optimized_seq = seqToTest.get_codon_optimized_sequence(codon_frequency_table)
            seqrec = SeqRecord(optimized_seq, id=seq_id, description="codon optimized sequence")
            optimized_sequences.append(seqrec)
            valid_sequence_ids.append(seq_id)

        logging.info(seq_id + " : " + str(seqToTest.qc_msg))

    if len(valid_sequence_ids) == 0:
        print("No valid sequences. No sequence optimized.")
    elif output != "":
        with open(output, 'w') as f_out:
            SeqIO.write(optimized_sequences, f_out, 'fasta')
    else:
        for sequence in optimized_sequences:
            print(sequence.format('fasta'))

    logging.info("--- Program ends at: %s ---" % asctime(localtime(time())))

if __name__ == '__main__':
    try:
        codon_optimizer()
    except IOError as e:
        print(e.args)
