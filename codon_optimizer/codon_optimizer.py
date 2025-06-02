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
    parser = argparse.ArgumentParser(description="Codon optimizer is a simple tool to exchange codons that have a low "
                                                 "frequency to codons with the highest frequency in the organism of "
                                                 "reference like Homo sapiens ")

    parser.add_argument('-input', action='store', dest='input', required=True,
                        default="",
                        help='input path to the fasta file containing sequences to optimize (default: local)')

    parser.add_argument('-reffreqtab', action='store', dest='reference_frequency_table',
                        default="",
                        help='input path to the reference frequency table, can be found: http://www.kazusa.or.jp/codon/'
                             ' (default: Homo sapiens)')

    parser.add_argument('-output', action='store', dest='output',
                        default="",
                        help='output fasta file with codon optimized sequences (default: local)')

    parser.add_argument('-nolog', action='store_true', dest='nolog',
                        default=False,
                        help='disable log (default: False)')

    parser.add_argument('--version', action='version', version='codon optimizer 0.3')
    return parser


class Sequence(Seq):
    def __init__(self, sequence):
        Seq.__init__(self, sequence.upper())
        self.sequence = str(self)
        self.qc_msg = []

    def get_optimized_codon(self, codon, codon_frequency_table):
        """
        Retrieve most frequent codon corresponding to input codon amino acid from the input codon frequency table
        :param codon: codon
        :param codon_frequency_table: codon frequency table
        :return: optimized codon
        """
        codon_indexes = where(codon_frequency_table['Codon'] == codon.upper())[0]
        codon_row = codon_frequency_table.iloc[codon_indexes]

        amino_acid = codon_row['AmAcid'].iloc[0]

        amino_acid_indexes = where(codon_frequency_table['AmAcid'] == amino_acid)[0]
        amino_acid_rows = codon_frequency_table.iloc[amino_acid_indexes]  # Fix: Use `.iloc[]` instead of `.ix[]`
        optimized_codon_index = where(amino_acid_rows["/1000"] == max(amino_acid_rows["/1000"]))[0]

        return amino_acid_rows["Codon"].iloc[optimized_codon_index].iloc[0]

    def get_codon(self, start):
        return Sequence(self.sequence[start:start + 3])

    def get_codon_optimized_sequence(self, codon_frequency_table):
        """
        Replace each of codon input sequence with its most frequently represented triplet sequence
        :param codon_frequency_table: csv formatted codon frequency table
        :return: codon optimized sequence
        """
        optimized_codons = []
        for index in arange(0, len(self.sequence), 3):
            codon = self.get_codon(index)
            if codon.isx3():
                optimized_codons.append(self.get_optimized_codon(codon, codon_frequency_table))
            else:
                self.qc_msg.append("Sequence " + codon + " not multiple of 3")
                optimized_codons.append(codon)
        return Sequence("".join(optimized_codons))

    def is_longer_than(self, ln):
        return len(self.sequence) > ln

    def is_stop_codon_at_its_3prime(self):
        lngth = len(self.sequence)
        last_codon = self.sequence[(lngth-3):lngth]
        if last_codon not in ["TGA", "TAG", "TAA"]:
            self.qc_msg.append("No Stop codon at its 3 prime: " + last_codon)
            return False
        return True

    def is_met_codon_at_its_5prime(self):
        start_codon = self.sequence[0:3]
        if start_codon not in ["ATG"]:
            self.qc_msg.append("No Methionine codon at its 5 prime: " + start_codon)
            return False
        return True

    def is_nucleotide(self):
        if self.sequence == "":
            return False
        for nt in self.sequence:
            if nt not in ["A", "T", "G", "C", "U"]:
                self.qc_msg.append("Not nucleotide: " + nt)
                return False
        return True

    def is_dna(self):
        for nt in self.sequence:
            if nt == "U":
                return False
            if nt not in ["A", "T", "G", "C"]:
                self.qc_msg.append("Not nucleotide: " + nt)
        return True

    def is_rna(self):
        for nt in self.sequence:
            if nt == "U":
                return True
        return False

    def convert_rna_to_dna(self, inplace=True):
        converted_seq = self.sequence.replace("U", "T")
        if not inplace:
            return converted_seq
        self.sequence = converted_seq

    def is_stop_codon_in_coding_frame(self, coding_frame=0):
        sequence = self.sequence[coding_frame:]
        for pos in arange(0, len(sequence)-3, 3):
            my_seq_to_compare = sequence[pos:(pos + 3)]
            if Sequence(my_seq_to_compare).is_longer_than(2) and (my_seq_to_compare in ["TGA", "TAG", "TAA"]):
                return True
        return False

    def general_info(self):
        self.is_stop_codon_at_its_3prime()
        self.is_met_codon_at_its_5prime()

    def is_valid(self):
        """
        Verify if sequence is valid for processing
        :return: True if sequence can be processed, False if not
        """
        if not self.is_nucleotide():
            self.qc_msg.append("Sequence has unrecognized character.")
            return False
        if not self.is_longer_than(2):
            self.qc_msg.append("Sequence is too short.")
            return False
        if self.is_stop_codon_in_coding_frame():
            self.qc_msg.append("Sequence contains stop codon.")
            return False
        return True

    def isx3(self):
        return (len(self.sequence) % 3) == 0


def check_path(args):
    """
    Check if input file path exists
    :param args: input data to codon optimizer
    """
    if not os.path.exists(args.input):
        logging.error(f"Input file does not exist: {args.input}")
        raise IOError(f"Please provide a fasta file containing sequences to analyze: file path={args.input}")


def get_reference_frequency_table(args, execution_path):
    """
    Verify existence of a reference codon frequency table and load it in CSV format
    :param args: codon optimizer input data
    :param execution_path: path where program is executed
    :return: reference codon frequency table in CSV format
    """
    if args.reference_frequency_table != "":
        if not os.path.exists(args.reference_frequency_table):
            logging.error(f"Reference frequency table does not exist: {args.reference_frequency_table}")
            raise IOError("Please provide a correct path to the reference file.")
        return read_csv(args.reference_frequency_table)
    else:
        logging.info("No reference frequency table provided, using Homo sapiens as default table.")
        logging.info(execution_path + "/reference_seq/Homo_sapiens_codon_frequency.csv")
        return read_csv(execution_path + "/reference_seq/Homo_sapiens_codon_frequency.csv")


def process_sequences(not_optimized, codon_frequency_table):
    """
    Optimize sequence codons:
        - Converts it to RNA if input sequence is DNA
        - Check that sequence can be optimizable
        - Replace codons by there most frequent synonymous codon
    :param not_optimized: Sequence to optimize
    :param codon_frequency_table: codon frequency table
    :return: optimized sequence and their ids
    """
    valid_sequence_ids = []
    optimized_sequences = []
    for seq_id in not_optimized.keys():
        seq_to_test = Sequence(str(not_optimized[seq_id].seq))
        if seq_to_test.is_rna():
            seq_to_test.convert_rna_to_dna()
        if seq_to_test.is_valid():
            seq_to_test.general_info()
            optimized_seq = seq_to_test.get_codon_optimized_sequence(codon_frequency_table)
            seqrec = SeqRecord(optimized_seq, id=seq_id, description="codon optimized sequence")
            optimized_sequences.append(seqrec)
            valid_sequence_ids.append(seq_id)
        if not seq_to_test.is_valid():
            logging.warning(f"{seq_id} : {str(seq_to_test.qc_msg)}")
        else:
            logging.info(f"{seq_id} : Sequence is valid and optimized.")
    return valid_sequence_ids, optimized_sequences


def codon_optimizer():
    """
    Optimize sequence codons according to reference frequency codon table and format input sequence in fasta format
    :return: fasta formatted sequences
    """
    parser = program_setting()
    args = parser.parse_args()

    execution_path = os.path.dirname(__file__)
    os.chdir(execution_path)

    output = args.output
    check_path(args)

    # Ensure output path is correct (set it to the current directory if not specified)
    if not output:
        output = os.path.join(os.getcwd(), 'default_output.fasta')

    # Set up logging if not disabled
    if not args.nolog:
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)

        # StreamHandler for logging to terminal
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

        # FileHandler for logging to file
        log_filename = os.path.join(os.getcwd(), 'codon_optimizer.log')
        file_handler = logging.FileHandler(log_filename)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logging.info("--- Program starts at: %s ---" % asctime(localtime(time())))
    logging.info("Codon optimizer")
    logging.info(f"Arguments: {args}")

    codon_frequency_table = get_reference_frequency_table(args, execution_path)
    not_optimized = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    valid_sequence_ids, optimized_sequences = process_sequences(not_optimized, codon_frequency_table)

    if len(valid_sequence_ids) == 0:
        logging.warning("No valid sequences. No sequence optimized.")
    elif output != "":
        with open(output, 'w') as f_out:
            SeqIO.write(optimized_sequences, f_out, 'fasta')
        logging.info(f"Optimized sequences written to {output}")
    else:
        for sequence in optimized_sequences:
            logging.info(sequence.format('fasta'))

    logging.info("--- Program ends at: %s ---" % asctime(localtime(time())))


if __name__ == '__main__':
    try:
        codon_optimizer()
    except IOError as e:
        logging.error(f"IOError: {e.args}")
