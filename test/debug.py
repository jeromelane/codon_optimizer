import unittest
from codon_optimizer.codon_optimizer import codon_optimizer, Sequence
import filecmp
import sys
import os


class TestCodonOptimizer(unittest.TestCase):
    MIX_CODING_CASES = "/cases/mix_coding_not_coding_optimized_formatted.fasta"
    CODING_CASES = "/cases/input_mix_coding_not_coding.fasta"
    OUTPUT_FOR_CODING_CASES = "/output/mix_coding_not_coding_optimized.fasta"

    def test_input_output_coding_sequence(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + "/cases/coding_seq_codon_optimized_formatted.fasta"
        input_fasta = execution_path + "/cases/input_coding_sequences.fa"
        output = execution_path + "/output/coding_seq_codon_optimized.fasta"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()
        a = open(reference_output_fasta, 'r')
        b = open(output, 'r')
        a_read = a.read()
        b_read = b.read()
        a.close()
        b.close()
        self.assertTrue(a_read == b_read)

    def test_input_output_not_coding_sequence(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        input_fasta = execution_path + "/cases/input_not_coding_sequences.fa"
        output = execution_path + "/output/not_coding_seq_codon_optimized.fasta"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        self.assertFalse(os.path.exists(output))

    def test_input_output_mix_coding_not_coding_sequence(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + self.MIX_CODING_CASES
        input_fasta = execution_path + self.CODING_CASES
        output = execution_path + self.OUTPUT_FOR_CODING_CASES

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        a = open(reference_output_fasta, 'r')
        b = open(output, 'r')
        a_read = a.read()
        b_read = b.read()
        a.close()
        b.close()
        self.assertTrue(a_read == b_read)

    def test_input_output_mix_coding_not_coding_sequence_different_ref_freqtab(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + self.MIX_CODING_CASES
        input_fasta = execution_path + self.CODING_CASES
        output = execution_path + self.OUTPUT_FOR_CODING_CASES
        reffreqtab = execution_path + "/cases/Mus_musculus_codon_frequency.csv"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output, "-reffreqtab", reffreqtab]
        codon_optimizer()

        a = open(reference_output_fasta, 'r')
        b = open(output, 'r')
        a_read = a.read()
        b_read = b.read()
        a.close()
        b.close()
        self.assertFalse(a_read == b_read)

    def test_log_creation(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        input_fasta = execution_path + self.CODING_CASES
        output = execution_path + self.OUTPUT_FOR_CODING_CASES
        log = execution_path + "/output/codon_optimizer.log"
        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        self.assertTrue(os.path.exists(log))

    def test_input_file_does_not_exists(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        input_fasta = execution_path + "/file_that_does_not_exists.fasta"
        output = execution_path + self.MIX_CODING_CASES

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]

        with self.assertRaises(IOError):
            codon_optimizer()

    def test_input_rna_sequences(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + self.MIX_CODING_CASES
        input_fasta = execution_path + "/cases/input_coding_rna.fa"
        output = execution_path + self.OUTPUT_FOR_CODING_CASES

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        a = open(reference_output_fasta, 'r')
        b = open(output, 'r')
        a_read = a.read()
        b_read = b.read()
        a.close()
        b.close()
        self.assertTrue(a_read == b_read)


class TestSequenceControl(unittest.TestCase):
    def test_sequence_validation(self):
        self.assertTrue(Sequence("ATGCCCTGGGT").is_valid())
        self.assertFalse(Sequence("ATGCCCTGGGT*)*%").is_valid())
        self.assertFalse(Sequence("ATGCCCTGGGTATGA%").is_valid())
        self.assertTrue(Sequence("ATGCCCTGGGAAATGA").is_valid())

    def test_sequence_characters(self):
        self.assertFalse(Sequence("").is_nucleotide())
        self.assertFalse(Sequence("ATGC+%^&()*)*%").is_nucleotide())
        self.assertTrue(Sequence('ATGCCCTGG').is_nucleotide())

    def test_sequence_length(self):
        self.assertFalse(Sequence("AA").is_longer_than(2), 2)
        self.assertTrue(Sequence("AAAAAA").is_longer_than(2), 2)

    def test_multipleofthree(self):
        self.assertTrue(Sequence('ATG').isx3())
        self.assertFalse(Sequence('A').isx3())

    def test_stop_codon(self):
        self.assertTrue(Sequence('TGAATG').is_stop_codon_in_coding_frame())
        self.assertTrue(Sequence('TgAATG').is_stop_codon_in_coding_frame())
        self.assertFalse(Sequence('ATGA').is_stop_codon_in_coding_frame())
        self.assertFalse(Sequence('aTgA').is_stop_codon_in_coding_frame())
        self.assertFalse(Sequence('A').is_stop_codon_in_coding_frame())
        self.assertTrue(Sequence('ATGTGAtTGA').is_stop_codon_at_its_3prime())

    def test_methionine(self):
        self.assertTrue(Sequence('ATGTGA').is_met_codon_at_its_5prime())


def tests():
    unittest.main()