import unittest
from codon_optimizer import *
import filecmp
import sys
import os

class TestCodonOptimizer(unittest.TestCase):
    def test_input_output_coding_sequence(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + "/cases/coding_seq_codon_optimized_formatted.fasta"
        input_fasta = execution_path + "/cases/input_coding_sequences.fa"
        output = execution_path + "/output/coding_seq_codon_optimized.fasta"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        self.assertTrue(filecmp.cmp(reference_output_fasta, output))

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

        reference_output_fasta = execution_path + "/cases/mix_coding_not_coding_optimized_formatted.fasta"
        input_fasta = execution_path + "/cases/input_mix_coding_not_coding.fasta"
        output = execution_path + "/output/mix_coding_not_coding_optimized.fasta"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        self.assertTrue(filecmp.cmp(reference_output_fasta, output))

    def test_input_output_mix_coding_not_coding_sequence_different_ref_freqtab(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + "/cases/mix_coding_not_coding_optimized_formatted.fasta"
        input_fasta = execution_path + "/cases/input_mix_coding_not_coding.fasta"
        output = execution_path + "/output/mix_coding_not_coding_optimized.fasta"
        reffreqtab = execution_path + "/cases/Mus_musculus_codon_frequency.csv"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output, "-reffreqtab", reffreqtab]
        codon_optimizer()

        self.assertFalse(filecmp.cmp(reference_output_fasta, output))

    def test_log_creation(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + "/cases/mix_coding_not_coding_optimized_formatted.fasta"
        input_fasta = execution_path + "/cases/input_mix_coding_not_coding.fasta"
        output = execution_path + "/output/mix_coding_not_coding_optimized.fasta"
        log = execution_path + "/output/codon_optimizer.log"
        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        self.assertTrue(os.path.exists(log))

    def test_input_file_does_not_exists(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + "/cases/mix_coding_not_coding_optimized_formatted.fasta"
        input_fasta = execution_path + "/file_that_does_not_exists.fasta"
        output = execution_path + "/output/mix_coding_not_coding_optimized.fasta"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]

        with self.assertRaises(IOError):
            codon_optimizer()

    def test_input_rna_sequences(self):
        execution_path = os.path.dirname(__file__)
        os.chdir(execution_path)

        reference_output_fasta = execution_path + "/cases/coding_seq_codon_optimized_formatted.fasta"
        input_fasta = execution_path + "/cases/input_coding_rna.fa"
        output = execution_path + "/output/mix_coding_not_coding_optimized.fasta"

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]

        sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output]
        codon_optimizer()

        self.assertTrue(filecmp.cmp(reference_output_fasta, output))


class TestSequenceControl(unittest.TestCase):

    # def test_input_format(self):
    #     self.assertTrue(File().isFastaFormat())

    def test_sequence_validation(self):
        self.assertTrue(Sequence("ATGCCCTGGGT").isValid())
        self.assertFalse(Sequence("ATGCCCTGGGT*)*%").isValid())
        self.assertFalse(Sequence("ATGCCCTGGGTATGA%").isValid())
        self.assertTrue(Sequence("ATGCCCTGGGAAATGA").isValid())

    def test_sequence_characters(self):
        self.assertFalse(Sequence("").isNucleotide())
        self.assertFalse(Sequence("ATGC+%^&()*)*%").isNucleotide())
        self.assertTrue(Sequence('ATGCCCTGG').isNucleotide())

    def test_sequence_length(self):
        self.assertFalse(Sequence("AA").isLongerThan(2), 2)
        self.assertTrue(Sequence("AAAAAA").isLongerThan(2), 2)

    def test_multipleofthree(self):
        self.assertTrue(Sequence('ATG').isx3())
        self.assertFalse(Sequence('A').isx3())

    def test_stop_codon(self):
        self.assertTrue(Sequence('TGAATG').isStopCodonInCodingFrame())
        self.assertTrue(Sequence('TgAATG').isStopCodonInCodingFrame())
        self.assertFalse(Sequence('ATGA').isStopCodonInCodingFrame())
        self.assertFalse(Sequence('aTgA').isStopCodonInCodingFrame())
        self.assertFalse(Sequence('A').isStopCodonInCodingFrame())
        self.assertTrue(Sequence('ATGTGAtTGA').isStopCodonAtIts3prime())

    def test_methionine(self):
        self.assertTrue(Sequence('ATGTGA').isMetCodonAtIts5prime())

if __name__ == '__main__':
    unittest.main()