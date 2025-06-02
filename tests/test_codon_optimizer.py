import pytest
import os
from codon_optimizer.codon_optimizer import codon_optimizer, Sequence
import filecmp
import sys


@pytest.fixture
def setup_test_data():
    execution_path = os.path.dirname(__file__)
    input_fasta = os.path.join(execution_path, "cases", "input_coding_sequences.fa")
    output_fasta = os.path.join(execution_path, "output", "coding_seq_codon_optimized.fasta")
    reference_output_fasta = os.path.join(execution_path, "cases", "coding_seq_codon_optimized_formatted.fasta")
    
    # Clean the output file before starting the test
    if os.path.exists(output_fasta):
        os.remove(output_fasta)

    return input_fasta, output_fasta, reference_output_fasta


def test_input_output_coding_sequence(setup_test_data):
    input_fasta, output_fasta, reference_output_fasta = setup_test_data
    
    # Execute codon_optimizer
    sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output_fasta]
    codon_optimizer()

    # Compare output with the expected reference output
    assert filecmp.cmp(output_fasta, reference_output_fasta), "Output does not match expected result"


def test_input_output_not_coding_sequence():
    execution_path = os.path.dirname(__file__)
    input_fasta = os.path.join(execution_path, "cases", "input_not_coding_sequences.fa")
    output_fasta = os.path.join(execution_path, "output", "not_coding_seq_codon_optimized.fasta")

    # Execute codon_optimizer
    sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output_fasta]
    codon_optimizer()

    assert not os.path.exists(output_fasta), "Output file should not exist for non-coding sequences"


def test_input_output_mix_coding_not_coding_sequence(setup_test_data):
    input_fasta, output_fasta, reference_output_fasta = setup_test_data

    # Execute codon_optimizer
    sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output_fasta]
    codon_optimizer()

    # Compare output with the expected reference output
    assert filecmp.cmp(output_fasta, reference_output_fasta), "Output does not match expected result"


def test_log_creation():
    execution_path = os.path.dirname(__file__)
    log_file = os.path.join(execution_path, "output", "codon_optimizer.log")
    input_fasta = os.path.join(execution_path, "cases", "input_coding_sequences.fa")
    output_fasta = os.path.join(execution_path, "output", "coding_seq_codon_optimized.fasta")
    
    # Execute codon_optimizer
    sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output_fasta]
    codon_optimizer()

    assert os.path.exists(log_file), "Log file not created"


def test_input_file_does_not_exist():
    execution_path = os.path.dirname(__file__)
    input_fasta = os.path.join(execution_path, "file_that_does_not_exist.fasta")
    output_fasta = os.path.join(execution_path, "output", "coding_seq_codon_optimized.fasta")

    sys.argv = [sys.argv[0], "-input", input_fasta, "-output", output_fasta]
    
    with pytest.raises(IOError):
        codon_optimizer()


def test_input_rna_sequences(setup_test_data):
    input_fasta, output_fasta, reference_output_fasta = setup_test_data
    
    # RNA input case
    input_fasta_rna = os.path.join(os.path.dirname(__file__), "cases", "input_coding_rna.fa")

    # Execute codon_optimizer
    sys.argv = [sys.argv[0], "-input", input_fasta_rna, "-output", output_fasta]
    codon_optimizer()

    # Compare output with the expected reference output
    assert filecmp.cmp(output_fasta, reference_output_fasta), "Output does not match expected result"


class TestSequenceControl:
    def test_sequence_validation(self):
        assert Sequence("ATGCCCTGGGT").is_valid()
        assert not Sequence("ATGCCCTGGGT*)*%").is_valid()
        assert not Sequence("ATGCCCTGGGTATGA%").is_valid()
        assert Sequence("ATGCCCTGGGAAATGA").is_valid()

    def test_sequence_characters(self):
        assert not Sequence("").is_nucleotide()
        assert not Sequence("ATGC+%^&()*)*%").is_nucleotide()
        assert Sequence('ATGCCCTGG').is_nucleotide()

    def test_sequence_length(self):
        assert not Sequence("AA").is_longer_than(2)
        assert Sequence("AAAAAA").is_longer_than(2)

    def test_multipleofthree(self):
        assert Sequence('ATG').isx3()
        assert not Sequence('A').isx3()

    def test_stop_codon(self):
        assert Sequence('TGAATG').is_stop_codon_in_coding_frame()
        assert Sequence('TgAATG').is_stop_codon_in_coding_frame()
        assert not Sequence('ATGA').is_stop_codon_in_coding_frame()
        assert not Sequence('aTgA').is_stop_codon_in_coding_frame()
        assert not Sequence('A').is_stop_codon_in_coding_frame()
        assert Sequence('ATGTGAtTGA').is_stop_codon_at_its_3prime()

    def test_methionine(self):
        assert Sequence('ATGTGA').is_met_codon_at_its_5prime()
