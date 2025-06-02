# codon_optimizer
A simple tool to produce DNA sequences with the most abundant codons in a specific organism.

## Parameters
*  `-h, --help`            Show this help message and exit.
*  `-input FILENAME`       Input path to the fasta file containing sequences to optimize (default: local).
*  `-reffreqtab FILENAME`  Input path to the reference frequency table. The table should be in a specific format (see the `reference_seq` folder) (default: Homo sapiens).
*  `-output FILENAME`      Output fasta file with codon optimized sequences (default: local).
*  `-nolog`                Disable log output (default: False).
*  `--version`             Show program's version number and exit.

## Installation

1. **Clone the repository**:
    ```bash
    git clone https://github.com/jeromelane/codon_optimizer.git
    cd codon_optimizer
    ```

2. **Create a virtual environment and install dependencies**:
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    pip install -e .
    ```

## Usage

1. **Most simple execution**:
    ```bash
    codon_optimizer -input myfile.fa 
    ```
    For example, `myfile.fa` could be `/your/absolute/path/to/tests/cases/input_coding_rna.fa`.


2. **Disable logging**:
    ```bash
    codon_optimizer -input myfile.fa -nolog
    ```

3. **Use a specific codon usage table**:
    ```bash
    codon_optimizer -input myfile.fa -reffreqtab codonusage.txt
    ```

## Testing

1. **Install `pytest`**:
    ```bash
    source venv/bin/activate
    pip install pytest
    ```

2. **Run tests** using `pytest`:
    ```bash
    pytest tests
    ```

    This will run the test suite and display the results in the terminal.

3. **Log file creation**: During execution, a log file will be created at `output/codon_optimizer.log` if logging is enabled.
