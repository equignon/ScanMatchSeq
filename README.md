# ScanMatchSeq
ScanMatchSeq is a simple python script to look for regions reaching a certain value between two RNAs.

It runs through the first RNA with a sliding window and test the values of the positions inside it. If enough values meet the criteria,
a complementary sequence is searched on the other RNA (*Note: for the complementarity, G-U base pairing is taken in account*). If a complementary sequence
is found on the other RNA and if it also matches the criteria.
This script produces multiple output.
- "output.csv" contains all the sequences that matches at least one sequence of the other RNA with the corresponding positions and sequence on both RNAs.
- "chord_diag_A.csv" and "chord_diag_B.csv" are meant to be used with R to draw circular representation of the potential interactions.


## Requirements

- numpy
- pandas
- Biopython

*Note for Mac ARM users: installing Biopython might install the x64 version, if it's the case, it's necessary to install Biopython from source.*

## How to use

You can run the script by typing `python scanmatchseq.py` in your console/terminal.

The following inputs are required:
- FASTA file of your first RNA sequence (.fasta)
- FASTA file of your second RNA sequence (.fasta)
- Values of your first RNA (.txt with positions and values separated by a tabulation)
- Values of your second RNA (.txt with positions and values separated by a tabulation)
- Value searched in the sliding window

Some parameters can be changed:
- `--window` allows to modify the size of the sliding window
- `--threshold` allows to change the threshold at which the window considered matches criteria
- `--test_type` allows to choose the type of test conducted with the value searched. Possible tests are: 'equal', 'superior', 'inferior', 'superior or equal' and 'inferior or equal'
