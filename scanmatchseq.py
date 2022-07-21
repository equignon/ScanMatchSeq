#!/usr/bin/python

import sys
import os
import numpy as np
import pandas as pd
import csv
import pathlib
import argparse

from Bio import SeqIO # For arm Mac user: might need to install Biopython from source as pip seems to install x64 version

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("fasta_A", help="FASTA file of first sequence")
    parser.add_argument("fasta_B", help="FASTA file of second sequence")
    parser.add_argument("values_A", help="Values corresponding to first sequence")
    parser.add_argument("values_B", help="Values corresponding to second sequence")
    parser.add_argument("--window", help="Window size", type=int, default=4)
    parser.add_argument("--value_searched", help="Value searched in the two datasets", type=float, default=1.0)
    parser.add_argument("--threshold", help="Threshold for the sliding window", type=int, default=2)
    parser.add_argument("--test_type", help="Type of test to realize ('equal', 'superior', 'inferior', 'superior or equal', 'inferior or equal')", type=str, default="equal")

    args = parser.parse_args()

    # Filling dataframes with files given

    df_seq_A = pd.DataFrame()
    df_seq_A = add_dataset(df_seq_A, args.fasta_A, args.values_A)

    df_seq_B = pd.DataFrame()
    df_seq_B = add_dataset(df_seq_B, args.fasta_B, args.values_B)

    # Sliding window on first sequence to look for sequences containing at least x times the searched value (default = 2)

    df_results = pd.DataFrame(columns=['Position A', 'Sequence A', 'Sequence B', 'Position B'])

    for i, row in df_seq_A.iterrows():
        sequence_A = ""
        j = i

        if i == ((len(df_seq_A) - args.window) + 2):
            break

        parse_A = check_value(args, df_seq_A, sequence_A, i, j)
        sequence_A = parse_A[0]
        nb_found = parse_A[1]

        if(nb_found >= args.threshold):
            for k, row_b in df_seq_B.iterrows():
                sequence_B = ""
                j = k

                if k == ((len(df_seq_B) - args.window) + 2):
                    break

                parse_B = check_value(args, df_seq_B, sequence_B, k, j)
                sequence_B = parse_B[0]
                nb_found_B = parse_B[1]

                if nb_found_B >= args.threshold:
                    seq_equal = compare_sequence(args, sequence_A, sequence_B)

                    if seq_equal == True:
                        pos_A = str(i) + " - " + str(i + args.window - 1)
                        pos_B = str(k) + " - " + str(k + args.window - 1)
                        new_row = pd.DataFrame([{'Position A': pos_A, 'Sequence A': sequence_A, 'Sequence B': sequence_B, 'Position B': pos_B}])
                        df_results = pd.concat([df_results, new_row], axis = 0)
                    
                    sequence_B = sequence_B[::-1]
                    seq_equal = compare_sequence(args, sequence_A, sequence_B)

                    if seq_equal == True:
                        pos_A = str(i) + " - " + str(i + args.window - 1)
                        pos_B = str(k) + " - " + str(k + args.window - 1)
                        new_row = pd.DataFrame([{'Position A': pos_A, 'Sequence A': sequence_A, 'Sequence B': sequence_B, 'Position B': pos_B}])
                        df_results = pd.concat([df_results, new_row], axis = 0)

    df_results.to_csv(r'output.csv', sep='\t', index = False, header = True)

    df_chord_A = pd.DataFrame()
    df_chord_B = pd.DataFrame()
    for i, row in df_results.iterrows():
        pos_A = row[0].split(" - ")
        start_A = int(pos_A[0])
        end_A = int(pos_A[1])
        new_row = pd.DataFrame([{'Sequence':"SeqA", 'Sequence_A': start_A, 'Sequence_B': end_A}])
        df_chord_A = pd.concat([df_chord_A, new_row], axis = 0)

        pos_B = row[3].split(" - ")
        start_B = int(pos_B[0])
        end_B = int(pos_B[1])
        new_row = pd.DataFrame([{'Sequence':"SeqB", 'Sequence_A': start_B, 'Sequence_B': end_B}])
        df_chord_B = pd.concat([df_chord_B, new_row], axis = 0)

    df_chord_A.to_csv(r'chord_diag_A.csv', sep='\t', index = False, header = True)
    df_chord_B.to_csv(r'chord_diag_B.csv', sep='\t', index = False, header = True)

def add_dataset(df, file, value):
    for seq_record in SeqIO.parse(file, "fasta"):
        df = pd.DataFrame(index = np.arange(1, len(seq_record.seq)+1), columns = ['Nucleotide', 'Reactivity'])
        x = 1
        for p in seq_record.seq:
            df.at[x, 'Nucleotide'] = p
            x += 1
        with open(value) as f:
            for line in csv.reader(f, dialect = "excel-tab"):
                if(line[1] == "" or line[1] == "#N/A"):
                    df.at[int(line[0]), 'Reactivity'] = float('NaN')
                else:
                    df.at[int(line[0]), 'Reactivity'] = float(line[1])
    return df

def check_value(args, df, seq, i, j):
    nb_found = 0
    if args.test_type == "equal":
        while j < (i + args.window):
            # Generating sequence contained in the windows for sequence A
            seq += df.at[j, 'Nucleotide']

            # Counting the number of nucleotides that reaches the threshold given
            if(pd.isna(df.at[j, 'Reactivity'])):
                pass
            else:
                if(float(df.at[j, 'Reactivity']) == args.value_searched):
                    nb_found += 1
            j += 1

    if args.test_type == "superior":
        while j < (i + args.window):
            # Generating sequence contained in the windows for sequence A
            seq += df.at[j, 'Nucleotide']

            # Counting the number of nucleotides that reaches the threshold given
            if(pd.isna(df.at[j, 'Reactivity'])):
                pass
            else:
                if(float(df.at[j, 'Reactivity']) > args.value_searched):
                    nb_found += 1
            j += 1

    if args.test_type == "inferior":
        while j < (i + args.window):
            # Generating sequence contained in the windows for sequence A
            seq += df.at[j, 'Nucleotide']

            # Counting the number of nucleotides that reaches the threshold given
            if(pd.isna(df.at[j, 'Reactivity'])):
                pass
            else:
                if(float(df.at[j, 'Reactivity']) < args.value_searched):
                    nb_found += 1
            j += 1

    if args.test_type == "superior or equal":
        while j < (i + args.window):
            # Generating sequence contained in the windows for sequence A
            seq += df.at[j, 'Nucleotide']

            # Counting the number of nucleotides that reaches the threshold given
            if(pd.isna(df.at[j, 'Reactivity'])):
                pass
            else:
                if(float(df.at[j, 'Reactivity']) >= args.value_searched):
                    nb_found += 1
            j += 1


    if args.test_type == "inferior or equal":
        while j < (i + args.window):
            # Generating sequence contained in the windows for sequence A
            seq += df.at[j, 'Nucleotide']

            # Counting the number of nucleotides that reaches the threshold given
            if(pd.isna(df.at[j, 'Reactivity'])):
                pass
            else:
                if(float(df.at[j, 'Reactivity']) <= args.value_searched):
                    nb_found += 1
            j += 1

    return seq, nb_found

def compare_sequence(args, seq_A, seq_B):
    c = 1
    length = 0
    while c <= args.window:
        if(seq_A[c-1] == "A"):
            if(seq_B)[c-1] == "U":
                length += 1
        else:
            pass
        if(seq_A[c-1] == "U"):
            if seq_B[c-1] == "A" or seq_B[c-1] == "G":
                length += 1
        else:
            pass
        if(seq_A[c-1] == "G"):
            if seq_B[c-1] == "C" or seq_B[c-1] == "U":
                length += 1
            else:
                pass
        if(seq_A[c-1] == "C"):
            if seq_B[c-1] == "G":
                length += 1
            else:
                pass
        c += 1
    if(length == args.window):
        return True
    else:
        return False

if __name__ == '__main__':
    main(sys.argv)