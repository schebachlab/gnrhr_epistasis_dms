#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script used to sort, filter, and get average across replicates for Laura's mGNRHR DMS scan
@author: Charles Kuntz :: cpkuntz@iu.edu
"""

import sys
import os
import re
import textwrap
import string
from itertools import groupby
import numpy as np
import pandas as pd
import scipy.stats as sp


def import_file(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    return data_raw

def import_list(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.rstrip())
    return data_raw


def import_text_file(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    seqs_str_list = []
    for line in data_raw:
        item = line[0]
        seqs_str_list.append(item)
    return seqs_str_list


def print_list_to_text(output_array, outname):
    output_txt = open(outname+'.txt', 'w')
    for line in output_array:
        outline = str(line)
        print(outline, file=output_txt)
    output_txt.close()


def extract_fasta_data(data_raw):
    protein_sequence = []
    idx = 1
    while idx < len(data_raw):
        line_seq_str = data_raw[idx][0]
        line_seq_list = list(line_seq_str)
        for residue in line_seq_list:
            protein_sequence.append(residue)
        idx += 1
    return protein_sequence


def import_csv(path_to_file, wt_seq):
    file_obj = open(path_to_file)
    csv_dict = {}
    data_raw = []
    for line in file_obj:
        line = line.rstrip()
        line = line.split(',')
        variant_label = line[0].strip()
        csv_dict[variant_label] = line[1:]
        data_raw.append(line[1:])
    len_data = len(data_raw[1])
    #csv_ndarray = np.genfromtxt(path_to_file, delimiter=',')
    csv_df = pd.read_csv(path_to_file)
    return csv_dict, len_data, csv_df

def get_ints_ranks(csv_df):
    # Base percentile on non-normalized intensity score
    csv_df['rank'] = csv_df['score_ints'].rank(ascending=True)
    num_rows = csv_df.shape[0]
    csv_df['percentile'] = (csv_df['rank'] - 1) / num_rows * 100
    ranks_dict = csv_df.set_index('var').T.to_dict('list')
    return ranks_dict

def csv_average(csv_dict1, csv_dict2, csv_df1, csv_df2, threshold, outname):
    logfile = open(outname+'_spearman.txt', 'w')
    ints_ndarray1 = csv_df1[['score_ints']].values
    ints_ndarray2 = csv_df2[['score_ints']].values
    before_spearman_r, before_spearman_pvalue = sp.spearmanr(ints_ndarray1, ints_ndarray2)
    print('Spearman\'s r before data filtering: ', before_spearman_r, file=logfile)
    print('Spearman\'s p-value before data filtering: ', before_spearman_pvalue, file=logfile)
    print('Number of variants before data filtering: ', csv_df1.shape[0], file=logfile)
    variant_labels1 = [key for key in csv_dict1]
    variant_labels2 = [key for key in csv_dict2]
    set_variant_labels = set(variant_labels1) | set(variant_labels2)
    set_variant_labels = set_variant_labels - set(['var', 'variant'])
    list_variant_labels = list(set_variant_labels)
    csv_dict_avg = {}
    ranks_dict1 = get_ints_ranks(csv_df1)
    ranks_dict2 = get_ints_ranks(csv_df2)
    passing_score_ints1 = []
    passing_score_ints2 = []
    for variant_label in list_variant_labels:
        csv_data1 = csv_dict1[variant_label]
        csv_data2 = csv_dict2[variant_label]
        position = int(variant_label[1:-1])
        score_norm1 = float(csv_data1[1])
        score_norm2 = float(csv_data2[1])
        score_ints1 = float(csv_data1[2])
        score_ints2 = float(csv_data2[2])
        cts_avg1 = float(csv_data1[3])
        cts_avg2 = float(csv_data2[3])
        ints_rank1 = ranks_dict1[variant_label][-2]
        ints_rank2 = ranks_dict2[variant_label][-2]
        perc_rank1 = ranks_dict1[variant_label][-1]
        perc_rank2 = ranks_dict2[variant_label][-1]
        d = ints_rank1 - ints_rank2
        d_squared = d**2
        diff_perc = abs(perc_rank1 - perc_rank2)
        score_quality_test = 1 if all([diff_perc < threshold, cts_avg1 > 50 and cts_avg2 > 50]) else 0
        if score_quality_test == 1:
            passing_score_ints1.append(score_ints1)
            passing_score_ints2.append(score_ints2)
        score_norm_avg = np.mean([score_norm1, score_norm2])
        score_ints_avg = np.mean([score_ints1, score_ints2])
        rar_cts_avg = np.mean([cts_avg1, cts_avg2])
        csv_data_avg = [position, score_norm1, score_norm2, score_norm_avg, score_ints1, ints_rank1, perc_rank1, score_ints2, ints_rank2, perc_rank2, d_squared, diff_perc, score_ints_avg, cts_avg1, cts_avg2, rar_cts_avg, score_quality_test]
        csv_dict_avg[variant_label] = csv_data_avg
    passing_scores_array1 = np.array(passing_score_ints1)
    passing_scores_array2 = np.array(passing_score_ints2)
    after_spearman_r, after_spearman_pvalue = sp.spearmanr(passing_scores_array1, passing_scores_array2)
    print('Spearman\'s r after data filtering: ', after_spearman_r, file=logfile)
    print('Spearman\'s p-value after data filtering: ', after_spearman_pvalue, file=logfile)
    print('Number of variants after data filtering: ', len(passing_score_ints1), file=logfile)
    logfile.close()
    len_csv_avg = len(csv_dict_avg[list_variant_labels[0]])
    return csv_dict_avg, len_csv_avg

def expand_sequence(sequence, alphabet_order):
    mut_list = list(alphabet_order)
    expanded_position_list = []
    labels = []
    expanded_wt_seq = []
    expanded_mut_seq = []
    for idx in range(len(sequence)):
        position = idx + 1
        position_label = str(position)
        wt_label = sequence[idx]
        for mut_idx in range(len(mut_list)):
            mut_label = mut_list[mut_idx]
            label_list = [wt_label, position_label, mut_label]
            variant_label = ''.join(label_list)
            expanded_position_list.append(position_label)
            labels.append(variant_label)
            expanded_wt_seq.append(wt_label)
            expanded_mut_seq.append(mut_label)
    return expanded_position_list, labels, expanded_wt_seq, expanded_mut_seq


def get_csv_data(tm_positions, expanded_position_list, expanded_wt_seq, expanded_mut_seq, label_list, csv_dict, len_data):
    empty_line = ['nan' for x in range(len_data)]
    #empty_line = ['nan', 'nan']
    csv_data = [['Variant', 'position', 'wt_res', 'mut_res', 'idx', 'tm_res?', 'position', 'score_norm1', 'score_norm2', 'score_norm_avg', 'score_ints1', 'score_ints1_rank', 'score_ints1_percentile', 'score_ints2', 'score_ints2_rank', 'score_ints2_percentile', 'diff_rank_squared', 'diff_percentile', 'score_ints_avg', 'cts_avg1', 'cts_avg2', 'rar_cts_avg', 'pass_qual?']]
    for idx in range(len(expanded_position_list)):
        label = label_list[idx]
        position = label[1:-1]
        if position in tm_positions:
            tm_status = 'TMD'
        else:
            tm_status = 'SOL'
        data_list = [label_list[idx], expanded_position_list[idx], expanded_wt_seq[idx], expanded_mut_seq[idx], idx+1, tm_status]
        if label in csv_dict:
            data_output = csv_dict[label]
        else:
            data_output = empty_line
        for entry in data_output:
            data_list.append(entry)
        csv_data.append(data_list)
    return csv_data


def print_results(output_data, outname):
    myfile = open(outname+'.csv', 'w')
    for line in output_data:
        line_string = [str(x) for x in line]
        csv_line = ','.join(line_string)
        print(csv_line, file = myfile)
    myfile.close()


def interactive_configuration():
    path_fasta = str(input("Please enter a path to the FASTA file: "))
    path_csv1 = str(input("Please enter a path to the CSV file containing the data for replicate 1: "))
    path_csv2 = str(input("Please enter a path to the CSV file containing the data for replicate 2: "))
    path_tm = str(input("Please enter a path to the TMD positions list: "))
    threshold = float(input("Please enter a percentile difference threshold that the replicate values must fall within in order to pass quality filters: "))
    print("How do you want to order the entries of the one-hot vector?")
    print("Enter a string of 20 amino acids, using single-letter labels.")
    print("Default string for us is generally von Heijne's order of polarity:\n")
    print("    ILFVCMAWTYGSNHPQREKD*\n")
    print("Another choice is the default order used in most PSSM files, alphabetical order:\n")
    print("    ARNDCQEGHILKMFPSTWYV*\n")
    new_order_str = str(input("Enter the sequence: "))
    jobname = str(input("Please enter a name for this job: "))
    return path_fasta, path_csv1, path_csv2, path_tm, threshold, new_order_str, jobname


if __name__ == '__main__':


    path_fasta, path_csv1, path_csv2, path_tm, threshold, new_order_str, jobname = interactive_configuration()

    tm_positions = import_list(path_tm)

    fasta_raw = import_file(path_fasta)

    protein_seq = extract_fasta_data(fasta_raw)

    expanded_position_list, labels, expanded_wt_seq, expanded_mut_seq = expand_sequence(protein_seq, new_order_str)

    csv_dict1, len_data1, csv_df1 = import_csv(path_csv1, protein_seq)
    csv_dict2, len_data2, csv_df2 = import_csv(path_csv2, protein_seq)

    csv_dict_avg, len_csv_avg = csv_average(csv_dict1, csv_dict2, csv_df1, csv_df2, threshold, jobname)


    csv_data_avg = get_csv_data(tm_positions, expanded_position_list, expanded_wt_seq, expanded_mut_seq, labels, csv_dict_avg, len_csv_avg) 


    print_results(csv_data_avg, jobname)

