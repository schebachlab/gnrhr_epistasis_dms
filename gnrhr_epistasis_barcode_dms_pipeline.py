#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import re
import textwrap
import string
import random
import numpy as np
import pandas as pd
import scipy.stats as stat
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
import multiprocessing as mp
import difflib
from datetime import datetime
from collections import Counter
from collections import defaultdict
from configparser import ConfigParser


def qualityfilter3_original(data, threshold=30, errors=1):
    filtered_data = []
    avg_scores = []
    passed = 0
    failed = 0
    i = 2
    while i < len(data):
        domain_seq = data[i][0:]
        domain_qual = data[i+2][0:]
        domain_qual_num = [ord(q_symbol) - 33 for q_symbol in domain_qual]
        prob_error_read = [10**-(base_q/10) for base_q in domain_qual_num]
        sum_prob_error_read = np.sum(prob_error_read)
        avg_score_domain = np.mean(domain_qual_num)
        pass_conditions = [sum_prob_error_read < errors,
                        avg_score_domain >= threshold,
                        'N' not in domain_seq]
        if all(pass_conditions):
            filtered_data.append(domain_seq)
            avg_scores.append(avg_score_domain)
            passed += 1
            i += 5
        else:
            failed += 1
            i += 5
    percent_passed = (passed/(passed + failed))*100
    mean_scores = np.mean(avg_scores)
    total_reads = passed+failed
    return filtered_data, mean_scores, percent_passed, total_reads


def qualityfilter3_dict(data, threshold, errors=1):
    filtered_data = {}
    avg_scores = []
    passing_coords = []
    passed = 0
    failed = 0
    i = 2
    while i < len(data):
        domain_coord = data[i-2]
        domain_seq = data[i][0:]
        domain_qual = data[i+2][0:]
        domain_qual_num = [ord(q_symbol) - 33 for q_symbol in domain_qual]
        prob_error_read = [10**-(base_q/10) for base_q in domain_qual_num]
        sum_prob_error_read = np.sum(prob_error_read)
        avg_score_domain = np.mean(domain_qual_num)
        pass_conditions = [sum_prob_error_read <= errors,
                           avg_score_domain >= threshold,
                           'N' not in domain_seq]
        if all(pass_conditions):
            filtered_data.update({domain_coord:domain_seq})
            passing_coords.append(domain_coord)
            avg_scores.append(avg_score_domain)
            passed += 1
            i += 5
        else:
            failed += 1
            i += 5
    percent_passed = (passed/(passed + failed))*100
    mean_scores = np.mean(avg_scores)
    total_reads = passed+failed
    return filtered_data, mean_scores, percent_passed, total_reads, passing_coords

def qualityfilter3b(data_r1, data_r2, threshold, errors=1):
    
    filtered_data_r1, avg_scores_r1, passed_r1, num_reads_r1, passing_coords_r1 = qualityfilter3_dict(data_r1, threshold, errors)
    filtered_data_r2, avg_scores_r2, passed_r2, num_reads_r2, passing_coords_r2 = qualityfilter3_dict(data_r2, threshold, errors)

    set_s1_coords = set(passing_coords_r1)
    set_s2_coords = set(passing_coords_r2)
    intersection_sets = set_s1_coords.intersection(set_s2_coords)
    list_intersection = list(intersection_sets)
    filtered_data_combined = defaultdict(list)
    for mydict in (filtered_data_r1, filtered_data_r2):
        for key, value in mydict.items():
            filtered_data_combined[key].append(value)
    
    filtered_data_dict = dict(filtered_data_combined)
    # return list_intersection for debugging in interactive sessions
    return filtered_data_dict, avg_scores_r1, avg_scores_r2, passed_r1, passed_r2, num_reads_r1, num_reads_r2

def readlength_counter(fastq_data_raw):
    read_lengths = []
    idx = 0
    while idx < len(fastq_data_raw):
        read_len = len(fastq_data_raw[idx+2])
        read_lengths.append(read_len)
        idx += 5
    read_lengths_count = Counter(read_lengths)
    return read_lengths_count

def fastq_parser(path_to_fastq):
    filehandle = open(path_to_fastq, 'r')
    fastq_datastring = filehandle.read()
    delimiter = re.compile(r'[\S]+')
    fastq_datalist = delimiter.findall(fastq_datastring)
    filehandle.close()
    return fastq_datalist


def reverse_complement_func(forward_sequence):
    forward_alphabet = 'AGCTagct'
    revcomp_alphabet = 'TCGAtcga'
    reverse_sequence = forward_sequence[::-1]
    reverse_complement_sequence = reverse_sequence.translate({ord(x):
        y for (x, y) in zip(forward_alphabet, revcomp_alphabet)})
    return reverse_complement_sequence


def write_text_file(data, outname):
    text_file = open(outname+'_out.txt', 'w')
    for entry in data:
        print(entry, file=text_file)
    text_file.close()


def flanking_seq_finder(consensus_seqs, flanking_sequence):
    sequence_start_list = []
    failures = 0
    successes = 0
    for sequence_read in consensus_seqs:
        if flanking_sequence in sequence_read:
            sequence_start = sequence_read.find(flanking_sequence)
            sequence_start_list.append(sequence_start)
            successes += 1
        else:
            failures += 1
    seq_register_count = Counter(sequence_start_list)
    idx_flanking_seq = seq_register_count.most_common(1)[0][0]
    return sequence_start_list, successes, failures, seq_register_count, idx_flanking_seq


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


def import_file(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    return data_raw

def longest_match_debug(s1, s2):
    match = [[0]*(1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1+ len(s2)):
            if s1[x-1] == s2[y-1]:
                match[x][y] = match[x-1][y-1] + 1
                if match[x][y] > longest:
                    longest = match[x][y]
                    x_longest = x
            else:
                match[x][y] = 0
    overlap_seq = s1[x_longest - longest: x_longest]
    
    return overlap_seq

def get_overlap_sequencematcher(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size_overlap = s.find_longest_match(0, len(s1), 0, len(s2))
    if pos_b < pos_a:
        unique_to_s1 = s1[0:pos_a]
        unique_to_s2 = s2[size_overlap:]
        overlap_seq = s1[pos_a:pos_a+size_overlap]
        combined_seq = unique_to_s1+overlap_seq+unique_to_s2
    if pos_a < pos_b:
        unique_to_s1 = s1[size_overlap:]
        unique_to_s2 = s2[0:pos_b]
        overlap_seq = s2[pos_b:pos_b+size_overlap]
        combined_seq = unique_to_s2+overlap_seq+unique_to_s1
    if pos_a == pos_b:
        unique_to_s1 = ''
        unique_to_s2 = ''
        overlap_seq = s1
        combined_seq = s1
    return unique_to_s1, unique_to_s2, overlap_seq, combined_seq, size_overlap


def make_consensus_seq(filtered_data_dict, size_overlap_threshold):
    consensus_seqs = []
    overlap_sizes = []
    both_passed_qscore = 0
    for coord, seq_list in filtered_data_dict.items():
        if len(seq_list) == 2:
            both_passed_qscore += 1
            s1 = seq_list[0]
            s2 = reverse_complement_func(seq_list[1])
            unique_to_s1, unique_to_s2, overlap_seq, combined_seq, size_overlap = get_overlap_sequencematcher(s1, s2)
            overlap_sizes.append(size_overlap)
            if size_overlap >= size_overlap_threshold:
                consensus_seqs.append(combined_seq)
    len_dict = len(filtered_data_dict)
    len_filtered = len(consensus_seqs)
    return consensus_seqs, both_passed_qscore, overlap_sizes, len_dict, len_filtered

def make_consensus_seq_debug(filtered_data_dict, size_overlap_threshold):
    consensus_seqs = []
    s1_list = []
    s2_list = []
    overlap_sizes = []
    both_passed_qscore = 0
    for coord, seq_list in filtered_data_dict.items():
        if len(seq_list) == 2:
            both_passed_qscore += 1
            s1 = seq_list[0]
            s2 = reverse_complement_func(seq_list[1])
            unique_to_s1, unique_to_s2, overlap_seq, combined_seq, size_overlap = get_overlap_sequencematcher(s1, s2)
            overlap_sizes.append(size_overlap)
            if size_overlap >= size_overlap_threshold:
                consensus_seqs.append(combined_seq)
                s1_list.append(s1)
                s2_list.append(s2)
    len_dict = len(filtered_data_dict)
    len_filtered = len(consensus_seqs)
    return consensus_seqs, s1_list, s2_list


def isolate_barcodes(consensus_seqs, flanking_5p, flanking_3p, barcode_length):
    barcode_list = []
    reads_count = len(consensus_seqs)
    for read in consensus_seqs:
        if all([flanking_5p in read, flanking_3p in read]):
            re_string = flanking_5p+'(.*)'+flanking_3p
            re_test = re.search(re_string, read)
            if re_test:
                barcode_candidate = re_test.group(1)
                if len(barcode_candidate) in barcode_length:
                    barcode_list.append(barcode_candidate)
    barcode_count = Counter(barcode_list)
    unique_barcodes_count = len(barcode_count)
    csv_data = []
    for barcode, count in barcode_count.items():
        line = [barcode, count]
        csv_data.append(line)
    return barcode_list, unique_barcodes_count, csv_data


def isolate_barcodes_debug(consensus_seqs, flanking_5p, flanking_3p, barcode_length):
    barcode_list = []
    barcode_lengths = []
    reads_count = len(consensus_seqs)
    for read in consensus_seqs:
        if all([flanking_5p in read, flanking_3p in read]):
            idx_5p = read.find(flanking_5p)
            idx_3p = read.find(flanking_3p)
            barcode_candidate = read[idx_5p+len(flanking_5p):idx_3p]
            barcode_lengths.append(len(barcode_candidate))
            if len(barcode_candidate) == barcode_length:
                barcode_list.append(barcode_candidate)
    barcode_count = Counter(barcode_list)
    barcode_lengths_histo = Counter(barcode_lengths)
    unique_barcodes_count = len(barcode_count)
    csv_data = []
    for barcode, count in barcode_count.items():
        line = [barcode, count]
        csv_data.append(line)
    return barcode_list, unique_barcodes_count, csv_data


def print_results(output_data, outname):
    myfile = open(outname+'.csv', 'w')
    for line in output_data:
        line_string = [str(x) for x in line]
        csv_line = ','.join(line_string)
        print(csv_line, file = myfile)
    myfile.close()

def new_illumina_filter(all_barcodes_bin, set_dict_barcodes):
    filtered_barcodes_bin = []
    for barcode in all_barcodes_bin:
        if barcode in set_dict_barcodes:
            filtered_barcodes_bin.append(barcode)
    return filtered_barcodes_bin

def analyze_illumina(barcode_variant_dict, all_illumina_barcodes, variant_list):
    #variant_list = [variant_raw[x][0] for x in range(len(variant_raw))]
    reads_per_mutant_dict = {variant:0 for variant in variant_list}
    reads_per_mutant_dict["None"] = 0
    master_reads_pool = []
    for barcode in all_illumina_barcodes:
        reads_per_mutant_dict[barcode_variant_dict[barcode]] += 1
        if barcode_variant_dict[barcode] != "None":
            master_reads_pool.append(barcode)
    return reads_per_mutant_dict, master_reads_pool

def print_variant_counts(reads_per_mutant_dict, outname):
    variant_counts = []
    for variant, count in reads_per_mutant_dict.items():
        variant_counts.append([variant, count])
    myfile = open(outname+'.csv', 'w')
    for line in variant_counts:
        line_string = [str(x) for x in line]
        csv_line = ','.join(line_string)
        print(csv_line, file = myfile)
    myfile.close()

def convert_reads_barcodes_to_labels(barcode_variant_dict, all_illumina_barcodes):
    illumina_reads_labels_bin = []
    wtcount_bin = 0
    varcount_bin = 0
    none_count_bin = 0
    for barcode in all_illumina_barcodes:
        variant = barcode_variant_dict[barcode]
        if variant == "WT":
            wtcount_bin += 1
            illumina_reads_labels_bin.append(variant)
        elif variant == "None":
            none_count_bin += 1
        else:
            varcount_bin += 1
            illumina_reads_labels_bin.append(variant)
    
    return illumina_reads_labels_bin, wtcount_bin, varcount_bin, none_count_bin



def weighted_average_func(var_df_polarity_bins, intensity_bins, wt_mfi_exp):
    sum_over_bins__varcounts = sum(var_df_polarity_bins)
    weighted_varcounts = [intensity * var_df for intensity, var_df in zip(
            intensity_bins, var_df_polarity_bins)]
    sum_over_bins__weighted_varcounts = sum(weighted_varcounts)
    varcount_intensity_weighted = sum_over_bins__weighted_varcounts / sum_over_bins__varcounts
    var_intensity_over_exp_wt_mfi = varcount_intensity_weighted / wt_mfi_exp
    return varcount_intensity_weighted, var_intensity_over_exp_wt_mfi, sum_over_bins__varcounts


def create_pd_series(rarefied_illumina_labels, variant_list):
    #variant_list = [variant_raw[x][0] for x in range(len(variant_raw))]
    #variant_list_no_wt = variant_list[:-1]
    reads_per_mutant_dict = {variant:0 for variant in variant_list}
    #wtcount = 0
    for label in rarefied_illumina_labels:
        reads_per_mutant_dict[label] += 1
            
    output_series = pd.Series(reads_per_mutant_dict, name='Score')

    return output_series

def parse_config_file(config_file):
    parser = ConfigParser()
    parser.read(config_file)
    lib_diag_choice = parser.getboolean('Basic analysis parameters', 'Run library?')
    reverse_read_choice = parser.getboolean('Basic analysis parameters', 'Use reverse reads?')
    if reverse_read_choice not in [True, False]:
        print("Error. Must choose True to use forward and reverse reads and False to use only forward reads. Check configuration file.")
        sys.exit()
    num_bins = parser.getint('Data from DMS experiment', 'Number of bins')
    if (reverse_read_choice, lib_diag_choice) == (True, True):
        path_lib_r1 = parser.get('Library and NGS parameters', 'R1 library file')
        path_lib_r2 = parser.get('Library and NGS parameters', 'R2 library file')
    if (reverse_read_choice, lib_diag_choice) == (False, True):
        path_lib_r1 = parser.get('Library and NGS parameters', 'R1 library file')
        path_lib_r2 = 'NA'
    if lib_diag_choice == False:
        path_lib_r1 = 'NA'
        path_lib_r2 = 'NA'
    if num_bins < 2:
        print("Error. You must enter at least two bins. Check configuration file.")
        sys.exit()
    num_alpha_dict = dict(zip(range(0, 26), string.ascii_uppercase))
    path_fastq_raw_r1 = parser.get('Data from DMS experiment', 'R1 FastQ files')
    path_fastq_r1 = path_fastq_raw_r1.split()
    if reverse_read_choice == True:
        path_fastq_raw_r2 = parser.get('Data from DMS experiment', 'R2 FastQ files')
        path_fastq_r2 = path_fastq_raw_r2.split()
    else:
        path_fastq_r2 = ['NA' for bin in range(num_bins)]
    if len(path_fastq_r1) != num_bins:
        print("Error. Must choose R1 and R2 FastQ files for each bin. Check configuration file.")
        sys.exit()
    if len(path_fastq_r2) != num_bins:
        print("Error. Must choose R1 and R2 FastQ files for each bin. Check configuration file.")
        sys.exit()
    intensity_raw = parser.get('Data from DMS experiment', 'Bin intensities')
    intensity_string = intensity_raw.split()
    intensity = [float(x) for x in intensity_string]
    if len(intensity) != num_bins:
        print("Error. Must enter a bin intensity for each bin. Check configuration file.")
        sys.exit()
    MFI_list = intensity
    wt_mfi_exp = parser.getfloat('Data from DMS experiment', 'WT MFI')
    barcode_size_raw = parser.get('Library and NGS parameters', 'Barcode length')
    barcode_size_string = barcode_size_raw.split()
    barcode_size = [int(x) for x in barcode_size_string]
    flanking_5p = parser.get('Library and NGS parameters', 'Upstream flanking sequence')
    flanking_3p = parser.get('Library and NGS parameters', 'Downstream flanking sequence')
    path_dict = parser.get('Library and NGS parameters', 'Path to dictionary')
    jobname = parser.get('Basic analysis parameters', 'Job name')
    seed_switch = parser.getboolean('Basic analysis parameters', 'Set seed')
    if seed_switch == True:
        seed = parser.getint('Basic analysis parameters', 'Seed')
    else:
        seed = random.randrange(sys.maxsize)
    qscore_threshold = parser.getfloat('Quality filter', 'Quality score threshold')
    errors = parser.getfloat('Quality filter', 'Maximum expected errors in read')
    return (num_bins,
            num_alpha_dict,
            reverse_read_choice,
            lib_diag_choice,
            path_lib_r1,
            path_lib_r2,
            path_fastq_r1,
            path_fastq_r2,
            MFI_list,
            wt_mfi_exp,
            barcode_size,
            flanking_5p,
            flanking_3p,
            jobname,
            path_dict,
            seed,
            qscore_threshold,
            errors)





if __name__ == '__main__':


    if len(sys.argv) == 2:
        (num_bins,
         num_alpha_dict,
         reverse_read_choice,
         lib_diag_choice,
         path_lib_r1,
         path_lib_r2,
         path_fastq_r1,
         path_fastq_r2,
         MFI_list,
         wt_mfi_exp,
         barcode_size,
         flanking_5p,
         flanking_3p,
         jobname,
         path_dict,
         seed,
         qscore_threshold,
         errors) = parse_config_file(sys.argv[1])
    else:
        print("Error! Please include config file.")


    
    logfile = open(jobname+'_log.txt', 'w')

    overlap_size = 25

    
    print(datetime.now(), file=logfile)
    print('Job name: ', jobname, file=logfile)
    print("Parameters used in this analysis...", file=logfile)
    for facs_bin in range(num_bins):
        print("Path to R1 FASTQ file, bin", num_alpha_dict[facs_bin], ":\n", path_fastq_r1[facs_bin], file=logfile)
        if reverse_read_choice == True:
            print("Path to R2 FASTQ file, bin", num_alpha_dict[facs_bin], ":\n", path_fastq_r2[facs_bin], file=logfile)
        print("Mean fluoresence intensiy, bin", num_alpha_dict[facs_bin], ":", MFI_list[facs_bin], file=logfile)
    print('Mean fluorescence intensity, WT eGFP: ', wt_mfi_exp, file=logfile)
    print('Path to dictionary: ', path_dict, file=logfile)
    print('Barcode length: ', barcode_size, file=logfile)
    print('Flanking sequence upstream of barcode: ', flanking_5p, file=logfile)
    print('Flanking sequence downstream of barcode: ', flanking_3p, file=logfile)

    # Load data
    
    barcode_lists = []
    if reverse_read_choice == True:
        for facs_bin in range(num_bins):
            print('Using R1 and R2 reads from bin', num_alpha_dict[facs_bin], file=logfile)
            r1_raw = fastq_parser(path_fastq_r1[facs_bin])
            r2_raw = fastq_parser(path_fastq_r2[facs_bin])

            r1_readlengths = readlength_counter(r1_raw)
            r2_readlengths = readlength_counter(r2_raw)

            print('R1 read lengths from bin', num_alpha_dict[facs_bin], ':\n', r1_readlengths, file=logfile)
            print('R2 read lengths from bin', num_alpha_dict[facs_bin], ':\n', r2_readlengths, file=logfile)

            filtered_data_dict, avg_scores_r1, avg_scores_r2, passed_r1, passed_r2, raw_reads_r1, raw_reads_r2 = qualityfilter3b(r1_raw, r2_raw, qscore_threshold, errors)

            print('Quality score threshold: ', qscore_threshold, file=logfile)

            print('Number of R1 reads before filtering, bin', num_alpha_dict[facs_bin], ': ', raw_reads_r1, file=logfile)
            print('Number of R2 reads before filtering, bin', num_alpha_dict[facs_bin], ': ', raw_reads_r2, file=logfile)
            print('Percent passing reads, R1 reads, bin', num_alpha_dict[facs_bin], ': ', passed_r1, file=logfile)
            print('Percent passing reads, R2 reads, bin', num_alpha_dict[facs_bin], ': ', passed_r2, file=logfile)
            print('Average quality score, passing R1 reads, bin', num_alpha_dict[facs_bin], ': ', avg_scores_r1, file=logfile)
            print('Average quality score, passing R2 reads, bin', num_alpha_dict[facs_bin], ': ', avg_scores_r2, file=logfile)

            print('Generating consensus sequences from R1 and R2 reads for bin', num_alpha_dict[facs_bin], file=logfile)
            consensus_seqs, both_passed_qscore, overlap_sizes, len_dict, len_filtered = make_consensus_seq(filtered_data_dict, overlap_size)

            print('Number of reads passing either R1 or R2, bin', num_alpha_dict[facs_bin], ': ', len_dict, file=logfile)
            print('Number of reads passing both R1 and R2, bin', num_alpha_dict[facs_bin], ': ', both_passed_qscore, file=logfile)
            print('Number of consensus reads generated with bin', num_alpha_dict[facs_bin], ' R1-R2 overlap of', overlap_size, 'nucleotides: ', len_filtered, file=logfile)


            # bin diagnostic analysis
        
            sequence_start_list_5p, successes_5p, failures_5p, seq_register_count_5p, idx_flanking_seq_5p = flanking_seq_finder(consensus_seqs, flanking_5p)
            sequence_start_list_3p, successes_3p, failures_3p, seq_register_count_3p, idx_flanking_seq_3p = flanking_seq_finder(consensus_seqs, flanking_3p)

            print('Percent of bin', num_alpha_dict[facs_bin], ' consensus sequences containing', flanking_5p, ': ', (successes_5p/(successes_5p+failures_5p)) * 100, file=logfile)
            print('Start positions in bin', num_alpha_dict[facs_bin], ' of ', flanking_5p, ':\n ', seq_register_count_5p, file=logfile)

            print('Percent of bin ', num_alpha_dict[facs_bin], ' consensus sequences containing', flanking_3p, ': ', (successes_3p/(successes_3p+failures_3p)) * 100, file=logfile)
            print('Start positions in bin', num_alpha_dict[facs_bin], ' of ', flanking_3p, ':\n ', seq_register_count_3p, file=logfile)

            barcode_list_bin, unique_barcodes_count, csv_diag_data = isolate_barcodes(consensus_seqs, flanking_5p, flanking_3p, barcode_size)

            print('Number of barcode reads identified in bin', num_alpha_dict[facs_bin], ' consensus sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', len(barcode_list_bin), file=logfile)
            print('Number of unique barcodes identified in bin', num_alpha_dict[facs_bin], ' consensus sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', unique_barcodes_count, file=logfile)

            print_results(csv_diag_data, jobname+'_bin_'+num_alpha_dict[facs_bin]+'_genomic_flanks_diag')
            barcode_lists.append(barcode_list_bin)

    else:
        for facs_bin in range(num_bins):
            print('Using only R1 reads from bin', num_alpha_dict[facs_bin], file=logfile)
            r1_raw = fastq_parser(path_fastq_r1[facs_bin])
            r1_readlengths = readlength_counter(r1_raw)
            print('R1 read lengths from bin', num_alpha_dict[facs_bin], ':\n', r1_readlengths, file=logfile)
            filtered_data, avg_scores_r1, passed_r1, raw_reads_r1 = qualityfilter3_original(r1_raw, qscore_threshold, errors)
            print('Quality score threshold: ', qscore_threshold, file=logfile)
            print('Number of R1 reads before filtering, bin', num_alpha_dict[facs_bin], ': ', raw_reads_r1, file=logfile)
            print('Percent passing reads, R1 reads, bin', num_alpha_dict[facs_bin], ': ', passed_r1, file=logfile)
            print('Average quality score, passing R1 reads, bin', num_alpha_dict[facs_bin], ': ', avg_scores_r1, file=logfile)
        
            print('Omitting consensus sequence generation since only R1 reads are being used for bin', num_alpha_dict[facs_bin], file=logfile)

            sequence_start_list_5p, successes_5p, failures_5p, seq_register_count_5p, idx_flanking_seq_5p = flanking_seq_finder(filtered_data, flanking_5p)
            sequence_start_list_3p, successes_3p, failures_3p, seq_register_count_3p, idx_flanking_seq_3p = flanking_seq_finder(filtered_data, flanking_3p)

            print('Percent of bin', num_alpha_dict[facs_bin], ' passing R1 sequences containing', flanking_5p, ': ', (successes_5p/(successes_5p+failures_5p)) * 100, file=logfile)
            print('Start positions in bin', num_alpha_dict[facs_bin], ' of ', flanking_5p, ':\n ', seq_register_count_5p, file=logfile)

            print('Percent of bin ', num_alpha_dict[facs_bin], ' passing R1 sequences containing', flanking_3p, ': ', (successes_3p/(successes_3p+failures_3p)) * 100, file=logfile)
            print('Start positions in bin', num_alpha_dict[facs_bin], ' of ', flanking_3p, ':\n ', seq_register_count_3p, file=logfile)


            barcode_list_bin, unique_barcodes_count, csv_diag_data = isolate_barcodes(filtered_data, flanking_5p, flanking_3p, barcode_size)

            print('Number of barcode reads identified in bin', num_alpha_dict[facs_bin], ' passing R1 sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', len(barcode_list_bin), file=logfile)
            print('Number of unique barcodes identified in bin', num_alpha_dict[facs_bin], ' passing R1 sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', unique_barcodes_count, file=logfile)

            print_results(csv_diag_data, jobname+'_bin_'+num_alpha_dict[facs_bin]+'_genomic_flanks_diag')
            barcode_lists.append(barcode_list_bin)

    if (reverse_read_choice, lib_diag_choice) == (True, True):
        print('Using R1 and R2 reads from lib file', file=logfile)
        r1_lib_raw = fastq_parser(path_lib_r1)
        r2_lib_raw = fastq_parser(path_lib_r2)
        r1_lib_readlengths = readlength_counter(r1_lib_raw)
        r2_lib_readlengths = readlength_counter(r2_lib_raw)
        print('R1 library file read lengths from bin', num_alpha_dict[facs_bin], ':\n', r1_lib_readlengths, file=logfile)
        print('R2 library file read lengths from bin', num_alpha_dict[facs_bin], ':\n', r2_lib_readlengths, file=logfile)
        filtered_data_dict, avg_scores_r1, avg_scores_r2, passed_r1, passed_r2, raw_reads_r1, raw_reads_r2 = qualityfilter3b(r1_lib_raw, r2_lib_raw, qscore_threshold, errors)
        print('Quality score threshold: ', qscore_threshold, file=logfile)
        print('Number of R1 reads before filtering, lib file: ', raw_reads_r1, file=logfile)
        print('Number of R2 reads before filtering, lib file: ', raw_reads_r2, file=logfile)
        print('Percent passing reads, R1 reads, lib file: ', passed_r1, file=logfile)
        print('Percent passing reads, R2 reads, lib file: ', passed_r2, file=logfile)
        print('Average quality score, passing R1 reads, lib file: ', avg_scores_r1, file=logfile)
        print('Average quality score, passing R2 reads, lib file: ', avg_scores_r2, file=logfile)
        print('Generating consensus sequences from R1 and R2 reads for lib file', file=logfile)
        consensus_seqs_lib, both_passed_qscore, overlap_sizes, len_dict, len_filtered = make_consensus_seq(filtered_data_dict, overlap_size)
        print('Number of reads passing either R1 or R2, lib file: ', len_dict, file=logfile)
        print('Number of reads passing both R1 and R2, lib file: ', both_passed_qscore, file=logfile)
        print('Number of consensus reads generated with lib file R1-R2 overlap of', overlap_size, 'nucleotides: ', len_filtered, file=logfile)
        sequence_start_list_5p, successes_5p, failures_5p, seq_register_count_5p, idx_flanking_seq_5p = flanking_seq_finder(consensus_seqs, flanking_5p)
        sequence_start_list_3p, successes_3p, failures_3p, seq_register_count_3p, idx_flanking_seq_3p = flanking_seq_finder(consensus_seqs, flanking_3p)
        print('Percent of lib file consensus sequences containing', flanking_5p, ': ', (successes_5p/(successes_5p+failures_5p)) * 100, file=logfile)
        print('Start positions in lib file of ', flanking_5p, ':\n ', seq_register_count_5p, file=logfile)
        print('Percent of lib file consensus sequences containing', flanking_3p, ': ', (successes_3p/(successes_3p+failures_3p)) * 100, file=logfile)
        print('Start positions in lib file of ', flanking_3p, ':\n ', seq_register_count_3p, file=logfile)
        barcode_list_lib, unique_barcodes_count_lib, csv_diag_data_lib = isolate_barcodes(consensus_seqs_lib, flanking_5p, flanking_3p, barcode_size)
        print('Number of barcode reads identified in lib file consensus sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', len(barcode_list_lib), file=logfile)
        print('Number of unique barcodes identified in lib file consensus sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', unique_barcodes_count_lib, file=logfile)
        print_results(csv_diag_data_lib, jobname+'_lib_genomic_flanks_diag')


    if (reverse_read_choice, lib_diag_choice) == (False, True):
        r1_lib_raw = fastq_parser(path_lib_r1)
        r1_lib_readlengths = readlength_counter(r1_lib_raw)
        print('R1 library file read lengths from bin', num_alpha_dict[facs_bin], ':\n', r1_lib_readlengths, file=logfile)
        filtered_data_lib, avg_scores_r1, passed_r1, raw_reads_r1 = qualityfilter3_original(r1_lib_raw, qscore_threshold, errors)
        print('Quality score threshold: ', qscore_threshold, file=logfile)
        print('Number of R1 reads before filtering, lib file: ', raw_reads_r1, file=logfile)

        print('Percent passing reads, R1 reads, lib file: ', passed_r1, file=logfile)

        print('Average quality score, passing R1 reads, lib file: ', avg_scores_r1, file=logfile)

        print('Omitting consensus sequences generation since only R1 reads for lib are used', file=logfile)

        sequence_start_list_5p, successes_5p, failures_5p, seq_register_count_5p, idx_flanking_seq_5p = flanking_seq_finder(filtered_data_lib, flanking_5p)
        sequence_start_list_3p, successes_3p, failures_3p, seq_register_count_3p, idx_flanking_seq_3p = flanking_seq_finder(filtered_data_lib, flanking_3p)
        print('Percent of lib file consensus sequences containing', flanking_5p, ': ', (successes_5p/(successes_5p+failures_5p)) * 100, file=logfile)
        print('Start positions in lib file of ', flanking_5p, ':\n ', seq_register_count_5p, file=logfile)
        print('Percent of lib file consensus sequences containing', flanking_3p, ': ', (successes_3p/(successes_3p+failures_3p)) * 100, file=logfile)
        print('Start positions in lib file of ', flanking_3p, ':\n ', seq_register_count_3p, file=logfile)
        barcode_list_lib, unique_barcodes_count_lib, csv_diag_data_lib = isolate_barcodes(filtered_data_lib, flanking_5p, flanking_3p, barcode_size)
        print('Number of barcode reads identified in lib file consensus sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', len(barcode_list_lib), file=logfile)
        print('Number of unique barcodes identified in lib file consensus sequences using 5p ', flanking_5p, 'and 3p ', flanking_3p, ': ', unique_barcodes_count_lib, file=logfile)
        print_results(csv_diag_data, jobname+'_lib_genomic_flanks_diag')


    # Import dictionary data
    variant_barcode_list = import_file(path_dict)
    dict_variant_list = [barcode[1] for barcode in variant_barcode_list]
    set_dict_variants = set(dict_variant_list)
    ordered_list_variants = []
    [ordered_list_variants.append(x) for x in dict_variant_list if x not in ordered_list_variants]
    dict_barcode_list = [barcode[0] for barcode in variant_barcode_list]
    set_dict_barcodes = set(dict_barcode_list)
    variant_barcode_dict = {barcode[0]:barcode[1] for barcode in variant_barcode_list}


    all_barcodes_lib = []
    for facs_bin in range(num_bins):
        all_barcodes_lib += barcode_lists[facs_bin]


    if lib_diag_choice == True:
        all_barcodes_lib += barcode_list_lib
        barcode_count_lib = Counter(barcode_list_lib)

    set_all_barcodes = list(set(all_barcodes_lib))

    print('Number of unique barcodes in all bins: ', len(set_all_barcodes), file=logfile)
    
    barcode_count_list = []
    for facs_bin in range(num_bins):
        barcode_count_bin = Counter(barcode_lists[facs_bin])
        barcode_count_list.append(barcode_count_bin)

    diagnostic_data = []
    header_line = ['Barcode', 'In dict']
    if lib_diag_choice == True:
        header_line.append('Lib')
    for facs_bin in range(num_bins):
        header_line.append(num_alpha_dict[facs_bin])
    diagnostic_data.append(header_line)

    for barcode in set_all_barcodes:
        in_rho_dict = 1 if barcode in set_dict_barcodes else 0
        output_line = [barcode, in_rho_dict]
        if lib_diag_choice == True:
            output_line.append(barcode_count_lib[barcode])
        for facs_bin in range(num_bins):
            bin_count = barcode_count_list[facs_bin][barcode]
            output_line.append(bin_count)
        diagnostic_data.append(output_line)
    print_results(diagnostic_data, jobname+'_diag_combined')


    # use dict to filter out unused barcodes
    filtered_barcodes = []
    for facs_bin in range(num_bins):
        filtered_barcodes_bin = new_illumina_filter(barcode_lists[facs_bin], set_dict_barcodes)
        print('Number of barcodes in bin', num_alpha_dict[facs_bin], 'after filtering with barcodes in dictionary: ', len(filtered_barcodes_bin), file=logfile)
        filtered_barcodes.append(filtered_barcodes_bin)

    # These routines actually search the Illumina reads for labels using the dictionary
    master_reads_pool = []
    for facs_bin in range(num_bins):
        reads_per_mutant_dict_bin, master_pool_bin = analyze_illumina(variant_barcode_dict, filtered_barcodes[facs_bin], ordered_list_variants)
        print_variant_counts(reads_per_mutant_dict_bin, jobname+'_illumina_counts_'+num_alpha_dict[facs_bin])
        master_reads_pool.append(master_pool_bin)

    # This routine converts barcode reads to variant labels
    master_label_reads = []
    for facs_bin in range(num_bins):
        reads_labels_bin, wtcount_bin, varcount_bin, none_count_bin = convert_reads_barcodes_to_labels(variant_barcode_dict, master_reads_pool[facs_bin])
        master_label_reads.append(reads_labels_bin)
        print("Number of barcodes in master pool in bin", num_alpha_dict[facs_bin], ":", len(master_reads_pool[facs_bin]), file=logfile)
        print("Number of wt barcodes in master pool in bin", num_alpha_dict[facs_bin], ":", wtcount_bin, file=logfile)
        print("Number of variant barcodes in master pool in bin", num_alpha_dict[facs_bin], ":", varcount_bin, file=logfile)
        print("Number of unassigned barcodes in master pool in bin", num_alpha_dict[facs_bin], ":", none_count_bin, file=logfile)
        print("Number of reads in master label pool in bin", num_alpha_dict[facs_bin], ":", len(reads_labels_bin), file=logfile)


    len_data = []
    for facs_bin in range(num_bins):
        len_data_bin = len(master_label_reads[facs_bin])
        len_data.append(len_data_bin)
    minimum = min(len_data)
    index_min = np.argmin(len_data)
    print("Bin with the smallest number of reads: ", num_alpha_dict[index_min], file=logfile)
    print("Number of reads in bin", num_alpha_dict[index_min], ":", minimum, file=logfile)

    
    # Perform the rarefactions
    print("Rarefaction sampling...", file=logfile)
    random.seed(seed)
    print("Seed for rarefaction sampling: ", seed, file=logfile)
    rarefaction_01_reads = []
    rarefaction_02_reads = []
    for facs_bin in range(num_bins):
        rarefaction_01_reads_bin = random.sample(master_label_reads[facs_bin], minimum)
        rarefaction_01_reads.append(rarefaction_01_reads_bin)
        rarefaction_02_reads_bin = random.sample(master_label_reads[facs_bin], minimum)
        rarefaction_02_reads.append(rarefaction_02_reads_bin)


    # Variant identification from each rarefaction
    rarefaction_ser_01 = []
    #wtcount_01 = []
    for facs_bin in range(num_bins):
        rarefaction_ser_01_bin = create_pd_series(rarefaction_01_reads[facs_bin], ordered_list_variants)
        rarefaction_ser_01.append(rarefaction_ser_01_bin)
        #wtcount_01.append(wtcount_01_bin)
        
    rarefaction_ser_02 = []
    #wtcount_02 = []
    for facs_bin in range(num_bins):
        rarefaction_ser_02_bin = create_pd_series(rarefaction_02_reads[facs_bin], ordered_list_variants)
        rarefaction_ser_02.append(rarefaction_ser_02_bin)
        #wtcount_02.append(wtcount_02_bin)

    # Get weighted intensities
    print("Calculating surface trafficking score", file=logfile)
    varcount_ints_weight_01, var_wt_norm_01, varcount_dfsum_01 = weighted_average_func(rarefaction_ser_01, MFI_list, wt_mfi_exp)
    # For rarefaction 02:
    varcount_ints_weight_02, var_wt_norm_02, varcount_dfsum_02 = weighted_average_func(rarefaction_ser_02, MFI_list, wt_mfi_exp)


    var_wt_norm_avg = (var_wt_norm_01 + var_wt_norm_02) / 2
    varcount_ints_weight_avg = (varcount_ints_weight_01 + varcount_ints_weight_02) / 2
    
    # Save Excel workbooks
    print("Saving Excel workbooks...", file=logfile)
    # Excel workbook for rarefaction 01
    workbook_counts_rarefaction_01 = pd.ExcelWriter(jobname+'_rarefaction_01.xlsx')
    var_wt_norm_01.to_excel(workbook_counts_rarefaction_01, 'Var wt avg norm to WT')
    varcount_ints_weight_01.to_excel(workbook_counts_rarefaction_01, 'Var wt avg')
    for facs_bin in range(num_bins):
        rarefaction_ser_01[facs_bin].to_excel(workbook_counts_rarefaction_01, 'Rarefied cts from bin {}'.format(num_alpha_dict[facs_bin]))
    varcount_dfsum_01.to_excel(workbook_counts_rarefaction_01, 'Total rarefied cts')
    workbook_counts_rarefaction_01.save()

    # Excel workbook for rarefaction 02
    workbook_counts_rarefaction_02 = pd.ExcelWriter(jobname+'_rarefaction_02.xlsx')
    var_wt_norm_02.to_excel(workbook_counts_rarefaction_02, 'Var wt avg norm to WT')
    varcount_ints_weight_02.to_excel(workbook_counts_rarefaction_02, 'Var wt avg')
    for facs_bin in range(num_bins):
        rarefaction_ser_02[facs_bin].to_excel(workbook_counts_rarefaction_02, 'Rarefied cts from bin {}'.format(num_alpha_dict[facs_bin]))
    varcount_dfsum_02.to_excel(workbook_counts_rarefaction_02, 'Total rarefied cts')
    workbook_counts_rarefaction_02.save()

    # Excel workbook for average of rarefactions, variant count weighted by intensity only
    workbook_var_wt_norm_avg = pd.ExcelWriter(jobname+'_rarefaction_avg.xlsx')
    var_wt_norm_avg.to_excel(workbook_var_wt_norm_avg, 'Var wt avg norm to WT')
    varcount_ints_weight_avg.to_excel(workbook_var_wt_norm_avg, 'Var wt avg')
    workbook_var_wt_norm_avg.save()



    logfile.close()

