#!/usr/bin/env python
"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Hamim Zafar and Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import sys
import os
import copy_reg
import types
import multiprocessing as mp
import numpy as np
from functools import partial

from utils import Utils_Functions
from mp_genotype import MP_single_cell_genotype
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from calc_variant_prob import Calc_Var_Prob
from hzvcf import VCFDocument, VRecord

# Required for using multiprocessing
def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

U = Utils_Functions()
M = MP_single_cell_genotype()

# Default values
pe = 0.002
pad = 0.2
thr = 0.05
theta = 0.001  # Heterozygosity rate
m_thread = 1
CF_flag = 1
max_depth = 10000
debug = 0

# Process the inputs
argc = len(sys.argv)
i = 1
input_args = {}
while (i < argc):
    if (sys.argv[i] == '-i'):
        pileup = sys.argv[i + 1]         # input pileup file
    elif (sys.argv[i] == '-n'):
        n_cells = int(sys.argv[i + 1])      # Number of input bam files
    elif (sys.argv[i] == '-p'):
        pe = float(sys.argv[i + 1])    	  # probability of error
    elif (sys.argv[i] == '-d'):
        pd = float(sys.argv[i + 1])    	  # probability of deamination error
    elif (sys.argv[i] == '-a'):
        pad = float(sys.argv[i + 1])  	  # Probability of ADO
    elif (sys.argv[i] == '-f'):
        ref_file = sys.argv[i + 1]          # Reference Genome File
        input_args['-f'] = 'Provided'
    elif (sys.argv[i] == '-b'):
        bam_file_list = sys.argv[i + 1]	  # File containing list of bam files
        input_args['-b'] = 'Provided'
    elif (sys.argv[i] == '-o'):
        outfile = sys.argv[i + 1]      	  # Output File
        input_args['-o'] = 'Provided'
    elif (sys.argv[i] == '-t'):
        # threshold to use for calling variant
        thr = float(sys.argv[i + 1])
    elif (sys.argv[i] == '-c'):
        CF_flag = int(sys.argv[i + 1])      # Flag for using Consensus Filter
    elif (sys.argv[i] == '-m'):
        # Number of threads to use in multiprocessing
        m_thread = int(sys.argv[i + 1])
    elif (sys.argv[i] == '-debug'):
        # Number of threads to use in multiprocessing
        debug = int(sys.argv[i + 1])
    i = i + 2

 
try:
    b = input_args['-f']
except KeyError:
    print('Error: Reference genome file not provided. ' \
        'Use "-f"" for reference genome file.\n')
    exit(3)
try:
    b = input_args['-b']
except KeyError:
    print('Error: List of Bam files not provided. ' \
        'Use "-b"" for list of Bam files.\n')
    exit(3)
try:
    b = input_args['-o']
except KeyError:
    print('Error: Output file not provided. Use "-o" for Output file.\n')
    exit(3)
try:
    assert CF_flag <= 1
except AssertionError:
    print('CF_flag can have value 0 or 1. Use "-c" with proper value.\n')
    exit(3)

# Obtain the RG IDs from the bam files
bam_id_list = []
with open(bam_file_list, 'r') as f:
    f_bam_list = f.read().strip('\n').split('\n')
    for f_bam in f_bam_list:
        bam_file = f_bam.strip('\n')
        if not os.path.exists(bam_file):
            bam_file = os.path.join(os.getcwd(), bam_file)
        bam_id_list.append(U.Get_BAM_RG(bam_file))

n_cells = len(bam_id_list)

# Initialize the pool of multiprocessing
pool = mp.Pool(processes=m_thread)

# Constants to be used later
cell_no_threshold = n_cells / 2
# no of possible alternate alleles {0, 1, 2, ..., 2m}
max_allele_cnt = 2 * n_cells + 1
Base_dict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
genotype_dict = {0: '0/0', 1: '0/1', 2: '1/1'}

# Table for all the required nCr
nCr_matrix = U.Create_nCr_mat(max_allele_cnt)

# Dictionary for holding all the priors for different values of n
prior_variant_dict = {i: U.calc_prior(theta, i, 1) for i in range(n_cells + 1)}

# Open VCF file and print header
f_vcf = open(outfile, 'w')
vcf = VCFDocument(f_vcf)
vcf.populate_fields(bam_id_list)
vcf.populate_reference(ref_file)
vcf.print_header()

# List of all single_cell_ftr_pos object
all_single_cell_ftrs_list = n_cells * [None]
# Global list for storing which cell contains read support
read_flag_row = np.zeros(n_cells)
# Global list for storing which cell has alternate allele support
alt_allele_flag_row = np.zeros(n_cells)

if sys.stdin.isatty():
    with open(pileup, 'r') as f:
        lines = f.read().split('\n')
else:
    lines = sys.stdin

for line in lines:
    line = line.strip('\n')
    row = line.split('\t')
    if line == '':
        continue

    contig = row[0]
    pos = int(row[1])
    refBase = U.refineBase(row[2])
    
    total_depth = 0
    total_ref_depth = 0
    for i in range(1, n_cells + 1):
        curr_cell_pos_ftrs = Single_Cell_Ftrs_Pos(refBase, row[3*i: 3*i + 3])
        total_depth += curr_cell_pos_ftrs.depth
        total_ref_depth += curr_cell_pos_ftrs.refDepth
        all_single_cell_ftrs_list[i - 1] = curr_cell_pos_ftrs

    if total_depth <= 10 or total_depth == total_ref_depth:
        continue

    alt_count = total_depth - total_ref_depth
    alt_freq = float(alt_count) / total_depth

    # No reads supporting alternate allele, so no operations needed
    if alt_freq <= 0.01 or refBase not in ['A', 'T', 'G', 'C']:
        continue
    # Cases that are to be prefiltered
    if total_depth > 30 and (alt_count <= 2 or alt_freq <= 0.001):
        continue
    
    # List for storing the sngl_cell_objs that have read support, will be
    # further used in the model
    read_supported_cell_list = []
    # Gloabal list for storing the alternate allele counts
    total_alt_allele_count = np.zeros(4, dtype=int)
    # Global list for storing the indices of the cells having read support
    info_list = ['GT:AD:DP:GQ:PL']

    # Traverse through all the sngl_cell_ftr_obj and if has read support
    # further calculate the other quantities
    for j, sngl_cell_ftr_obj in enumerate(all_single_cell_ftrs_list):
        read_flag = U.checkReadPresence(sngl_cell_ftr_obj)
        if read_flag == 1:
            sngl_cell_ftr_obj.Get_base_call_string_nd_quals(refBase)
            alt_allele_flag = U.CheckAltAllele(sngl_cell_ftr_obj)
            if alt_allele_flag == 1:
                # Update the list of total_alt_allele_count
                total_alt_allele_count += sngl_cell_ftr_obj \
                    .get_Alt_Allele_Count()
            # Populate the list of read supported cells
            read_supported_cell_list.append(sngl_cell_ftr_obj)
        else:
            alt_allele_flag = 0
        read_flag_row[j] = read_flag
        alt_allele_flag_row[j] = alt_allele_flag

    # Operations on the single cells with read support
    # Number of cells with read support
    read_supported_n_cells = len(read_supported_cell_list)
    if read_supported_n_cells == 0:
        continue

    # Number of cells having read support
    read_smpl_count = read_flag_row.sum()
    # Number of cells having alternate allele support
    alt_smpl_count = alt_allele_flag_row.sum()
    # Update alt count
    alt_count = total_alt_allele_count.max()
    # Get the altBase
    if alt_count == 0:
        continue
    altBase = Base_dict[total_alt_allele_count.argmax()]

    # Calculate prior_allele_mat
    prior_allele_mat = U.Get_prior_allele_mat(read_smpl_count, alt_smpl_count, 
        cell_no_threshold, total_depth, alt_freq, pe)

    # Get prior_variant_number distribution (Eq. 11)
    prior_var_no = prior_variant_dict[read_supported_n_cells]

    for cell in read_supported_cell_list:
        cell.store_addl_info(refBase, altBase, alt_freq, prior_allele_mat)

    # Obtain the value of probability of SNV
    var_prob_obj = Calc_Var_Prob(read_supported_cell_list)
    zero_var_prob, denominator = var_prob_obj \
        .calc_zero_var_prob(n_cells, max_depth, nCr_matrix, pad, prior_var_no)

    # Probability of SNV passes the threshold
    if zero_var_prob <= thr:

        func = partial(M.get_info_string, read_supported_cell_list, n_cells,
            nCr_matrix, prior_var_no, denominator, genotype_dict)

        if debug:
            output = [func(i) for i in range(read_supported_n_cells)]
        else:
            output = pool.map(func, range(read_supported_n_cells))
        read_supported_info_list = [p[0] for p in output]
        read_supported_barcodes = [p[1] for p in output]

        barcode = '<'
        for single_cell_ftrs_list in all_single_cell_ftrs_list:
            if single_cell_ftrs_list.depth == 0:
                info_list.append('./.')
                barcode += 'X'
            else:
                info_list.append(read_supported_info_list[0])
                del read_supported_info_list[0]
                barcode += read_supported_barcodes[0]
                del read_supported_barcodes[0]
        barcode += '>'

        if zero_var_prob == 0:
            qual = 99
        else:
            qual = min(99, -10 * np.log10(zero_var_prob))

        AC, AF, AN = U.calc_chr_count(barcode)
        if total_ref_depth > 0:
            baseQranksum = U.calc_base_q_rank_sum(read_supported_cell_list)
        else:
            baseQranksum = 0.0
        QD = U.calc_qual_depth(barcode, all_single_cell_ftrs_list, qual)
        SOR = U.calc_strand_bias(read_supported_cell_list, alt_count)
        max_prob_ratio = U.find_max_prob_ratio(var_prob_obj.matrix)
        PSARR = U.calc_per_smpl_alt_ref_ratio(
            total_ref_depth, alt_count, read_smpl_count, alt_smpl_count)

        # Write line to vcf output
        vcf_record = VRecord(contig, pos)
        info_record = [str(AC), "%.2f" % AF, str(AN), "%.2f" % baseQranksum,
            str(total_depth), str(QD), "%.2f" % SOR,
            "%.2f" % max_prob_ratio, "%.2f" % PSARR]
        info_str = ';'.join(['{}={}'.format(j[0], info_record[i]) \
            for i, j in enumerate(vcf.info_fields)])

        if CF_flag == 1 and U.consensus_filter(barcode):
            filter_str = 'PASS'
        else:
            filter_str = '.'
        vcf_record.set_fields(refBase, altBase, '.', qual, filter_str, info_str,
            info_list, barcode)
        vcf.print_my_record(vcf_record)
