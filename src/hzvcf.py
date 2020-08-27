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

import argparse
import os
import sys
import gzip
from subprocess import Popen, PIPE, check_call
import re
import time

VCF_meta_template = """##fileformat=VCFv4.1
##fileDate={_t.tm_year}-{_t.tm_mon}-{_t.tm_mday}
##source=MonoVar
{_d.FILTER_META}{_d.INFO_META}{_d.FORMAT_META}{_d.REF_META}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{_d.FILES_META}
"""

VCF_record_template = "{_r.CHROM}\t{_r.POS}\t{_r.ID}\t{_r.REF}\t{_r.ALT}\t{_r.QUAL}\t{_r.FILTER}\t{_r.INFO}\t{_r.FORMAT}\t{_r.PASSCODE}\n"


class VCFDocument():

    def __init__(self, outf):

        self.time = time.ctime()

        self.info_fields = []
        self.filter_fields = []
        self.format_fields = []
        self.outf = outf

    def populate_fields(self, bam_id_list):
        self.filter_fields.append(('LowQual', 'Low quality'))
        self.format_fields.append(
            ('AD', '.', 'Integer', 'Allelic depths for the ref and alt alleles in the order listed'))
        self.format_fields.append(
            ('DP', '1', 'Integer', 'Approximate read depth (reads with MQ=255 or with bad mates are filtered)'))
        self.format_fields.append(('GQ', '1', 'Integer', 'Genotype Quality'))
        self.format_fields.append(('GT', '1', 'String', 'Genotype'))
        self.format_fields.append(
            ('PL', 'G', 'Integer', 'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification'))
        self.info_fields.append(
            ('AC', 'A', 'Integer', 'Allele count in genotypes, for each ALT allele, in the same order as listed'))
        self.info_fields.append(
            ('AF', 'A', 'Float', 'Allele Frequency, for each ALT allele, in the same order as listed'))
        self.info_fields.append(
            ('AN', '1', 'Integer', 'Total number of alleles in called genotypes'))
        self.info_fields.append(
            ('BaseQRankSum', '1', 'Float', 'Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities'))
        self.info_fields.append(
            ('DP', '1', 'Integer', 'Approximate read depth; some reads may have been filtered'))
        self.info_fields.append(
            ('QD', '1', 'Float', 'Variant Confidence/Quality by Depth'))
        self.info_fields.append(
            ('SOR', '1', 'Float', 'Symmetric Odds Ratio of 2x2 contingency table to detect strand bias'))
        self.info_fields.append(
            ('MPR', '1', 'Float', 'Log Odds Ratio of maximum value of probability of observing non-ref allele to the probability of observing zero non-ref allele'))
        self.info_fields.append(
            ('PSARR', '1', 'Float', 'Ratio of per-sample Alt allele supporting reads to Ref allele supporting reads'))
        self.files_list = bam_id_list

    def populate_reference(self, ref_file):
        self.ref_file = ref_file

    def print_header(self):

        self.FILTER_META = ''.join(
            '##FILTER=<ID=%s,Description="%s">\n' % _ for _ in self.filter_fields)
        self.FORMAT_META = ''.join(
            '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % _ for _ in self.format_fields)
        self.INFO_META = ''.join(
            '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % _ for _ in self.info_fields)
        self.FILES_META = '\t'.join(self.files_list)
        self.REF_META = '##reference=file:{0}\n'.format(self.ref_file)
        self.outf.write(VCF_meta_template.format(_d=self, _t=time.localtime()))

    def print_record(self, record):

        record.INFO = ';'.join("%s=%s" % (
            _[0], str(record.info[_[0]])) for _ in self.info_fields if _[0] in record.info)
        self.outf.write(VCF_record_template.format(_r=record))

    def print_my_record(self, record):
        self.outf.write(VCF_record_template.format(_r=record))

    def close(self):

        self.outf.close()


class VRecord:
    def __init__(self, chrm, pos):
        self.CHROM = chrm
        self.POS = pos


    def set_fields(self, ref, alt, id_in, qual, filter_in, info, samples,
                barcode):
        self.ID = id_in
        self.REF = ref
        self.ALT = alt
        self.QUAL = str(qual)
        self.FILTER = filter_in
        self.INFO = info
        self.FORMAT = "\t".join(samples)
        self.PASSCODE = barcode


if __name__ == '__main__':
    print('Here be dragons...')
