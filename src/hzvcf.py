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

import time

VCF_meta_template = """##fileformat=VCFv4.1
##fileDate={_d.time.tm_year}:{_d.time.tm_mon}:{_d.time.tm_mday}-{_d.time.tm_hour}:{_d.time.tm_min}:{_d.time.tm_sec}
##source=MonoVar_NB
{_d.FILTER_META}
{_d.INFO_META}
{_d.FORMAT_META}
{_d.REF_META}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{_d.FILES_META}
"""

class VCFDocument():

    def __init__(self, outf):
        self.time = time.localtime()
        self.outf = outf

        self.info_fields = [
            ('AC', 'A', 'Integer',
                'Allele count in genotypes, for each ALT allele, in the same order as listed'),
            ('AF', 'A', 'Float',
                'Allele Frequency, for each ALT allele, in the same order as listed'),
            ('AN', '1', 'Integer',
                'Total number of alleles in called genotypes'),
            ('BaseQRankSum', '1', 'Float',
                'Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities'),
            ('DP', '1', 'Integer',
                'Approximate read depth; some reads may have been filtered'),
            ('QD', '1', 'Float',
                'Variant Confidence/Quality by Depth'),
            ('SOR', '1', 'Float',
                'Symmetric Odds Ratio of 2x2 contingency table to detect strand bias'),
            ('MPR', '1', 'Float',
                'Log Odds Ratio of maximum value of probability of observing non-ref allele to the probability of observing zero non-ref allele'),
            ('PSARR', '1', 'Float',
                'Ratio of per-sample Alt allele supporting reads to Ref allele supporting reads')
        ]
        self.INFO_META = '\n'.join([
            '##INFO=<ID={},Number={},Type={},Description="{}">'.format(*i) \
                for i in self.info_fields])

        self.filter_fields = [
            ('LowQual', 'Low quality'),
            ('NoConsensus', 'Only 1 sample contains called SNV')
        ]
        self.FILTER_META = '\n'.join([
            '##FILTER=<ID={},Description="{}">'.format(*i) \
                for i in self.filter_fields])

        self.format_fields = [
            ('AD', '.', 'Integer', 
                'Allelic depths for the ref and alt alleles in the order listed'),
            ('DP', '1', 'Integer',
                'Approximate read depth (reads with MQ=255 or with bad mates are filtered)'),
            ('GQ', '1', 'Integer', 'Genotype Quality'),
            ('GT', '1', 'String', 'Genotype'),
            ('PL', 'G', 'Integer',
                'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification')
        ]
        self.FORMAT_META = '\n'.join([
            '##FORMAT=<ID={},Number={},Type={},Description="{}">'.format(*i) \
                for i in self.format_fields])
        
        self.ref_file = None
        self.REF_META = ''
        self.files_list = []
        self.FILES_META = ''


    def set_files(self, bam_id_list):
        self.files_list.extend(bam_id_list)
        self.FILES_META = '\t'.join(self.files_list)
        

    def set_reference(self, ref_file):
        self.ref_file = ref_file
        self.REF_META = '##reference={}'.format(ref_file)


    def print_header(self):
        self.outf.write(VCF_meta_template.format(_d=self))


    def print_record(self, rec_data):
        rec_str = '\t'.join(rec_data) + '\n'
        self.outf.write(rec_str)


    def close(self):
        self.outf.close()


if __name__ == '__main__':
    print('Here be dragons...')
