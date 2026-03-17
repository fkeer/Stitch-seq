import argparse
import pysam
from tqdm import tqdm
from collections import defaultdict


def add_grna_tags(input_bam_fname, output_bam_fname, grna_len, counts_report_fname=None, valid_barcodes=None):
    input_bam = pysam.AlignmentFile(input_bam_fname, "rb")
    #output_bam = pysam.AlignmentFile(output_bam_fname, "wb", template=input_bam)
    counts = defaultdict(int)
    for read in tqdm(input_bam):
        # parse read name to extract the gRNA sequence
        # assumed format: BX-[gRNA]....
        grna = read.qname[3:grna_len+3]
        read.qname = read.qname[grna_len+3:]
        read.set_tag('BX', grna)
        #output_bam.write(read)
        if counts_report_fname:
            if not valid_barcodes or grna in valid_barcodes:
                counts[(grna, read.reference_name)] += 1
    input_bam.close()
    #output_bam.close()
    if counts_report_fname:
        with open(counts_report_fname, 'w') as f:
            for k, v in counts.items():
                f.write("%s: %s\n" % (k, v))

def load_barcode_whitelist(barcode_wl_fname):
    with open(barcode_wl_fname) as f:
        return set([line.strip() for line in f])

def main():
    parser = argparse.ArgumentParser(description='Reformats gRNA BX tags in BAM, reports gRNA/gene counts')
    parser.add_argument('--input_bam',
                        help='Input BAM file; gRNAs are given in the read name prefix with format: BX-[gRNA]',
                        required=True)
    parser.add_argument('--output_bam', help='Output BAM file with gRNAs in the BX tag')
    parser.add_argument('--grna_len', help='gRNA length', type=int, default=20)
    parser.add_argument('--counts_report', help='Report file with contig+gRNA counts', default=None)
    parser.add_argument('--grna_whitelist', help='gRNA whitelist to use for count reporting', default=None)
    args = parser.parse_args()
    valid_barcodes = None
    if args.grna_whitelist:
        valid_barcodes = load_barcode_whitelist(args.grna_whitelist)
    add_grna_tags(args.input_bam, args.output_bam, args.grna_len, args.counts_report, valid_barcodes)


if __name__ == '__main__':
    main()
