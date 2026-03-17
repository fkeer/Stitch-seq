import argparse
from tqdm import tqdm
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from polyleven import levenshtein
import logging
import gzip
from multiprocessing import Pool, cpu_count, current_process
import sys


def correct_barcode(barcode, valid_barcodes, dist_thr, read_name):
    if barcode in valid_barcodes:  # there is an exact match
        return barcode
    # find all the matches within the specified Levenshtein distance
    dist_to_targets = [(levenshtein(barcode, target_barcode), target_barcode) for target_barcode in valid_barcodes]
    matches = [entry[1] for entry in dist_to_targets if entry[0] <= dist_thr]

    if len(matches) == 0 or len(matches) > 1:
        # if no matches are found or if more than one match is found within this distance threshold =>
        # no correction is performed
        logging.info("Could not correct %s in %s, number of matches: %d" % (barcode, read_name, len(matches)))
        return None
    logging.info("Corrected %s to %s in %s" % (barcode, matches[0], read_name))
    return matches[0]

def load_barcode_whitelist(barcode_wl_fname):
    with open(barcode_wl_fname) as f:
        return set([line.strip() for line in f])

def chunk_fastq_iterator(iterator, chunk_size):
    chunk=[]
    for entry in iterator:
        chunk.append(entry)
        if len(chunk) == chunk_size:
            yield chunk
            chunk=[]
    if chunk:
        yield chunk

def process_read_chunk(chunk_data):
    chunk, valid_barcodes, dist_thr, barcode_len = chunk_data

    chunk_n_exact_matches = 0
    chunk_n_total_matches = 0
    output_lines = []

    for read_name, seq, qual in chunk:
        barcode = seq[:barcode_len]
        target_barcode = correct_barcode(barcode, valid_barcodes, dist_thr, read_name)
        if target_barcode:
            chunk_n_exact_matches += target_barcode == barcode
            chunk_n_total_matches += 1
            barcode = target_barcode
        read_name_w_barcode = "BX-%s:%s" % (barcode, read_name)
        output_lines.append("@%s\n%s\n+\n%s\n" % (read_name_w_barcode, seq[barcode_len:], qual[barcode_len:]))
    return output_lines, chunk_n_exact_matches, chunk_n_total_matches


def generate_pool_args(in_handle, valid_barcodes, dist_thr, barcode_len, CHUNK_SIZE):
    """
    Generator function: Yields arguments lazily for the multiprocessing pool.
    """
    # This loop reads the file sequentially, but only yields one chunk at a time.
    for chunk in chunk_fastq_iterator(FastqGeneralIterator(in_handle), CHUNK_SIZE):
        # Yield the argument tuple: (chunk data, constants)
        yield (chunk, valid_barcodes, dist_thr, barcode_len)


def extract_barcodes(input_fastq_fname, output_fastq_fname, barcode_len, valid_barcodes=None, dist_thr=None, threads=cpu_count()):

    n_reads = 0
    n_exact_matches = 0
    n_total_matches = 0

    CHUNK_SIZE=50000

    if input_fastq_fname == "-":
        # Input is being piped from pigz
        in_handle = sys.stdin
    else:
        # Fallback for unpiped testing (will still be slow)
        in_handle = gzip.open(input_fastq_fname, 'rt') 

    # 2. Determine Output Handle (Write to stdout if piped)
    if output_fastq_fname == "-":
        out_handle = sys.stdout
    else:
        out_handle = gzip.open(output_fastq_fname, 'wt') # Or simple open if unzipped


    try:
        
        with Pool(processes=threads) as pool:

            arg_iter = generate_pool_args(in_handle, valid_barcodes, dist_thr, barcode_len, CHUNK_SIZE)
                
            for output_lines, chunk_exact, chunk_total in pool.imap(process_read_chunk, arg_iter):
                n_reads += len(output_lines) 
                n_exact_matches += chunk_exact
                n_total_matches += chunk_total
                    
                out_handle.write("".join(output_lines))

    finally:
        # Ensure file handles are closed
        if in_handle is not sys.stdin:
            in_handle.close()
        if out_handle is not sys.stdout:
            out_handle.close()


    logging.info("# Reads processed: %d" % (n_reads))
    if valid_barcodes:
        logging.info("# gRNAs with an exact match: %d" % (n_exact_matches))
        logging.info("# gRNAs with a match: %d" % (n_total_matches))


def main():
    parser = argparse.ArgumentParser(description='Script for gRNA extraction/correction')
    parser.add_argument('--input_fastq',
                        help='Input FASTQ file where the gRNA is given by the read prefix of specified length')
    parser.add_argument('--output_fastq', help='Output FASTQ file with extracted gRNAs')
    parser.add_argument('--log', help='Log file', default=None)
    parser.add_argument('--grna_len', help='Barcode length', type=int, default=20)
    parser.add_argument('--grna_whitelist', help='gRNA whitelist', default=None)
    parser.add_argument('--dist_thr', help='Maximum Levenshtein distance threshold for gRNA error correction',
                        type=int, default=2)
    parser.add_argument('--threads', help='Number of threads to use for parallel processing',
                        type=int, default=cpu_count()) 
    parser.add_argument('--input_stdin', action='store_true', help='Read input from stdin.')
    parser.add_argument('--output_stdout', action='store_true', help='Write output to stdout.')
    

    args = parser.parse_args()

    handlers = []
    input_file = "-" if args.input_stdin else args.input_fastq
    output_file = "-" if args.output_stdout else args.output_fastq
    if args.log:
        handlers.append(logging.FileHandler(args.log, mode='w'))
    else:
        handlers.append(logging.NullHandler())
    logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.DEBUG, handlers=handlers)
    valid_barcodes = None
    if args.grna_whitelist:
        valid_barcodes = load_barcode_whitelist(args.grna_whitelist)
    extract_barcodes(input_file, output_file, args.grna_len, valid_barcodes, args.dist_thr, args.threads)


if __name__ == '__main__':
    main()
