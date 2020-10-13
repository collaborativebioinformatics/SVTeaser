import os
import math
import shutil
import inspect
import logging
import argparse
import subprocess
from collections import OrderedDict
from random import randint

from acebinf import cmd_exe
from truvari import setup_logging
from svteaser.vcfeditor import update_vcf
from svteaser.utils import vcf_compress, add_fasta_entry
import pandas as pd
import pysam

def edit_surv_params(fn):
    """
    Edit the SURVIVOR simSV parameters file
    """
    params = OrderedDict()
    with open(fn, 'r') as fh:
        for line in fh:
            data = line.strip().split(":")
            name = data[0]
            try:
                val = int(data[1].strip())
            except ValueError:
                val = data[1].strip()
            params[name] = val
    
    # We can control things here. 
    # might make sense to expose all of these to command line?
    params["TRANSLOCATION_number"] = 0
    params["INVERSION_number"] = 0
    params["DUPLICATION_number"] = 0
    params["INV_del_number"] = 0
    params["INV_dup_number"] = 0
    params["INDEL_minimum_length"] = 50
    params["INDEL_maximum_length"] = 1000
    params["DUPLICATION_maximum_length"] = 1000
    with open(fn, 'w') as fout:
        for key, val in params.items():
            fout.write(f"{key}: {val}\n")

def generate_surv_params(param_file):
    logging.debug(f"Running SURVIVOR")
    ret = cmd_exe("SURVIVOR simSV {}".format(param_file))
    logging.debug(ret.stderr)
    logging.debug(ret.stdout)
    if ret.ret_code != 0:
        logging.error("Problem running SURVIVOR")
        logging.error(ret.stderr)
        exit(ret.ret_code)

def generate_regions_from_file(regions_file):
    regions = pd.read_csv(regions_file, sep=", ")
    region_list = []
    for index, row in regions.iterrows():
        # 1. Fetch the region from reference
        chrm = row['chr']
        start = row['start']
        end = row['end']
        region_list.append((chrm, start, end))
    return region_list

def verify_requested_regions(ref, num_regions, length):
    total_regions = 0
    for chrom in ref.references:
        chrom_length = ref.get_reference_length(chrom)
        total_regions = total_regions + math.floor(chrom_length/length)

    if total_regions < num_regions:
        logging.critical("Unable to generate %d non-overlapping regions. Reference too short?", num_regions)
        logging.critical("Generating %d non-overlapping regions", total_regions)
        return total_regions
    else:
        return num_regions

def generate_random_regions(ref_file, region_length, num_regions):
    def generate_region(ref, length):
        chridx = randint(0, len(ref.references)-1)
        chrom = ref.references[chridx]
        chrom_length = ref.get_reference_length(chrom)
        num_non_overlap_regions = math.floor(chrom_length/length)
        randidx = randint(0,num_non_overlap_regions)

        start = randidx*length
        end = start + length
        return randidx, chrom, start, end

    ref = pysam.Fastafile(ref_file)

    region_list = []
    chrom_randidx = {}

    num_regions = verify_requested_regions(ref, num_regions, region_length)

    # We let the while loop run for num_tries before exiting. 
    num_tries = 2*num_regions
    loop_count = 0
    while len(region_list) != num_regions:
        loop_count = loop_count + 1
        randidx, chrom, start, end = generate_region(ref, region_length) 

        # If the region contains "N", then discard this turn.
        reg_string = ref.fetch(chrom, start, end)
        if "N" in reg_string:
            continue

        if chrom in chrom_randidx:
            if randidx not in chrom_randidx[chrom]:
                region_list.append((chrom, start, end))
                chrom_randidx[chrom].append(randidx)
            if loop_count == num_tries:
                logging.critical("Unable to generate %d non-overlapping regions. Tried %d times", (num_regions, num_tries))
                logging.error("Exiting")
                exit(1)
        else:
            region_list.append((chrom, start, end))
            chrom_randidx[chrom] = [randidx]

    return region_list

def update_altered_fa(ref_seq, altered_ref_seq, padding):
    begin_seq = ref_seq[0:padding]
    end_seq = ref_seq[len(ref_seq)-padding: ]
    return begin_seq + altered_ref_seq + end_seq


def process_regions(ref_file, regions, out_dir, param_file):
    out_vcf_path = os.path.join(out_dir, "svteaser.sim.vcf")
    out_ref_fa_path = os.path.join(out_dir, "svteaser.ref.fa")
    out_altered_fa_path = os.path.join(out_dir, "svteaser.altered.fa")
    from io import StringIO
    out_vcf_fh = StringIO()

    header = None
    out_ref_fa_fh = open(out_ref_fa_path, "w+")
    out_altered_fa_fh = open(out_altered_fa_path, "w+")
    chr_header = None
    ref = pysam.FastaFile(ref_file)

    # Define padding in reference region where SVs are not to be inserted.
    padding = 800
    logging.debug("Processing regions")
    for i, (chrom, start, end) in enumerate(regions):
        logging.debug("%s %d %s", chrom, start, end)
        # Track status.
        if (i + 1) % 50 == 0:
            logging.info("Processed {}/{} regions...".format(i + 1, len(regions)))

        # Temporary dir.
        temp_dir = os.path.join(out_dir, "temp")
        os.mkdir(temp_dir)

        # Extract ref sequence.
        name = "{}_{}_{}".format(chrom, start, end)
        ref_seq = ref.fetch(chrom, start, end)

        # Remove some buffer from beginning and ending,
        # so that the tails do not contain SVs. These will be added
        # back later on.
        ref_seq_surv = ref_seq[padding:len(ref_seq)-padding]
        # Write ref sequence to temporary fa file.
        temp_ref_fa = os.path.join(temp_dir, "temp_ref.fa")
        with open(temp_ref_fa, "w") as fh:
            add_fasta_entry(name, ref_seq_surv, fh)

        # Run SURVIVOR.
        prefix = os.path.join(temp_dir, "simulated")
        ret = cmd_exe("SURVIVOR -h")
        survivor_cmd = " ".join(["SURVIVOR",
                                 "simSV",
                                 temp_ref_fa,
                                 param_file,
                                 "0.0",
                                 "0",
                                 prefix])
        for i in range(100):
            ret = cmd_exe(survivor_cmd)
        ret = cmd_exe(survivor_cmd)
        # should be checking here

        # Read output of SURVIVOR
        altered_fa_path = "{}.fasta".format(prefix)
        insertions_fa_path = "{}.insertions.fa".format(prefix)
        sim_vcf = "{}.vcf".format(prefix)
        # Update VCF
        temp_vcf = os.path.join(temp_dir, "temp.vcf")
        update_vcf(temp_ref_fa, insertions_fa_path, sim_vcf, temp_vcf, pos_padding=padding)

        # Merge seqs and variants entries into single FA/VCF files
        # Add the initial and last 800bp back to the altered fasta
        altered_seq = pysam.FastaFile(altered_fa_path).fetch(name)
        altered_seq = update_altered_fa(ref_seq, altered_seq, padding)
        add_fasta_entry(name, altered_seq, out_altered_fa_fh)

        add_fasta_entry(name, ref_seq, out_ref_fa_fh)

        vcf_reader = pysam.VariantFile(temp_vcf)
        if not header:
            header = str(vcf_reader.header).split('\n')
            chr_header = header[-2]
            header = header[:-2]
        else:
            for name, ctg in vcf_reader.header.contigs.items():
                header.append(f"##contig=<ID={name},length={ctg.length}>")

        for record in vcf_reader:
            out_vcf_fh.write(str(record))


        # Remove temporary files.
        shutil.rmtree(temp_dir)
    
    with open(out_vcf_path, 'w') as fout:
        fout.write("\n".join(header) + '\n' + chr_header + '\n')
        out_vcf_fh.seek(0)
        fout.write(out_vcf_fh.read())
    out_altered_fa_fh.close()
    out_ref_fa_fh.close()
    out_vcf_fh.close()
    vcf_compress(out_vcf_path)

def find_survivor():
    ret = cmd_exe("SURVIVOR -h")
    if ret.ret_code != 0:
        logging.error("Cannot find SURVIVOR in environment")
        exit(ret.ret_code)

def surv_sim_main(args):
    """
    Run survivor simSV commands. Output them into a directory
    """
    args = parseArgs(args)
    # check the SURVIVOR is in the environment
    find_survivor()

    try:
        os.mkdir(args.output)
    except FileExistsError:
        logging.error(f"Output directory {args.output} already exists")
        exit(1)

    # Generate SURVIVOR param file
    param_file = os.path.join(args.output, "surv_params")
    generate_surv_params(param_file)
    edit_surv_params(param_file)

    regions = None
    if args.sv_regions:
        # Read sv_regions file, if provided
        regions = generate_regions_from_file(args.sv_regions)
    elif args.num_sv_regions:
        #Choose a random chromosome, a random region of 10kb within the chromosome
        regions = generate_random_regions(args.reference,
                                          args.len_sv_region,
                                          args.num_sv_regions)

    assert(regions is not None), "No regions to process. Please provide at least 1 region."

    process_regions(args.reference, regions, args.output, param_file)

    logging.info("Finished")


def parseArgs(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="surv_sim", description=inspect.getdoc(surv_sim_main),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("reference", metavar="REF", type=str,
                        help="Reference file overwhich to simulate SVs")
    parser.add_argument("output", metavar="OUT", type=str, default="output",
                        help="SVTeaser output basename (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    parser.add_argument('--sv_regions', type=str,
                        help='Comma separated file containing (chr, region_start, region_end). \
                        For every row, an SV of length randint(50, mx_variation) is generated with the region \
                        specified by (chr, start, end).\
                        chr, start, end \
                        chr22, 1000, 20000 \
                        chr22, 50000, 80000', required=False)
    parser.add_argument('--num_sv_regions', type=int, default=10,
                        help='Alternatively to the csv file defined by --sv_regions, user can also \
                              provide number of regions to generate SVs for. The programme will randomly \
                              choose locations within the genome to introduce the SVs. --sv_regions will be given priority \
                              if both options are provided.',
                        required=False)
    parser.add_argument('--len_sv_region', type=int, default=10000,
                        help='The length of regions to create.',
                        required=False)
    args = parser.parse_args(args)
    args.reference = os.path.abspath(args.reference)
    args.output = args.output + ".svt"
    setup_logging(args.debug)
    return args
