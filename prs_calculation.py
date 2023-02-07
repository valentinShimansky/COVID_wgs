"""
    Module for prs calculation
"""
import argparse
from subprocess import Popen
import os

PLINK = "plink2"
PRS_SUFFIX = "prs"


def execute(args):
    if not os.path.isfile(args.gwas_file):
        raise Exception(f'File {args.gwas_file} does not exist.')
    model = args.model
    snp_id_col, ea_col, beta_col = args.snp_id, args.ea, args.beta
    bcf_files = args.bcf_file
    output_dir = args.output
    mem = args.mem
    threads = args.threads
    output_file = f'{os.path.splitext(os.path.basename(bcf_files))[0]}_{PRS_SUFFIX}'
    log_file = os.path.join(output_dir, f'{output_file}.log')
    with open(log_file, "a") as log_file:
        command = [PLINK,
                   '--bcf', bcf_files,
                   '--score', model, snp_id_col, ea_col, beta_col, 'header', 'list-variants', 'cols=scoresums',
                   '--out', os.path.join(output_dir, output_file)]
        if mem:
            command += ['--memory', mem]
        if threads:
            command += ['--threads', threads]
        process = Popen(command,
                        stdout=log_file,
                        stderr=log_file)
        process.communicate()


def parse():
    parser = argparse.ArgumentParser()
    input_group = parser.add_argument_group("Input parameters")
    input_group.add_argument('--model', dest="model", type=str, required=True,
                             help='Path to model.')
    input_group.add_argument('--snp-id', dest="snp_id", required=True, type=str,
                             help='Column number for ID of SNP (1-based).')
    input_group.add_argument('--ea', required=True, type=str,
                             help='Column number for effect allele (1-based).')
    input_group.add_argument('--beta', required=True, type=str,
                             help='Column number for beta (1-based).')
    input_group.add_argument('--bcf', dest="bcf_file", type=str, default=None,
                             help='Path to BCF file.')
    input_group.add_argument('--output', dest="output", default=None, type=str,
                             help='Path to result folder.')
    input_group.add_argument('--mem', default=None, type=str, dest='mem',
                             help='Upper bound on the memory pool in megabytes.')
    input_group.add_argument('--threads', default=None, type=str, dest='threads',
                             help='Number of threads for Plink to use.')
    return parser


def run_prs_calculation_main():
    args = parse().parse_args()
    try:
        execute(args)
    except Exception as error:
        raise error


if __name__ == '__main__':
    run_prs_calculation_main()
