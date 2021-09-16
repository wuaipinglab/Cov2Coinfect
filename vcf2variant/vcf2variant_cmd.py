import os
import sys
import argparse

from vcf2variant.vcfparser.vcfparser import vcf2pandas
from vcf2variant.scripts.get_candidate_lineages import get_candidate_lineages
from vcf2variant.scripts.identify_coinfection_lineages import identify_coinfection_lineages
from vcf2variant import __version__


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="VCF to SARS-CoV-2 variant.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"vcf2variant {__version__}")
    parser.add_argument("directory", help="Output directory")
    parser.add_argument("--vcf", default=None,
                        help="Input VCF directory when VCF files are not in the output directory")
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    args = parser.parse_args()

    output_dir = args.directory
    vcf_dir = args.vcf

    data_dir = os.path.join(os.path.dirname(__file__), "data")
    ngs_dir = vcf2pandas(output_dir, vcf_dir)
    candidate_dir = get_candidate_lineages(ngs_dir, output_dir, data_dir)
    identify_coinfection_lineages(candidate_dir, output_dir, data_dir)
