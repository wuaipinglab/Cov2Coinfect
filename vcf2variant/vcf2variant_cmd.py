import os
import sys
import argparse

from vcf2variant.vcfparser.vcfparser import vcf2pandas
from vcf2variant import __version__


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="VCF to SARS-CoV-2 variant.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"vcf2variant {__version__}")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("--vcf", default=None,
                        help="Input VCF directory when VCF files are not in the output directory")
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    args = parser.parse_args()

    output_dir = args.output
    vcf_dir = args.vcf

    DIRPATH = os.path.join(os.path.dirname(__file__), "data") + "/"
    csv_dir = vcf2pandas(output_dir, vcf_dir)
