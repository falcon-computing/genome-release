#!/usr/bin/env python2.7
# Chris Eisenhart chrisEisenhart1992@gmail.com
"""
Test the falcon toolkit, runs all pipelines on the input files.

The user can provide an output directory which has expected output files that will be 
compared to the generated files. 
"""
import sys
import operator
import argparse
import subprocess
import os
import falconCommonTest
import logging


def parseArgs(args): 
    """ 
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using args.option. For example args.a stores the alphabet provided. 
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument ("build",
            help = "Path to fcs-genome toolkit that should be tested",
            action = "store")
    parser.add_argument ("reference_genome",
            help = "Path to a reference.fasta file, it is assumed necessary indeces are there (.fai, .bwt, etc)",
            action = "store")
    parser.add_argument ("expected_output_dir",
            help = "Directory that holds expected output files the generated ones should be compared too",
            action = "store")
    parser.add_argument ("fastq_R1",
            help = "The R1 fastq file",
            action = "store")
    parser.add_argument ("fastq_R2",
            help = "The R2 fastq file",
            action = "store")
    parser.add_argument ("known_SNPs",
            help = "Known variants that should be used for base recalibration (vcf file)",
            action = "store")
    # Optional variables
    parser.add_argument("--verbose",
            # Sets warning tags to on
            help = "The program will send warning messages to stderr. This is very messy, use with caution.",
            action= "store_false")
    parser.add_argument("--no_comparison",
            help = "Do not compare to the expected output files",
            action= "store_true")
    parser.add_argument("--hard_remove_files",
            help = "Remove generated files even if they fail comparison",
            action= "store_true")

    parser.set_defaults(no_comparison = False)
    args = parser.parse_args() # Initializes the options specified above
    return args # Returns the parser as an object


def create_generated_file_dict(generic_fn):
    """
    Over the course of testing we will create many files. These files are tracked with a
    dictionary which is put together here.
    """
    bam_file = generic_fn + ".bam"
    bam_index_file = generic_fn + ".bam.bai"

    # GATK3 variant files
    bqsr_dir = generic_fn
    htc_vcf_file = generic_fn + ".htc.vcf.gz"
    htc_vcf_index_file = generic_fn + ".htc.vcf.gz.tbi"
    mutect2_vcf_file = generic_fn + ".mutect2.vcf.gz"
    joint_vcf_file = generic_fn + ".joint.vcf.gz"
    ug_vcf_file = generic_fn + ".ug.vcf.gz"
    germline_vcf_file = generic_fn + ".germline.vcf.gz"

    # GATK4 variant files
    generic_4_fn = generic_fn + "_gatk4"
    bqsr_4_dir = generic_4_fn
    htc_vcf_4_file = generic_4_fn + ".htc.vcf.gz"
    htc_vcf_4_index_file = generic_4_fn + ".htc.vcf.gz.tbi"
    mutect2_vcf_4_file = generic_4_fn + ".mutect2.vcf.gz"
    joint_vcf_4_file = generic_4_fn + ".joint.vcf.gz"
    ug_vcf_4_file = generic_4_fn + ".ug.vcf.gz"
    germline_vcf_4_file = generic_4_fn + ".germline.vcf.gz"

    # Make a dict of all files that will be generated (for comparisons and to remove them later)
    generated_file_dict = {"bam_file":bam_file, "bam_index_file":bam_index_file, "bqsr_dir":bqsr_dir,
                           "htc_vcf_file":htc_vcf_file, "htc_vcf_index_file":htc_vcf_index_file,
                           "mutect2_vcf_file":mutect2_vcf_file, "joint_vcf_file":joint_vcf_file,
                           "ug_vcf_file":ug_vcf_file, "germline_vcf_file":germline_vcf_file,
                           "bqsr_4_dir":bqsr_4_dir, "htc_vcf_4_file":htc_vcf_4_file,
                           "htc_vcf_4_index_file":htc_vcf_4_index_file, "mutect2_vcf_4_file":mutect2_vcf_4_file,
                           "joint_vcf_4_file":joint_vcf_4_file, "ug_vcf_4_file":ug_vcf_4_file,
                           "germline_vcf_4_file":germline_vcf_4_file}


def run_pipeline_components(fcs_genome, ref, fastq1, fastq2, gen_fdict):
    """
    Take in a path to the binary, reference genome path, two fastq files and a dictionary mapping
    string file names to their path.

    gen_fdict = generated file dictionary; maps string file names to their paths
    """
    align_time = falconCommonTest.test_align(fcs_genome, ref, fastq1, fastq2, gen_fdict["bam_file"])

    # GATK3
    bqsr_time = falconCommonTest.test_bqsr(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["bqsr_dir"])
    htc_time = falconCommonTest.test_htc(fcs_genome, ref, bqsr_dir, gen_fdict["htc_vcf_file"])
    mutect2_time = falconCommonTest.test_mutect2(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["mutect2_vcf_file"])
    joint_time = falconCommonTest.test_joint(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["joint_vcf_file"])
    ug_time = falconCommonTest.test_ug(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["ug_vcf_file"])
    # Falcon germline NOTE: this isn't working for the local testing
    germline_time = falconCommonTest.test_germline(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["germline_vcf_4_file"])

    # Run GATK4, it uses the same .bam file
    bqsr_4_time = falconCommonTest.test_bqsr(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["bqsr_dir"], use_GATK4=True)
    htc_4_time = falconCommonTest.test_htc(fcs_genome, ref, bqsr_4_dir, gen_fdict["htc_vcf_file"],  use_GATK4=True)
    mutect2_time = falconCommonTest.test_mutect2(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["mutect2_vcf_file"], use_GATK4=True)
    joint_time = falconCommonTest.test_joint(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["joint_vcf_file"], use_GATK4=True)
    ug_time = falconCommonTest.test_ug(fcs_genome, ref, gen_fdict["bam_file"], known_snps, ["ug_vcf_file"], use_GATK4=True)
    germline_time = falconCommonTest.test_germline(fcs_genome, ref, gen_fdict["bam_file"], known_snps, gen_fdict["germline_vcf_4_file"], use_GATK4=True)

    logging.info("Run time table;\nalign\tbqsr\thtc\tmutect2\tjoint\tug\tgermline\tbqsr_4\thtc4\tmutect2_4\tjoint_4\tug_4\tgermline_4\n" +
                 "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(align_time, bqsr_time, htc_time, mutect2_time, joint_time, ug_time, germline_time) +
                 "\t{}\t{}\t{}\t{}\t{}\t{}\n".format(bqsr_4_time, htc_4_time, mutect2_4_time, joint_4_time, ug_4_time, germline_4_time))


def compare_generated_files_to_expected(no_comparison, generated_file_list, expected_output_dir):
    """
    Take in a list of files that the pipeline generated and a path to a directory holding
    expected output files. It is assumed the files are named the exact same. 

    Compare the generated files to the expected files, complain is something is off.
    """
    result = False
    if not no_comparison:
        try:
            passed_comparison = falconCommonTest.compare_files(generated_file_list, options.expected_output_dir)
            result = passed_comparison
        except:
            logging.exception("Tragic failure in comparison")
    return result


def remove_files(passed_comparison, generated_file_list):
    """
    Remove files and directories if passed_comparison is True
    """
    if passed_comparison:
        for new_file in generated_file_list:
           falconCommonTest.remove_object(new_file)


def main(args):
    """
    This is mostly a wrapper script for functions found in falconCommonTest.  This script
    runs the best practices pipeline on a single sample and tracks the run times. 
    """
    options = parseArgs(args)

    # Input files
    fcs_genome = options.build
    ref = options.reference_genome
    known_snps = options.known_SNPs
    fastq1 = options.fastq_R1
    fastq2 = options.fastq_R2

    # Create a dictionary that maps variable names to a file name. All these files should
    # be generated by testing and need to be compared and removed later.
    generic_fn = fastq1.replace(".gz", "").replace(".fastq", "").split("/")[-1].split("_")[0]
    gen_fdict = create_generated_file_dict(generic_fn)

    # Run the various components of the pipeline and track their run times.
    run_pipeline_components(fcs_genome, ref, fastq1, fastq2, gen_fdict["bam_file"], gen_fdict)

    # Compare the generated files to expected ones.
    passed_comparison = compare_generated_files_to_expected(options.no_comparison, gen_fdict.values(), options.expected_output_dir)

    # Remove the files if they passed the comparison
    if options.hard_remove_files: passed_comparison = True
    remove_files(passed_comparison, gen_fdict.values())


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
