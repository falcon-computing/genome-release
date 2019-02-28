#!/usr/bin/env python2.7
# Chris Eisenhart chrisEisenhart1992@gmail.com
"""
"""
import sys
import operator
import argparse
import subprocess
import os
import time
import filecmp
import shutil
import gzip
import logging


def load_vcf_file(vcf_file):
    """
    Take a vcf file and return a list of sorted variants handles gzipped and
    regular vcf files
    """
    results = []
    if ".gz" in vcf_file:
        with gzip.open(vcf_file, 'rb') as f:
            for line in f:
                if line.startswith("#"): continue
                results.append(line.strip())
    else:
        for line in open(vcf_file, "r"):
            if line.startswith("#"): continue
            results.append(line.strip())
    return sorted(results)


def compare_vcf_files(file1, file2):
    """
    Load up both vcf files into memory as sorted lists of variants

    Go over each list one variant at a time and compare them.  Flag variants
    which are different. 
    """
    vcf_contents1 = load_vcf_file(file1)
    vcf_contents2 = load_vcf_file(file2)

    logging.info(vcf_contents1, vcf_contents2)
    differing_variants = []
    for var1, var2, in zip(vcf_contents1, vcf_contents2):
        dif = False
        for part1, part2 in zip(var1, var2):
            if part1 is not part2:
                dif = True
        if dif:
            differing_variants.append((var1, var2))
    if differing_variants:
        for pair in differing_variants:
            logging.info("Variants {} differ!".format(" ".join(pair)))
        return False
    return True


def niave_file_compare(file1, file2):
    """
    Do a niave file comparison by checking if they are they within 1% of each
    others file size
    """
    size1 = float(os.stat(file1).st_size[:-1])
    size2 = float(os.stat(file2).st_size[:-1])
    if size1 * .99 <= size2 <= size1 *1.01:
        logging.error("Files {} and {} are too different in size".format(file1, file2))
        return False
    if size2 * .99 <= size1 <= size2 *1.01:
        logging.error("Files {} and {} are too different in size".format(file1, file2))
        return False
    return True


def compare_directories(dir1, dir2):
    """
    Take two directories and compare their contents. 
    """
    # Validate new generated files
    dir1_file_list = []
    for fn in os.listdir(dir1):
        dir1_file_list.append(os.path.join(dir1, fn))
    dir2_file_list = []
    for fn in os.listdir(dir2):
        dir2_file_list.append(os.path.join(dir2, fn))

    failed_list = []
    for file1, file2 in zip(sorted(dir1_file_list), sorted(dir2_file_list)):
       if ".vcf" in file_object and ".tbi" not in file_object:
           compare_vcf_files(file_object, exp_object)
       # do a basic file size comparison here
       elif not niave_file_compare(file_object, exp_object):
           failed_list.append((file_object, exp_object))
    return failed_list


def compare_files_to_expected_files(output_file_list, expected_output_dir):
    """
    Take in two lists of files, go over the lists one index at a time and compare the files
    complain if they are different. 
    """
    failed_list = []
    for file_object in output_file_list:
        exp_object = os.path.join(expected_output_dir, file_object)
        # Differentiate between files and directories
        if os.path.isdir(file_object):
           failed_list += compare_directories(file_object, exp_object)
        else:
           if ".vcf" in file_object and ".tbi" not in file_object:
               compare_vcf_files(file_object, exp_object)
           # do a basic file size comparison here
           elif not niave_file_compare(file_object, exp_object):
               failed_list.append((file_object, exp_object))

    if failed_list:
        for pair in failed_list:
                logging.info("Files {} are not the same! The test has failed".format(pair))
        return False

    logging.info("Output files are as expected")
    return True


def remove_object(file_object):
    """
    File or directory, remove it. 
    """
    if os.path.isdir(file_object):
        shutil.rmtree(file_object)
    else:
        os.remove(file_object)


def run_command(command):
    """
    Take in a comand as a list, run it with the subprocess module, keep the stdout and error
    around and track how long it takes
    """
    logging.info("Running command {}".format(" ".join(command)))
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    start = time.time()
    output, error = process.communicate()
    end = time.time()
    elapsed_time = end - start
    logging.info("Command completed in {} seconds".format(elapsed_time))
    error_code = process.returncode
    return elapsed_time


def test_align(fcs_genome, ref, fastq1, fastq2, out_file):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    paths to two fastq files and an output file path.

    Run the fcs-genome align command (accelerated BWA) on the fastq files and write the output.
    """
    # Run the command and track the time
    command = [fcs_genome, "align", "-r", ref, "-1", fastq1, "-2", fastq2, "-o", out_file]
    elapsed_time = run_command(command)
    return elapsed_time


def test_bqsr(fcs_genome, ref, input_bam, vcf_compare, out_dir, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam file and a VCF file for comparison.
    
    Run the fcs-genome bqsr command on the bam file and generate the output directory.
    """
    # Run the command and track the time
    command = [fcs_genome, "bqsr", "-r", ref, "-i", input_bam, "-K", vcf_compare, "-o", out_dir]
    if use_GATK4:
        command.append("-g")
    elapsed_time = run_command(command)
    return elapsed_time


def test_htc(fcs_genome, ref, bam_dir, out_file, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome htc command (variant calling) on the bam file(s) and write the output.
    """
    # Run the command and track the time
    command = [fcs_genome, "htc", "-r", ref, "-i", bam_dir, "-v", "-o", out_file]
    if use_GATK4:
        command.append("-g")
    elapsed_time = run_command(command)
    return elapsed_time


def test_mutect2(fcs_genome, ref, bam_dir, known_snps, out_file, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome mutect2 command (variant calling) on the bam file(s) and write the output.
    """
    # Run the command and track the time
    command = [fcs_genome, "mutect2", "-r", ref, "-i", bam_dir, "-v", "-o", out_file]
    if use_GATK4:
        command.append("-g")
    elapsed_time = run_command(command)
    return elapsed_time


def test_joint(fcs_genome, ref, bam_dir, known_snps, out_file, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome joint command (variant calling) on the bam file(s) and write the output.
    """
    # Run the command and track the time
    command = [fcs_genome, "joint", "-r", ref, "-i", bam_dir, "-v", "-o", out_file]
    if use_GATK4:
        command.append("-g")
    elapsed_time = run_command(command)
    return elapsed_time


def test_ug(fcs_genome, ref, bam_dir, known_snps, out_file, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome ug command (variant calling) on the bam file(s) and write the output.
    """
    command = [fcs_genome, "ug", "-r", ref, "-i", bam_dir, "-v", "-o", out_file]
    if use_GATK4:
        command.append("-g")
    elapsed_time = run_command(command)
    return elapsed_time


def test_germline(fcs_genome, ref, fastq1, fastq2, known_snps, out_file, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome germline command (variant calling) on the bam file(s) and write the output.
    """
    command = [fcs_genome, "germline", "-r", ref, "-1", fastq1, "-2", fastq2, "-v", "-o", out_file]
    if use_GATK4:
        command.append("-g")
    elapsed_time = run_command(command)
    return elapsed_time
