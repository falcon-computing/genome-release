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
    try:
        vcf_contents1 = load_vcf_file(file1)
    except IOError:
        logging.warning("The output file {} is not present".format(file1))
        return False
    except:
        logging.exception("Comparing vcf files {} failed".format(file1))
        return False
    try:
        vcf_contents2 = load_vcf_file(file2)
    except IOError:
        logging.warning("The expected output file {} is not present".format(file2))
        return False
    except:
        logging.exception("Comparing vcf files {} failed".format(file2))
        return False
    differing_variants = []
    for var1, var2, in zip(vcf_contents1, vcf_contents2):
        dif = False
        if var1 != var2: dif = True
        if dif:
            differing_variants.append((var1, var2))
    if differing_variants:
        for pair in differing_variants:
            logging.warning("Variants {} differ!".format(" ".join(pair)))
        return False
    return True


def niave_file_compare(file1, file2):
    """
    Do a niave file comparison by checking if they are they within 2% of each
    others file size
    """
    try:
        size1 = float(os.stat(file1).st_size)
    except IOError:
        logging.warning("The output file {} is not present".format(file1))
        return False
    except OSError:
        logging.warning("The output file {} is not present".format(file1))
        return False
    except:
        logging.exception("Comparing files {} failed".format(file1))
        return False
    try:
        size2 = float(os.stat(file2).st_size)
    except IOError:
        logging.warning("The output file {} is not present".format(file2))
        return False
    except OSError:
        logging.warning("The output file {} is not present".format(file2))
        return False
    except:
        logging.exception("Comparing files {} failed".format(file2))
        return False
    lower_size = size1 * .98
    upper_size = size1 * 1.02
    if upper_size <= size2 <= lower_size:
        logging.error("Files {} and {} are too different in size".format(file1, file2))
        return False
    return True


def compare_directories(dir1, dir2):
    """
    Take two directories and compare their contents. 
    """
    # Validate new generated files
    dir1_file_list = []
    try:
        for fn in os.listdir(dir1):
            dir1_file_list.append(os.path.join(dir1, fn))
    except OSError:
        logging.warning("Could not find dir {}".format(dir1))
        return []
    dir2_file_list = []
    try:
        for fn in os.listdir(dir2):
            dir2_file_list.append(os.path.join(dir2, fn))
    except OSError:
        logging.warning("Could not find dir {}".format(dir1))
        return []

    failed_list = []
    for file1, file2 in zip(sorted(dir1_file_list), sorted(dir2_file_list)):
       if ".vcf" in file1 and ".tbi" not in file1:
           compare_vcf_files(file1, file2)
       # do a basic file size comparison here
       elif not niave_file_compare(file1, file2):
           failed_list.append((file1, file2))
    return failed_list


def compare_files_to_expected_files(output_file_list, expected_output_dir):
    """
    Take in one list of files and an output directory, go over the list one file at a time
    verify a file with the same name exists in the expected output directory and do a comparison.
    """
    failed_list = []
    for file_object in output_file_list:
        exp_object = os.path.join(expected_output_dir, file_object)
        # Differentiate between files and directories
        if os.path.isdir(file_object):
           failed_list += compare_directories(file_object, exp_object)
        else:
            if "vcf" in file_object and "tbi" not in file_object: 
               compare_vcf_files(file_object, exp_object)
           # do a basic file size comparison here
            elif not niave_file_compare(file_object, exp_object):
               failed_list.append((file_object, exp_object))

    if failed_list:
        for pair in failed_list:
            logging.error("Files {} are not the same! The test has failed".format(pair))
        return False

    logging.debug("Output files are as expected")
    return True


def remove_object(file_object):
    """
    File or directory, remove it. 
    """
    try:
        if os.path.isdir(file_object):
            
            shutil.rmtree(file_object)
            logging.debug("Removed directory {}".format(file_object))
        else:
            os.remove(file_object)
            logging.debug("Removed file {}".format(file_object))
    except IOError:
        logging.debug("The file {} is not present".format(file_object))
    except OSError:
        logging.debug("The output file {} is not present".format(file_object))
    except:
        logging.exception("Could not remove {}".format(file_object))


def run_command(command, timeout, expected_dir, out_file_list):
    """
    Take in a command as a list, run it with the subprocess module, keep the stdout and error
    around and track how long it takes
    """
    logging.debug("Running command {}".format(" ".join(command)))
    for new_file in out_file_list:
        remove_object(new_file)


    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        start = time.time()
        timer = 0
        while process.poll() is None:
            time.sleep(1)
            timer += 1
            sys.stdout.flush    
            if timer == timeout:
                raise Exception ("Command timed out")
        output, error = process.communicate()
    except:
        logging.exception("Command {} failed.".format(" ".join(command)))
    end = time.time()
    elapsed_time = '%.2f' % float(end - start) # Round to 2 decimal places
    if int(end - start) >= timeout: elapsed_time = "timed_out"
    logging.debug("Command completed in {} seconds".format(elapsed_time))
    error_code = process.returncode
    if error_code > 0:
        logging.error("{}".format(error))
        elapsed_time = "failed"

    if elapsed_time is not "failed":
        passed_comparison = compare_files_to_expected_files(out_file_list, expected_dir)
        if not passed_comparison:
            elapsed_time = "mangled"

    for new_file in out_file_list:
        logging.debug("Removing {}".format(new_file))
        remove_object(new_file)

    return elapsed_time


def test_align(fcs_genome, timeout, expected_dir, ref, fastq1, fastq2, out_file):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    paths to two fastq files and an output file path.

    Run the fcs-genome align command (accelerated BWA) on the fastq files and write the output.
    """
    # Run the command and track the time
    command = [fcs_genome, "align", "-r", ref, "-1", fastq1, "-2", fastq2, "-o", out_file]
    out_file_list = [out_file, out_file + ".bai"]
    elapsed_time = run_command(command, timeout, expected_dir, out_file_list)
    return elapsed_time


def test_bqsr(fcs_genome, timeout, expected_dir, ref, vcf_compare, generic_fn, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam file and a VCF file for comparison.
    
    Run the fcs-genome bqsr command on the bam file and generate the output directory.
    """
    # Run the command and track the time
    input_bam = os.path.join(expected_dir, generic_fn + ".bam")
    command = [fcs_genome, "bqsr", "-r", ref, "-i", input_bam, "-K", vcf_compare, "-o", generic_fn]
    if use_GATK4:
        command.append("-g")
    out_file_list = [generic_fn]
    elapsed_time = run_command(command, timeout, expected_dir, out_file_list)
    return elapsed_time


def test_htc(fcs_genome, timeout, expected_dir, ref, generic_fn, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome htc command (variant calling) on the bam file(s) and write the output.
    """
    # Run the command and track the time
    bam_dir = os.path.join(expected_dir, generic_fn)
    htc_vcf_file = generic_fn + ".htc.vcf"
    htc_vcf_index_file = generic_fn + ".htc.vcf.gz.tbi"
    command = [fcs_genome, "htc", "-r", ref, "-i", bam_dir, "-v", "-o", htc_vcf_file]
    out_file_list = [htc_vcf_file + ".gz", htc_vcf_index_file]
    if use_GATK4:
        command.append("-g")
    elapsed_time = run_command(command, timeout, expected_dir, out_file_list)
    return elapsed_time


def test_mutect2(fcs_genome, timeout, expected_dir, ref, generic_fn, tumor_bam, known_snps, use_GATK4=False,
                 normal_panel_vcf=None, gnomad_vcf=None):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome mutect2 command (variant calling) on the bam file(s) and write the output.
    """
    # Run the command and track the time
    bam_dir = os.path.join(expected_dir, generic_fn + ".bam")
    mutect2_vcf_file = generic_fn + ".mutect2.vcf"
    mutect2_vcf_index_file = generic_fn + ".mutect2.vcf.gz.tbi"
    command = [fcs_genome, "mutect2", "-r", ref, "-n", bam_dir, "-t", tumor_bam, "-o", mutect2_vcf_file]
    out_file_list = [mutect2_vcf_file + ".gz", mutect2_vcf_index_file]
    if use_GATK4:
        mutect2_outdir = generic_fn + ".mutect2"
        mutect2_vcf_4_filt_file = generic_fn + ".mutect2.filt"
        out_file_list = [mutect2_vcf_4_filt_file + ".gz", mutect2_vcf_4_filt_file + ".gz.tbi",
                         mutect2_outdir + ".gz", mutect2_outdir + ".gz.tbi"]
        command = [fcs_genome, "mutect2", "-r", ref, "-n", bam_dir, "-t", tumor_bam, "-o", mutect2_outdir,
                  "-g", "--normal_name", "sample", "--tumor_name", "sample", "--panels_of_normals", normal_panel_vcf,
                  "--germline", gnomad_vcf, "--filtered_vcf", mutect2_vcf_4_filt_file]

    elapsed_time = run_command(command, timeout, expected_dir, out_file_list)
    return elapsed_time


def test_joint(fcs_genome, timeout, expected_dir, ref, gvcf_dir, generic_fn, known_snps, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome joint command (variant calling) on the bam file(s) and write the output.
    """
    # Run the command and track the tim
    joint_vcf_file = generic_fn + ".joint.vcf"
    joint_vcf_index_file = generic_fn + ".joint.vcf.gz.tbi"
    command = [fcs_genome, "joint", "-r", ref, "-i", gvcf_dir, "-o", joint_vcf_file, "--database_name", "my_database"]
    out_file_list = [joint_vcf_file + ".gz" , joint_vcf_index_file]
    if use_GATK4:
        command.append("--gatk4")
    elapsed_time = run_command(command, timeout, expected_dir, out_file_list)
    return elapsed_time


def test_ug(fcs_genome, timeout, expected_dir, ref, generic_fn, known_snps):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome ug command (variant calling) on the bam file(s) and write the output.
    """
    bam_dir = os.path.join(expected_dir, generic_fn)
    ug_vcf_file = generic_fn + ".ug.vcf"
    ug_vcf_index_file = generic_fn + ".ug.vcf.gz.tbi"
    command = [fcs_genome, "ug", "-r", ref, "-i", bam_dir, "-o", ug_vcf_file]
    out_file_list = [ug_vcf_file + ".gz", ug_vcf_index_file]
    elapsed_time = run_command(command, timeout, expected_dir, out_file_list)
    return elapsed_time


def test_germline(fcs_genome, timeout, expected_dir, ref, fastq1, fastq2, known_snps, generic_fn, use_GATK4=False):
    """
    Take in a path to the binary, path to the reference fasta (assumed index files are in the same folder), 
    a path to a bam directory (or bam file) and an output vcf file name

    Run the fcs-genome germline command (variant calling) on the bam file(s) and write the output.
    """
    bam_dir = os.path.join(expected_dir, generic_fn)
    germline_vcf_file = generic_fn + ".germline.vcf"
    command = [fcs_genome, "germline", "-r", ref, "-1", fastq1, "-2", fastq2, "-v", "-o", germline_vcf_file]
    if use_GATK4:
        command.append("--gatk4")
    out_file_list = [germline_vcf_file + ".gz"]
    elapsed_time = run_command(command, timeout, expected_dir, out_file_list)
    return elapsed_time
