#!/usr/bin/env python2.7
# Chris Eisenhart chrisEisenhart1992@gmail.com
"""
"""
import sys
import operator
import argparse
import subprocess
import os

def parseArgs(args): 
    """ 
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using args.option. For example args.a stores the alphabet provided. 
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument ("input_mode",
            help = "Run the program",
            action = "store")
    parser.add_argument("--verbose",
            # Sets warning tags to on
            help = " The program will send warning messages to stderr. This is very messy, use with caution.",
            action= "store_false")

    parser.set_defaults(input_mode = True)
    args = parser.parse_args() # Initializes the options specified above
    return args # Returns the parser as an object


def which(name):
    found = 0 
    result = ""
    for path in os.getenv("PATH").split(os.path.pathsep):
        full_path = path + os.sep + name
        if os.path.exists(full_path):
            """
            if os.stat(full_path).st_mode & stat.S_IXUSR:
                found = 1
                print(full_path)
            """
            found = 1
            result = full_path
    # Return a UNIX-style exit code so it can be checked by calling scripts.
    # Programming shortcut to toggle the value of found: 1 => 0, 0 => 1.
    return result
    sys.exit(1 - found)


def main(args):
    """
    """
    options = parseArgs(args)
    sw_tb = os.path.isfile(os.environ['SW_TB'])
    pmm_tb = os.path.isfile(os.environ['PMM_TB'])
    if os.environ['CLOUD'] != "aws":
        smem_bit = os.path.isfile(os.environ['SMEM_TB'])

    sw_bit = os.path.isfile(os.environ['SW_BIT'])
    pmm_bit = os.path.isfile(os.environ['PMM_BIT'])
    if os.environ['CLOUD'] != "aws":
        smem_bit = os.path.isfile(os.environ['SMEM_BIT'])

    test = [os.environ['SW_TB'], os.environ['SW_BIT'], os.environ['ref_genome'],
            os.path.join(os.environ['tbdata_dir'], "sw", "input"),
            os.path.join(os.environ['tbdata_dir'], "sw", "golden_out")]

    proc = subprocess.Popen(test,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    output, error = proc.communicate()
    error_code = proc.returncode

    test2 = [os.environ['PMM_TB'], os.environ['PMM_BIT'], os.path.join(os.environ['tbdata_dir'], "pmm")]
    proc2 = subprocess.Popen(test2,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    output2, error2 = proc2.communicate()
    error_code2 = proc2.returncode

    xbutil_path = which("xbutil")
    if xbutil_path:
        if xbutil_path != "/dev/null":
            test3 = [xbutil_path, "dmatest"]
            proc3 = subprocess.Popen(test3, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            output3, error3 = proc3.communicate()
            error_code3 = proc3.returncode

###########################################################################################################
    # All tests have been attempted and results are stored in python variables, go over the variables and
    # report any that have failed. 


    if not sw_tb:
        print("Could not find the SW_TB file")
    if not pmm_tb:
        print("Could not find the PMM_TB file")
    if not sw_bit:
        print("Could not find the SW_BIT file")
    if not pmm_bit:
        print("Could not find the PMM_BIT file")

    if error_code != 0:
        print("Testing the SW_TB and SW_BIT files with command '{}' failed!".format(" ".join(test)))
        print("Stdout = {}".format(output))
        print("Stderr = {}".format(error))

    if error_code2 != 0:
        print("Testing the PMM_TB and PMM_BIT files with command '{}' failed!".format(" ".join(test2)))
        print("Stdout = {}".format(output2))
        print("Stderr = {}".format(error2))

    if xbutil_path:
        if xbutil_path != "/dev/null":
            if error_code3 !=0:
                print("Testing the xbutil file with command '{}' failed!".format(" ".join(test3)))
                print("Stdout = {}".format(output3))
                print("Stderr = {}".format(error3))


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
