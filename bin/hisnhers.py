#!/usr/bin/env python

'''
hisnhers.py
Harvard Informatics Script for Nextgen HiSeq Extraction and Reporting of Sequences

 Starts with wrong chmod and bad shebang
 The __main__ bit at the end is missing.

 Has a bad indent on a function

 Read the fastq files into a list of sequences
 The function for converting them into a list of sequences exists, but it requires 
 a file handle.  Script has file path hard coded to a file that doesn't exist. 
 Empty lines are returned as list element. Whitespace is not stripped.

 Exercise: Write a file-exists check and exit if it fails; use Google to find the answer
 Exercise: Use a 'with' block to get a file handle
 Exercise: Iterate over the sequences and report the length.  Add the sequence iteration key
 Exercise: Iterate over the sequences and report the A, T, C, and G counts.  Use string formatting. Make sure you've done "lower()"
 Exercise: Report length and counts after skipping the adapter in the first 20 bases
 Exercise: Get the file name as a command line argument
 Exercise: Fix the conversion function to create a tuple with the sequence id, bases, and quality scores
 Exercise: Fix the conversion function to get rid of empty lines and strip whitespace.


 Shells out to "megaAssembler", redirecting stderr to /dev/null,  
 and attempts to read the output file which fails.
 Turns out it is reporting "out of memory" on stderr and is missing a required argument.
 Useful information is written on standard out in addition to writing an assembled contig to a file.

 Exercise: Wrap file open in an exists check and a context manager block.
 Exercise: Use args array and "join" to make the command line
 Exercise: Write a popen.subprocess.  Report the return code, stdout, and stderr.  Exit with stderr message if
 return code is not 0
 Exercise: Replace megaAssembler with "hyperAssembler"
 Exercise: Parse the start time and end time reported by hyperAssembler into a python datetime and get the delta

 Interlude for PyPi, packages, virtual environments
 Exercise: Google install python dateutil.  Try it.
 Exercise: Create an Anaconda clone.  Try to install python-dateutil with pip; --upgrade.
 

 Reads hyperAssembler output, some space / tab-delimited lines, sorts them and prints them out.

 Exercise: split lines into elements with the "split" function and create a dictionary keyed by ??
 Exercise: created "sorted" dict using sorted and passing a sort function / lambda
 Exercise: use OrderedDict because you can't guarantee key sort order in a regular dict
 Exercise: convert split, sort, and dict creation into a dictionary comprehension with an OrderedDict?


 Import the annotator module and annotate each contig. 
 Exercise (repl): Use dir() on the annotation module to find out what functions are there.  Distinguish between
 those with underbars and those without
 Exercise (repl): Use help() to get the docstrings

 Convert annotation to parallel processing
 Exercise: Use multiprocessing to do a pool of parallel annotations

 Convert annotation to MPI
 Exercise: conda install mpi4py, test (it should fail); conda remove mpi4py
 Exercise: load gcc / mpi module, do pip install of mpi4py;

 Exercise: Use mpi4py to make MPI code.  Use a try block to catch import errors and default to multiprocessing.

'''

import os, traceback, sys, re
import json
import subprocess
from argparse import ArgumentParser, RawDescriptionHelpFormatter


# Logging setup
import logging

DEBUG = os.environ.get("HISNHERS_DEBUG",0)


def annotate(seqid, bases):
    '''
    Annotate the contig
    '''
    pass


def fastqToSequenceList(fileh):
    '''
    Takes a fastq file handle, returns a list tuples including
    seqid, sequence bases, and quality scores
    '''
    seqs = []
    seqid = bases = qscores = None


    for line in fileh:
        if line.strip() == '':
            continue

        if line.startswith('@'):
            # It's a new record
            seqs.append((seqid,bases,qscores))
            seqid = bases = qscores = None

            m = re.match(r'^@([^ ]+).*',line)
            if m is None:
                raise Exception('No sequence identifier found after the @ symbol')
            seqid = m.groups(1)

        elif seqid is not None and bases is None and line.strip() != '+':
            bases = line.strip()
        elif seqid is not None and bases is not None and line.strip() != '+':
            qscores = line.strip()

    return seqs


def main(argv = None):

    # Setup argument parser
    parser = ArgumentParser(description='Some code', formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
    parser.add_argument('FASTQ_FILE',help='Fastq file')
    parser.add_argument('--skip-adapter', dest='adapter_sequence', help='Skip 20 bases of adapter sequence.')
    args = parser.parse_args()

    # Read fastq file and report length, base counts
    seqs = []
    fqfilename = args.FASTQ_FILE
    with open(fqfilename,'r') as f:
        seqs = fastqToSequenceList(f)

    adapter_sequence = args.adapter_sequence
    for i,seqdata in enumerate(seqs):
        seqstr = seqdata[1]
        if adapter_sequence:
            seqstr = seqstr[20:]

        print 'Length %d: %d' % (i,len(seqstr))
        basecounts = 'Base counts- '
        for base in ['A','T','C','G']:
            basecounts += '%s: %d' % (base,seqstr.count(base))
        print basecounts


    # Run assembler with fastq file input and read the output contig
    contigfilename = '%s.contigs' % fqfilename
    assemblerargs = [
        'megaAssembler',
        '-i',
        fqfilename,
        '-o',
        contigfilename,
    ]
    cmd = ' '.join(assemblerargs)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    returncode, stdoutstr, stderrstr = proc.communicate()

    if returncode != 0:
        raise Exception('Error running assembler with cmd %s\nstdout: %s\nstderr: %s' % (cmd,stdoutstr,stderrstr))

    contigs = []
    with open(contigfilename,'r') as c:
        seqid = None
        for line in c:
            if line.strip() == '':
                continue
            # Check for seq id line
            m = re.match(r'^>\s*([^ ]+).*', line)
            if m is not None:
                seqid = m.groups(1)
            else:
                # Otherwise, it's the bases
                contigs.append((seqid, line.strip()))
                seqid = None

    # Using a multiprocessing Pool
    from multiprocessing import Pool
    numprocs = os.environ.get('ANNOTATION_PROC_NUM',4)
    pool = Pool(numprocs)
    annotations = pool.map(annotate,contigs)


    # Using MPI
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    # size = comm.Get_size()
    # rank = comm.Get_rank()
    
    contigs = comm.scatter(contigs,root=0)
    annots = annotate(contigs)

    annotations = comm.gather(annots,root=0)


    # # One at a time
    # for seqid, bases in contigs:
    #     annotations += annotate(seqid, bases)


    # Dump annotations in JSON form
    with open('%s.annotations' % fqfilename, 'w') as f:
        f.write(json.dumps(annotations))


if __name__ == '__main__':
    sys.exit(main())
