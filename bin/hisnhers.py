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
import time


# Logging setup
import logging

DEBUG = os.environ.get("HISNHERS_DEBUG",0)


def runcmd(cmd):
    '''
    Execute a command and return stdout, stderr, and the return code
    '''
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdoutstr, stderrstr = proc.communicate()
    return (proc.returncode,stdoutstr,stderrstr)



def fastqToSequenceList(fileh):
    '''
    Takes a fastq file handle, returns a list tuples including
    seqid, sequence bases, and quality scores
    '''
    seqs = []
    seqid = bases = qscores = None

    if fileh.closed:
        raise Exception('Fastq file is closed.')

    for line in fileh:
        if line.startswith('@'):
            # It's a new record
            if seqid is not None:
                seqs.append((seqid,bases,qscores))
            seqid = bases = qscores = None

            m = re.match(r'^@([^ ]+).*',line)
            if m is None:
                raise Exception('No sequence identifier found after the @ symbol')
            seqid = m.group(1)

        elif seqid is not None and bases is None and line.strip() != '+':
            bases = line.strip()
        elif seqid is not None and bases is not None and line.strip() != '+':
            qscores = line.strip()

    return seqs


def main():

    # Read fastq file and report length, base counts
    seqs = []

    if len(sys.argv) < 2:
        print 'Must supply a file name'
        return 1
    fqfilename = sys.argv[1]

    with open(fqfilename,'r') as f:
        seqs = fastqToSequenceList(f)

    adapter_sequence = None
    for i,seqdata in enumerate(seqs):
        seqstr = seqdata[1]
        if adapter_sequence:
            seqstr = seqstr[20:]

        print 'Length %d: %d' % (i,len(seqstr))
        basecounts = 'Base counts- '
        for base in ['A','T','C','G']:
            basecounts += '%s: %d\t' % (base,seqstr.count(base))
        print basecounts


    # Write out sequences in fasta format
    (path,ext) = os.path.splitext(fqfilename)
    fafilename = path + '.fa'
    print 'Writing to %s' % fafilename
    with open(fafilename,'w') as f:
        for seqdata in seqs:
            f.write('>%s\n%s\n' % (seqdata[0],seqdata[1]))


    # # Run megaAssembler with fastq file input and read the output contig
    # contigfilename = '%s.contigs' % fqfilename
    # assemblerargs = [
    #     'megaAssembler',
    #     fqfilename,
    # ]

    # Run hyperAssembler with fastq file input and read the output contig
    contigfilename = '%s.contigs' % fafilename
    assemblerargs = [
        'hyperAssembler',
        fafilename,
    ]

    cmd = ' '.join(assemblerargs)
    returncode, stdoutstr, stderrstr = runcmd(cmd)

    if returncode != 0:
        raise Exception('Error running assembler with cmd %s\nstdout: %s\nstderr: %s' % (cmd,stdoutstr,stderrstr))


    # Get the start and end time from stdout
    from dateutil import parser
    match = re.search(r'Start time: (.*)\n', stdoutstr, re.MULTILINE)
    if match:
        starttime = parser.parse(match.group(1))
    match = re.search(r'End time: (.*)\n', stdoutstr, re.MULTILINE)
    if match:
        endtime = parser.parse(match.group(1))
    if starttime and endtime:
        delta = endtime - starttime
        print 'Elapsed assembly time %d seconds' % delta.total_seconds()


    contigs = []
    with open(contigfilename,'r') as c:
        seqid = None
        for line in c:
            line = line.strip()
            if line == '':
                continue
            # Check for seq id line
            m = re.match(r'^>\s*([^ ]+).*', line)
            if m is not None:
                seqid = m.group(1)
            else:
                # Otherwise, it's the bases
                contigs.append((seqid, line))
                seqid = None


    from ha.annotate import annotateStartStopCodons, annotatePalindromes

    # Using a multiprocessing Pool
    starttime = time.time()
    from multiprocessing import Pool
    numprocs = os.environ.get('ANNOTATION_PROC_NUM',2)
    pool = Pool(numprocs)

    annotations = []
    results = []
    for contig in contigs:
        result = pool.apply_async(annotateStartStopCodons,contig)
        results.append(result)
        result = pool.apply_async(annotatePalindromes,contig)
        results.append(result)

    for result in results:
        annotations += result.get()

    endtime = time.time()
    print 'Elapsed parallel annotation time %d seconds' % int(endtime - starttime)


    # # Using MPI
    # from mpi4py import MPI
    # comm = MPI.COMM_WORLD
    # # size = comm.Get_size()
    # # rank = comm.Get_rank()
    
    # contigs = comm.scatter(contigs,root=0)
    # annots = annotate(contigs)

    # annotations = comm.gather(annots,root=0)


    # One at a time
    annotations = []
    starttime = time.time()
    for seqid, contig in contigs:
        annotations += annotateStartStopCodons(seqid, contig)
        annotations += annotatePalindromes(seqid, contig)
    endtime = time.time()
    print 'Elapsed serial annotation time %d seconds' % int(endtime - starttime)

    # Make a dictionary keyed by contig name
    annotatedcontigs = {}
    for annotation in annotations:
        annotatedcontigs.setdefault(annotation['seqid'],[]).append(annotation)

    # Sort the annotations by start location
    for contig,annotations in annotatedcontigs.iteritems():
        annotatedcontigs[contig].sort(key=lambda annot: annot['start'])

    # Dump annotations in JSON form
    with open('%s.annotations' % fafilename, 'w') as f:
        f.write(json.dumps(annotatedcontigs,indent=4))
if __name__ == '__main__':
    sys.exit(main())