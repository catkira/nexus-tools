#!/usr/bin/env python

import sys
import os.path
import subprocess
import argparse

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

flexbar_er = "0.2";
flexbar_ol = "4";
flexbar_fm = "22";
flexbar_ml = "18";
flexbar_tnum = "4";
flexbarAdapterFilename = "../data/adapters.fa";
flexbarBarcodeFilename = "../data/barcodes.fa";
bowtieIndexFilename = "P:/bowtie-1.1.2/dm3"
bowtieLocation = "P:/bowtie-1.1.2/bowtie"

parser = argparse.ArgumentParser(description="Preprocess fastq files and do mapping")
parser.add_argument('--output', type=str)
parser.add_argument('input_file')

results, leftovers = parser.parse_known_args()
print leftovers

print results.input_file
#print results.output

if results.output is not None:
 outFilenamePrefix, file_extension = os.path.splitext(results.output)
 outFilenamePrefix = os.path.dirname(results.input_file)+ "/" +outFilenamePrefix;
else:
 outFilenamePrefix, file_extension = os.path.splitext(results.input_file)
 outFilenamePrefix += ""
#print "output file: " + outFilenamePrefix + ".bam"

inputFile = results.input_file

inFilenamePrefix, inFileExtension = os.path.splitext(results.input_file)
inFilenamePrefixWithoutPath = os.path.basename(results.input_file);
inFilenamePrefixWithoutPath, temp = os.path.splitext(inFilenamePrefixWithoutPath)

outputDir = os.path.abspath("../data/" + inFilenamePrefixWithoutPath) 
flexbarOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + inFileExtension

args = ("seqan_flexbar", results.input_file, "-tl", "5", "-tt", "-t", "-ss", "-tnum", flexbar_tnum,"-er",flexbar_er, "-ol", flexbar_ol, "-fm", flexbar_fm, "-ml", flexbar_ml, "-b", flexbarBarcodeFilename, "-a", flexbarAdapterFilename,"-o", flexbarOutputFilename)
#args = ("seqan_flexbar", results.input_file, "-tl", "5","-b", flexbarBarcodeFilename, "-a", flexbarAdapterFilename,"-o", flexbarOutputFilename)
if not os.path.exists(outputDir):
 os.makedirs(outputDir)
print args
popen = subprocess.Popen(args + tuple(leftovers))
popen.wait()
if popen.returncode != 0:
 print "error"
 sys.exit()
print flexbarOutputFilename + " created"


bowtieInputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + "_matched_barcode" + inFileExtension
bowtieOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + ".sam"

args = ("python", bowtieLocation, "-S", "-p", "4", "--chunkmbs", "512", "-k", "1", "-m", "1", "-v", "2", "--strata", "--best", bowtieIndexFilename, bowtieInputFilename, bowtieOutputFilename)
popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output
if popen.returncode != 0:
 print "error"
 sys.exit()

args = ("python", "bam_indexer.py", bowtieOutputFilename)
popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output
if popen.returncode != 0:
 print "error"
 sys.exit()
 
 
#cleanup
os.remove(bowtieOutputFilename)
os.remove(outputDir + "/" + inFilenamePrefixWithoutPath + ".bam")
