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
bowtieLocation = ""
dataDir = os.path.abspath("../data/") + "/"

parser = argparse.ArgumentParser(description="Preprocess fastq files and do mapping")
parser.add_argument('--output', type=str)
parser.add_argument('input_file')
parser.add_argument('--data_dir', type=str)
parser.add_argument('genome', type=str);

results, leftovers = parser.parse_known_args()
print leftovers

print "Reads: " + results.input_file
print "Genome: " + results.genome
#print results.output
if results.data_dir is not None:
 dataDir = os.path.abspath(results.data_dir) + "/"

genomeFilename = results.genome;

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

outputDir = dataDir + inFilenamePrefixWithoutPath
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

head, tail = os.path.split(genomeFilename)
genomeIndex, file_extension = os.path.splitext(tail)
 
# check if bowtie index already exists
if (os.path.isfile(dataDir + genomeIndex + ".1.ebwt") == False):
 args = ("python",  bowtieLocation + "bowtie-build", "-o", "1", genomeFilename, dataDir+genomeIndex)
 popen = subprocess.Popen(args, stdout=subprocess.PIPE)
 popen.wait()
 output = popen.stdout.read()
 print output
 if popen.returncode != 0:
  print "error"
  sys.exit() 


bowtieInputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + "_matched_barcode" + inFileExtension
bowtieOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + ".sam"

args = ("python", bowtieLocation + "bowtie", "-S", "-p", "4", "--chunkmbs", "512", "-k", "1", "-m", "1", "-v", "2", "--strata", "--best", dataDir+genomeIndex, bowtieInputFilename, bowtieOutputFilename)
popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output
if popen.returncode != 0:
 print "error"
 sys.exit()
 
#nexus-pre
nexusOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + "_filtered.bam"

args = ("nexus-pre", bowtieOutputFilename,  "-p", "-b")
popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output
if popen.returncode != 0:
 print "error"
 sys.exit()

args = ("python", "bam_indexer.py", nexusOutputFilename)
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
