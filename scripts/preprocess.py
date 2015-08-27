#!/usr/bin/env python

import sys
import os.path
import subprocess
import argparse
import platform

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

flexbar_er = "0.2";
flexbar_ol = "4";
flexbar_fm = "22";
flexbar_ml = "18";
flexbar_tnum = "4";
flexbarAdapterFilename = os.path.dirname(os.path.realpath(__file__)) + "/../data/adapters.fa";
flexbarBarcodeFilename = os.path.dirname(os.path.realpath(__file__)) + "/../data/barcodes.fa";
dataDir = os.getcwd() + "/"

parser = argparse.ArgumentParser(description="Preprocess fastq files and do mapping")
parser.add_argument('--output', type=str)
parser.add_argument('--exo', action='store_true')
parser.add_argument('--bowtie_location', nargs='?', default = "")
parser.add_argument('input_file')
parser.add_argument('--data_dir', type=str)
parser.add_argument('genome', type=str)

results, leftovers = parser.parse_known_args()
print leftovers

print "Reads: " + results.input_file
print "Genome: " + results.genome
#print results.output
if results.data_dir is not None:
 dataDir = os.path.abspath(results.data_dir) + "/"

genomeFilename = results.genome
bowtieLocation = results.bowtie_location

if(platform.system() == "Windows" and results.bowtie_location == ""):
 print "Bowtie location is required under windows"
 sys.exit()

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
inFilenamePrefixWithoutPath, temp = os.path.splitext(inFilenamePrefixWithoutPath)

outputDir = dataDir + inFilenamePrefixWithoutPath
flexbarOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + inFileExtension

if results.exo:
 print "processing Chip-Exo data"
 args = ("seqan_flexbar", results.input_file, "-tt", "-t", "-ss", "-tnum", flexbar_tnum,"-er",flexbar_er, "-ol", flexbar_ol, "-fm", flexbar_fm, "-ml", flexbar_ml, "-a", flexbarAdapterFilename,"-o", flexbarOutputFilename)
else:
 print "processing Chip-Nexus data"
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
genomeIndexFile = os.path.dirname(genomeFilename) + "/" + genomeIndex;
if (os.path.isfile(genomeIndexFile + ".1.ebwt") == False):
 if(platform.system() == "Linux" or platform.system() == "Linux2"):
  args = (bowtieLocation + "bowtie-build", "-o", "1", genomeFilename, genomeIndexFile)
 else:
  args = ("python",  bowtieLocation + "bowtie-build", "-o", "1", genomeFilename, genomeIndexFile)
 popen = subprocess.Popen(args, stdout=subprocess.PIPE)
 popen.wait()
 output = popen.stdout.read()
 print output
 if popen.returncode != 0:
  print "error"
  sys.exit() 

if results.exo:
 bowtieInputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + inFileExtension
else:
 bowtieInputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + "_matched_barcode" + inFileExtension
bowtieOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + ".sam"

if(platform.system() == "Linux" or platform.system() == "Linux2"):
 args = (bowtieLocation + "bowtie", "-S", "-p", "4", "--chunkmbs", "512", "-k", "1", "-m", "1", "-v", "2", "--strata", "--best", genomeIndexFile, bowtieInputFilename, bowtieOutputFilename)
else:
 args = ("python", bowtieLocation + "bowtie", "-S", "-p", "4", "--chunkmbs", "512", "-k", "1", "-m", "1", "-v", "2", "--strata", "--best", genomeIndexFile, bowtieInputFilename, bowtieOutputFilename)
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

if(platform.system() == "Linux" or platform.system() == "Linux2"):
 args = (os.path.dirname(os.path.realpath(__file__)) + "/bam_indexer.py", nexusOutputFilename)
else:
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
os.remove(nexusOutputFilename)
