#!/usr/bin/env python

import sys
import os.path
import subprocess
import argparse
import platform

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

def indexBamFile(filename, script_path):
	args = ("python", script_path + "/bam_indexer.py", filename)
	print args
	print "Creating indexed bam file..."
	popen = subprocess.Popen(args, stdout=subprocess.PIPE)
	popen.wait()
	output = popen.stdout.read()
	if results.verbose == True:
		print output
	if popen.returncode != 0:
		print "error"
		sys.exit()
	return;

flexbar_er = "0.2";
flexbar_ol = "4";
flexbar_fm = "13";
flexbar_ml = "13";
flexbar_oh = "4";
flexbar_times = "5";
script_path = os.path.dirname(os.path.realpath(__file__))
flexbarAdapterFilename = os.path.dirname(os.path.realpath(__file__)) + "/../data/adapters_best_short.fa";
flexbarBarcodeFilename = os.path.dirname(os.path.realpath(__file__)) + "/../data/barcodes.fa";
dataDir = os.getcwd() + "/"

parser = argparse.ArgumentParser(description="Preprocess fastq files and do mapping")
parser.add_argument('--exo', action='store_true')
parser.add_argument('--clean', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--bowtie_location', nargs='?', default = "")
parser.add_argument('--num_threads', nargs='?', default = "4")
parser.add_argument('input_file')
parser.add_argument('--output_dir', type=str)
parser.add_argument('--filter_chromosomes', type=str, default="")
parser.add_argument('--random_split', action='store_true')
parser.add_argument('genome', type=str)

results, leftovers = parser.parse_known_args()

print "Reads: " + results.input_file
print "Genome: " + results.genome
#print results.output

genomeFilename = results.genome
bowtieLocation = results.bowtie_location
if len(bowtieLocation) > 0:
    bowtieLocation = bowtieLocation + "/"

if(platform.system() == "Windows" and results.bowtie_location == ""):
 print "Bowtie location is required under windows"
 sys.exit()

inputFile = os.path.abspath(results.input_file)

temp, inFileExtension = inputFile.split(os.extsep, 1)
inFileExtension, temp = os.path.splitext(inFileExtension) # remove .gz if its a fastq.gz file
inFileExtension = "." + inFileExtension
inFilenamePrefixWithoutPath = os.path.basename(inputFile)
inFilenamePrefixWithoutPath, temp = inFilenamePrefixWithoutPath.split(os.extsep, 1)

if results.output_dir is not None:
    outputDir = os.path.abspath(results.output_dir) 
    head, tail = os.path.split(results.output_dir)
    if len(tail) > 0: 
        inFilenamePrefixWithoutPath = tail;
    outputDir = outputDir + "/"
else:
    outputDir = inFilenamePrefixWithoutPath
#output has to be fastq format, because bowtie does not support fastq.gz
flexbarOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + inFileExtension

if results.exo:
 bowtieInputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + inFileExtension
else:
 bowtieInputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + "_matched_barcode" + inFileExtension
bowtieOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + ".sam"

flexbar_path = script_path+"/../bin/flexbar++"
if(platform.system() == "Windows"):
	flexbar_path += ".exe"
	
if results.exo:
 args = (flexbar_path, results.input_file, "-tt", "-t", "-ss", "-st", "-app", "-tnum", results.num_threads, "-times", flexbar_times, "-er",flexbar_er, "-ol", flexbar_ol, "-oh", flexbar_oh, "-fm", flexbar_fm, "-ml", flexbar_ml, "-a", flexbarAdapterFilename,"-o", flexbarOutputFilename)
else:
 args = (flexbar_path, results.input_file, "-tl", "5", "-tt", "-t", "-ss", "-st", "-app", "-tnum", results.num_threads, "-times", flexbar_times, "-er",flexbar_er, "-ol", flexbar_ol, "-oh", flexbar_oh, "-fm", flexbar_fm, "-ml", flexbar_ml, "-b", flexbarBarcodeFilename, "-a", flexbarAdapterFilename,"-o", flexbarOutputFilename)
 #args = ("flexbar++", results.input_file, "-tl", "5","-b", flexbarBarcodeFilename, "-a", flexbarAdapterFilename,"-o", flexbarOutputFilename)
if not os.path.exists(outputDir):
 os.makedirs(outputDir)
if (os.path.isfile(bowtieInputFilename) == False or results.overwrite == True):
    print "Filtering pre-mapping barcodes and trimming adapters..."
    if results.verbose == True:
        popen = subprocess.Popen(args + tuple(leftovers))
    else:
        popen = subprocess.Popen(args + tuple(leftovers), stdout=subprocess.PIPE)
    popen.wait()
    if popen.returncode != 0:
     print "error"
     sys.exit()
    if results.verbose == True:
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
 if results.verbose == True:
    print output
 if popen.returncode != 0:
  print "error"
  sys.exit() 


if (os.path.isfile(bowtieOutputFilename) == False or results.overwrite == True):
    args = (bowtieLocation + "bowtie", "-S", "-p", results.num_threads, "--chunkmbs", "512", "-k", "1", "-m", "1", "-v", "2", "--strata", "--best", genomeIndexFile, bowtieInputFilename, bowtieOutputFilename)
    if(platform.system() != "Linux" and platform.system() != "Linux2"):
     args = ("python",) + args
    print "Mapping reads..."
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    if results.verbose == True:
        print output
    if popen.returncode != 0:
     print "error"
     sys.exit()
 
#nexus-pre
nexusOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + "_filtered.bam"
if results.random_split == True:
	nexusOutputFilenameSplit1 = outputDir + "/" + inFilenamePrefixWithoutPath + "_filtered_split1.bam"
	nexusOutputFilenameSplit2 = outputDir + "/" + inFilenamePrefixWithoutPath + "_filtered_split2.bam"

nexus_pre_path = script_path+"/../bin/nexus-pre"
if(platform.system() == "Windows"):
	nexus_pre_path += ".exe"
	
	
if (os.path.isfile(nexusOutputFilename) == False or 
(results.random_split == True and (os.path.isfile(nexusOutputFilenameSplit1) == False or os.path.isfile(nexusOutputFilenameSplit2) == False)) or
results.overwrite == True):
    args = (nexus_pre_path, bowtieOutputFilename,  "-fc", results.filter_chromosomes)
    if results.random_split == True:
		args += ("-rs",)
    print "Filtering post-mapping barcodes..."
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    if results.verbose == True:
        print output
    if popen.returncode != 0:
        print "error"
        sys.exit()

indexBamFile(nexusOutputFilename, script_path)			
if results.random_split == True:
	indexBamFile(nexusOutputFilenameSplit2, script_path)		
	indexBamFile(nexusOutputFilenameSplit1, script_path)		
	
#cleanup
if results.clean:
    print "deleting intermediate files..."
    os.remove(bowtieOutputFilename)
    os.remove(nexusOutputFilename)
    if results.exo:
        os.remove(flexbarOutputFilename)
    else:
        os.remove(outputDir + "/" + inFilenamePrefixWithoutPath + "_matched_barcode" + inFileExtension)
        os.remove(outputDir + "/" + inFilenamePrefixWithoutPath + "_unidentified" + inFileExtension)
