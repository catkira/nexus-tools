/*
Author: Benjamin Menkuec
Copyright 2015 Benjamin Menkuec
License: LGPL
*/
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <random>

#include "readsim.h"

seqan::ArgumentParser buildParser(void)
{
    seqan::ArgumentParser parser;

    setCategory(parser, "Chip Nexus/Exo Read Simulator");
    setShortDescription(parser, "Chip Nexus/Exo Read Simulator");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "");

    addDescription(parser, "");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseOption genomeOpt = seqan::ArgParseOption(
        "g", "genome", "FastA file containing the reference genome.",
        seqan::ArgParseArgument::INPUT_FILE, "GENOME");
    setValidValues(genomeOpt, seqan::SeqFileIn::getFileExtensions());
    addOption(parser, genomeOpt);

    seqan::ArgParseOption adapterOpt = seqan::ArgParseOption(
        "a", "adapters", "FastA file containing adapters.",
        seqan::ArgParseArgument::INPUT_FILE, "ADAPTERS");
    setValidValues(adapterOpt, seqan::SeqFileIn::getFileExtensions());
    addOption(parser, adapterOpt);

    seqan::ArgParseOption numReadsOpt = seqan::ArgParseOption(
        "n", "num_peaks", "Number of peaks to generate.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setMinValue(numReadsOpt, "0");
    addOption(parser, numReadsOpt);

    seqan::ArgParseOption fixedBarcodeOpt = seqan::ArgParseOption(
        "fb", "barcode", "fixed barcode for Chip Nexus reads",
        seqan::ArgParseOption::STRING, "SEQUENCE");
    setDefaultValue(fixedBarcodeOpt, "CTAG");
    addOption(parser, fixedBarcodeOpt);

    seqan::ArgParseOption randomBarcodeOpt = seqan::ArgParseOption(
        "rb", "random_barcode", "number of bases for random barcode",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(randomBarcodeOpt, 5);
    addOption(parser, randomBarcodeOpt);

    seqan::ArgParseOption readLengthOpt = seqan::ArgParseOption(
        "rl", "read length", "number of bases per read",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(readLengthOpt, 44);
    addOption(parser, readLengthOpt);

    return parser;
}

void generateRandomBarcode(std::string& randomBarcode, const int n)
{
    randomBarcode.clear();
    if (n <= 0)
        return;
    for (unsigned int k = 0; k < n;k++)
    {
        switch (rand() % 4) {
        case 0:
            randomBarcode.append("A"); break;
        case 1:
            randomBarcode.append("C"); break;
        case 2:
            randomBarcode.append("G"); break;
        case 3:
            randomBarcode.append("T"); break;
        }
    }
}

void substituteAdapter(std::string& read, const std::vector<std::string> adapters)
{
    unsigned int minAdapterLength = 4;
    unsigned int maxAdapterLength = 20;

    std::default_random_engine generator;
    std::uniform_int_distribution<unsigned int> distribution(minAdapterLength, maxAdapterLength);
    unsigned int adapter = rand() % adapters.size();
    unsigned int adapterLength = distribution(generator);
    if (adapterLength > adapters[adapter].size())
        adapterLength = adapters[adapter].size();
    read.replace(read.size() - adapterLength, adapterLength, adapters[adapter].substr(0,adapterLength));
}

int main(int argc, char const ** argv)
{
    seqan::ArgumentParser parser = buildParser();
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if input was successfully parsed.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::string refGenome;
    if (isSet(parser, "g"))
    {
        std::string genomeFilePath;
        getOptionValue(genomeFilePath, parser, "g");

        seqan::SeqFileIn genomeFile;
        if (!open(genomeFile, genomeFilePath.c_str(), seqan::OPEN_RDONLY))
        {
            std::cerr << "Error while opening file'" << genomeFilePath << "'.\n";
            return 1;
        }
        std::cout << "reading reference genome... " << std::endl;
        while (!atEnd(genomeFile))
        {
            std::string id;
            std::string bases;
            readRecord(id, bases, genomeFile);
            refGenome.append(bases);
        }
        std::cout <<"    read "<<refGenome.size() << " bases" << std::endl;
    }
    std::vector<std::string> adapters;
    if (isSet(parser, "a"))
    {
        std::string adapterFilePath;
        getOptionValue(adapterFilePath, parser, "a");

        seqan::SeqFileIn adapterFile;
        if (!open(adapterFile, adapterFilePath.c_str(), seqan::OPEN_RDONLY))
        {
            std::cerr << "Error while opening file'" << adapterFilePath << "'.\n";
            return 1;
        }
        std::cout << "reading adapters... " << std::endl;
        while (!atEnd(adapterFile))
        {
            std::string id;
            std::string adapter;
            readRecord(id, adapter, adapterFile);
            adapters.emplace_back(adapter);
        }
        std::cout << "    read " << adapters.size() << " adapters" << std::endl;
    }

    int numReads = 0;
    if (isSet(parser, "n"))
        getOptionValue(numReads, parser, "n");

    std::string fixedBarcode;
    getOptionValue(fixedBarcode, parser, "fb");

    int numRandomBarcode = -1;
    getOptionValue(numRandomBarcode, parser, "rb");

    int readLength = 0;
    getOptionValue(readLength, parser, "rl");

    const unsigned int peakHalfWidth = 10;
    const unsigned int peakStdDev = 10;
    const unsigned int peakNumReads = 50;

    std::default_random_engine generator;
    std::normal_distribution<float> distribution((float)peakHalfWidth, (float)peakStdDev);

    unsigned int numPCRArtifacts = 0;
    unsigned int numAdapters = 0;

    for (unsigned int n = 0;n < numReads;++n)
    {
        unsigned int peakPos = rand() % (refGenome.size() - readLength - peakHalfWidth*2 + numRandomBarcode + fixedBarcode.size());
        unsigned int k = 0;
        while(k<peakNumReads)
        {
            unsigned int pos = (unsigned int)distribution(generator);
            if (pos < peakPos - peakHalfWidth || peakPos > peakPos + peakHalfWidth)
                continue;
            pos += peakPos;
            // add fixed and random barcode
            std::string randomBarcode;
            generateRandomBarcode(randomBarcode, numRandomBarcode);
            std::string read = randomBarcode + fixedBarcode + refGenome.substr(pos, readLength - numRandomBarcode - fixedBarcode.size());
            // add adapter
            substituteAdapter(read, adapters);
            ++numAdapters;
            // add PCR artifacts
            unsigned int PCRArtifactPercentage = 10;
            if ((rand() % 100) < PCRArtifactPercentage)
            {
                std::cout << pos << ": " << read << std::endl;
                ++numPCRArtifacts;
            }

            std::cout << pos << ": " << read << std::endl;
            ++k;
        }
    }


}
               