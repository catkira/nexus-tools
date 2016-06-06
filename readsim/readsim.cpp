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
        "n", "num_reads", "Number of reads to generate.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setMinValue(numReadsOpt, "0");
    addOption(parser, numReadsOpt);

    seqan::ArgParseOption fixedBarcodeOpt = seqan::ArgParseOption(
        "fb", "barcode", "fixed barcode for Chip Nexus reads",
        seqan::ArgParseOption::STRING, "SEQUENCE");
    addOption(parser, fixedBarcodeOpt);

    seqan::ArgParseOption randomBarcodeOpt = seqan::ArgParseOption(
        "rb", "random_barcode", "number of bases for random barcode",
        seqan::ArgParseOption::INTEGER, "VALUE");
    addOption(parser, randomBarcodeOpt);

    seqan::ArgParseOption readLengthOpt = seqan::ArgParseOption(
        "rl", "read length", "number of bases per read",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(readLengthOpt, 44);
    addOption(parser, readLengthOpt);

    return parser;
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
    if (isSet(parser, "fb"))
        getOptionValue(fixedBarcode, parser, "fb");

    int numRandomBarcode = -1;
    if (isSet(parser, "rb"))
        getOptionValue(numRandomBarcode, parser, "rb");

    int readLength = 0;
    getOptionValue(readLength, parser, "rl");


    for (unsigned int n = 0;n < numReads;++n)
    {

    }


}
