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
#include <boost/algorithm/string.hpp>

#include "readsim.h"

std::string getFilePrefix(const std::string& fileName, const bool withPath = true)
{
    std::size_t found = fileName.find_last_of(".");
    std::size_t found2 = std::string::npos;
    if (!withPath)
        found2 = fileName.find_last_of("/\\");
    if (found == std::string::npos)
        return std::string();
    if (found2 == std::string::npos)
        return fileName.substr(0, found);
    return fileName.substr(found2 + 1, found - found2);
}

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

    seqan::ArgParseOption erOpt = seqan::ArgParseOption(
        "er", "error-rate", "Error rate in generated reads",
        seqan::ArgParseOption::DOUBLE, "VALUE");
    setDefaultValue(erOpt, 0.2);
    addOption(parser, erOpt);


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

    seqan::ArgParseOption outPrefixOpt = seqan::ArgParseOption(
        "o", "output_filename", "Base filename of output file without suffix",
        seqan::ArgParseOption::STRING, "STRING");
    setDefaultValue(outPrefixOpt, "readsim_out");
    addOption(parser, outPrefixOpt);

    seqan::ArgParseOption refFileOpt = seqan::ArgParseOption(
        "fr", "reference_file", "File containing the ideal post preprocessing file",
        seqan::ArgParseOption::STRING, "STRING");
    addOption(parser, refFileOpt);

    seqan::ArgParseOption preFileOpt = seqan::ArgParseOption(
        "fp", "preprocessing_file", "File from preprocessing",
        seqan::ArgParseOption::STRING, "STRING");
    addOption(parser, preFileOpt);

    seqan::ArgParseOption qualityOpt = seqan::ArgParseOption(
        "q", "quality", "Simulate quality detoriation");
    addOption(parser, qualityOpt);

    return parser;
}

char getRandomBase()
{
    switch (rand() % 4) {
    case 0:
        return 'A'; break;
    case 1:
        return 'C'; break;
    case 2:
        return 'G'; break;
    case 3:
        return 'T'; break;
    }
    return 'A';
}

void generateRandomBarcode(std::string& randomBarcode, const unsigned int n)
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

unsigned int substituteAdapter(std::string& read, const std::vector<std::string> adapters)
{
    unsigned int minAdapterLength = 0;
    unsigned int maxAdapterLength = 20;

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<unsigned int> distribution(minAdapterLength, maxAdapterLength);
    unsigned int adapter = rand() % adapters.size();
    unsigned int adapterLength = distribution(generator);
    if (rand() % 2) // only add adapters in 50% of all cases
        adapterLength = 0;

    if (adapterLength > adapters[adapter].size())
        adapterLength = adapters[adapter].size();
    read.replace(read.size() - adapterLength, adapterLength, adapters[adapter].substr(0,adapterLength));
    return adapterLength;
}

unsigned int insertErrors(std::string& read, const double er)
{
    unsigned int nErrors = 0;
    auto it = read.begin();
    while (it != read.end())
    {
        if ((unsigned int)(rand() % 100) < (unsigned int)(er * 100))
        {
            *it = getRandomBase();
            ++nErrors;
        }
        ++it;
    }
    return nErrors;
}

std::vector<int> doQualities(seqan::Dna5QString& read)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    const unsigned int maxStartPos = length(read);
    const unsigned int minStartPos = maxStartPos / 2;
    std::uniform_int_distribution<unsigned int> distribution(minStartPos, maxStartPos);
    unsigned int startPos = distribution(generator);
    const int bestQuality = 40;
    unsigned char q = bestQuality;
    std::vector<int> qualities;
    qualities.assign(maxStartPos, bestQuality);
    for (auto pos = startPos; pos < maxStartPos; pos++)
    {
        if (q >= 5)
            q -= 5;
        else
            q = 0;
        const float percentCorrect = ((float)(q)/ bestQuality) * 100;
        if (rand() % 100 > percentCorrect)
        {
            read[pos] = getRandomBase();
        }
        qualities[pos] = q;
        seqan::assignQualityValue(read[pos], (int)q);
    }
    return qualities;
}

void assignQualities(seqan::Dna5QString& read, const std::vector<int>& qualities)
{
    const unsigned int length = seqan::length(read);
    for (unsigned int i = 0; i < length; i++)
    {
        seqan::assignQualityValue(read[i], qualities[i]);
    }
}

float Sen(unsigned int TP, unsigned int FN)
{
    return static_cast<float>(TP) / static_cast<float>(TP + FN);
}

float Spec(unsigned int TN, unsigned int FP)
{
    return static_cast<float>(TN) / static_cast<float>(TN + FP);
}

float PPV(unsigned int TP, unsigned int FP)
{
    return static_cast<float>(TP) / static_cast<float>(TP + FP);
}

float NPV(unsigned int TN, unsigned int FN)
{
    return static_cast<float>(TN) / static_cast<float>(TN + FN);
}

double MCC(unsigned int TN, unsigned int FN, unsigned int TP, unsigned int FP)
{
    return ((double)TP*TN - (double)FP*FN) / sqrt((double)(TP + FP)*(double)(TP + FN)*(double)(TN + FP)*(double)(TN + FN));
}

int main(int argc, char const ** argv)
{
    seqan::ArgumentParser parser = buildParser();
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if input was successfully parsed.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    if (isSet(parser, "fp") && isSet(parser, "fr")) {
        std::cout << "comparing files...\n";

        std::string refFilePath;
        getOptionValue(refFilePath, parser, "fr");
        seqan::SeqFileIn refFile;
        if (!open(refFile, refFilePath.c_str(), seqan::OPEN_RDONLY))
        {
            std::cerr << "Error while opening file'" << refFilePath << "'.\n";
            return 1;
        }
        std::string preFilePath;
        getOptionValue(preFilePath, parser, "fp");
        seqan::SeqFileIn preFile;
        if (!open(preFile, preFilePath.c_str(), seqan::OPEN_RDONLY))
        {
            std::cerr << "Error while opening file'" << preFilePath << "'.\n";
            return 1;
        }

        std::map<std::string, unsigned int> refReads;
        std::cout << "Reading reference file... ";
        unsigned int nBases = 0;
        unsigned int maxLen = 0;
        unsigned int n = 0;
        while (!atEnd(refFile))
        {
            std::string id;
            std::string bases;
            readRecord(id, bases, refFile);
            refReads[id] = bases.size();
            if (bases.size() > maxLen)
                maxLen = bases.size();
            nBases += bases.size();
            ++n;
        }
        std::cout << refReads.size() << " reads" << std::endl;
        std::cout << "Comparing reads... ";
        n = 0;
        unsigned int matches = 0;
        unsigned int TN = 0;
        unsigned int TP = 0;
        unsigned int TNb = 0;
        unsigned int TPb = 0;
        unsigned int overtrimmed = 0;
        unsigned int nOvertrimmed = 0;
        unsigned int undertrimmed = 0;
        unsigned int nUndertrimmed = 0;
        while (!atEnd(preFile))
        {
            std::string id;
            std::string bases;
            readRecord(id, bases, preFile);
            const auto it = refReads.find(id);
            if (it->second == bases.size())
            {
                if (it->second == maxLen)
                {
                    ++TN;
                    TNb += it->second;
                }
                else
                {
                    ++TP;
                    TPb += it->second;
                }
                ++matches;
            }
            else if (it->second > bases.size())
            {
                ++overtrimmed;
                nOvertrimmed += static_cast<unsigned int>(it->second - bases.size());
            }
            else if (it->second < bases.size())
            {
                ++undertrimmed;
                nUndertrimmed += static_cast<unsigned int>(bases.size() - it->second);
            }
            it->second = 0;
            ++n;
        }
        auto it = refReads.begin();
        while (it != refReads.end())
        {
            if (it->second > 0)
            {
                ++overtrimmed;
                nOvertrimmed += it->second;
            }
            ++it;
        }

        std::cout << n << " reads" << std::endl;

        std::cout << "\n\nStatistics\n";
        std::cout << "----------\n";
        std::cout << "Number of reads in reference file    : " << refReads.size() << std::endl;
        std::cout << "Number of reads in preprocessed file : " << n << std::endl;

        std::cout << "---------- per read ----------\n";
        assert(TN + TP == matches);
        const auto FP = overtrimmed;
        const auto FN = undertrimmed;
        assert(TN + TP + FN + FP == refReads.size());
        std::cout << "Number of correctly trimmed reads    : " << matches << " (" << (float)(TN+TP) / (float)refReads.size() * 100 << "%)" << std::endl;
        std::cout << "Number of over trimmed reads         : " << overtrimmed << " (" << (float)overtrimmed / (float)refReads.size() * 100 << "%)" << std::endl;
        std::cout << "Number of under trimmed reads        : " << undertrimmed << " (" << (float)undertrimmed / (float)refReads.size() * 100 << "%)" << std::endl;
        std::cout << "Sensitivity_r                        : " << Sen(TP, FN) << std::endl;
        std::cout << "Specificity_r                        : " << Spec(TN, FP) << std::endl;
        std::cout << "PPV_r                                : " << PPV(TP, FP) << std::endl;
        std::cout << "NPV_r                                : " << NPV(TN, FN) << std::endl;
        std::cout << "MCC_r                                : " << MCC(TN, FN, TP, FP) << std::endl;

        std::cout << "---------- per base ----------\n";
        const auto FPb = nOvertrimmed;
        const auto FNb = nUndertrimmed;
        assert(TNb + TPb + FNb + FPb == nBases);

        std::cout << "Number of correctly trimmed bases    : " << matches << " (" << (float)(TNb + TPb) / (float)nBases * 100 << "%)" << std::endl;
        std::cout << "Number of over trimmed bases         : " << nOvertrimmed << " (" << (float)nOvertrimmed / (float)nBases * 100 << "%)" << std::endl;
        std::cout << "Number of under trimmed bases        : " << nUndertrimmed << " (" << (float)nUndertrimmed / (float)nBases * 100 << "%)" << std::endl;
        std::cout << "Sensitivity_b                        : " << Sen(TPb, FNb) << std::endl;
        std::cout << "Specificity_b                        : " << Spec(TNb, FPb) << std::endl;
        std::cout << "PPV_b                                : " << PPV(TPb, FPb) << std::endl;
        std::cout << "NPV_b                                : " << NPV(TNb, FNb) << std::endl;
        std::cout << "MCC_b                                : " << MCC(TNb, FNb, TPb, FPb) << std::endl;

        close(refFile);
        close(preFile);
        return 0;
    }

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
        boost::to_upper(refGenome);
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
            boost::to_upper(adapter);
            adapters.emplace_back(adapter);
        }
        std::cout << "    read " << adapters.size() << " adapters" << std::endl;
    }
    std::cout << std::endl;

    unsigned int numReads = 0;
    if (isSet(parser, "n"))
        getOptionValue(numReads, parser, "n");

    const bool runQualities = isSet(parser, "q");

    std::string fixedBarcode;
    getOptionValue(fixedBarcode, parser, "fb");
    boost::to_upper(fixedBarcode);

    int numRandomBarcode = -1;
    getOptionValue(numRandomBarcode, parser, "rb");

    int readLength = 0;
    getOptionValue(readLength, parser, "rl");

    std::string outPrefix;
    getOptionValue(outPrefix, parser, "o");

    double er = 0;
    getOptionValue(er, parser, "er");

    const unsigned int peakHalfWidthMean = 10;
    const unsigned int peakHalfWidthStdDev = 10;
    const unsigned int peakCoverageMean = 30;
    const unsigned int peakCoverageStdDev = 15;


    std::default_random_engine generator;
    std::normal_distribution<float> peakWidthDistribution((float)peakHalfWidthMean, (float)peakHalfWidthStdDev);
    std::normal_distribution<float> peakCoverageDistribution((float)peakCoverageMean, (float)peakCoverageStdDev);

    unsigned int numPCRArtifacts = 0;
    unsigned int numAdapters = 0;

    seqan::SeqFileOut rawReads;
    open(rawReads, std::string(outPrefix + ".fq").c_str());
    seqan::SeqFileOut preprocessedReads;
    open(preprocessedReads, std::string(outPrefix + "_preprocessed.fq").c_str());
    unsigned int nRead = 0;
    unsigned int nErrors = 0;
    unsigned int nGeneratedBases = 0;
    while(nRead<numReads)
    {
        unsigned int peakPos = rand() % (refGenome.size() - readLength - peakHalfWidthMean*4 + numRandomBarcode + fixedBarcode.size());
        peakPos += peakHalfWidthMean*2;
        unsigned int k = 0;
        unsigned int peakNumReads = (unsigned int)peakCoverageDistribution(generator);
        while(k<peakNumReads && nRead < numReads)
        {
            int posWithinPeak = (unsigned int)peakWidthDistribution(generator); // get peak width
            if (abs(posWithinPeak) > 2*peakHalfWidthMean)
                continue;
            const unsigned int pos = peakPos + posWithinPeak;
            // add fixed and random barcode
            std::string randomBarcode;
            generateRandomBarcode(randomBarcode, numRandomBarcode);
            std::string read = randomBarcode + fixedBarcode + refGenome.substr(pos, readLength - numRandomBarcode - fixedBarcode.size());
            // add adapter
            const auto adapterLength = substituteAdapter(read, adapters);
            
            // insert errors
            nErrors += insertErrors(read, er);

            seqan::Dna5QString temp = read;
            seqan::Dna5QString temp2 = refGenome.substr(pos, readLength - numRandomBarcode - fixedBarcode.size() - adapterLength);
            std::vector<int> qualities;
            if (runQualities)
            {
                qualities = doQualities(temp);
                qualities = decltype(qualities)(qualities.begin() + numRandomBarcode - fixedBarcode.size(), qualities.end());
                assignQualities(temp2, qualities);
            }
            writeRecord(rawReads, std::to_string(nRead), temp);
            writeRecord(preprocessedReads, std::to_string(nRead), temp2);
            nGeneratedBases += read.size();
            ++nRead;
            ++numAdapters;
            // add PCR artifacts
            unsigned int PCRArtifactPercentage = 10;
            if ((unsigned int)(rand() % 100) < PCRArtifactPercentage  && nRead < numReads)
            {
                temp = read;
                if (runQualities)
                {
                    qualities = doQualities(temp);
                    qualities = decltype(qualities)(qualities.begin() + numRandomBarcode - fixedBarcode.size(), qualities.end());
                    assignQualities(temp2, qualities);
                }
                writeRecord(rawReads, std::to_string(nRead) + "_PCR_artifact", temp);
                writeRecord(preprocessedReads, std::to_string(nRead) + "_PCR_artifact", temp2);
                nGeneratedBases += read.size();
                ++numAdapters;
                ++numPCRArtifacts;
                ++nRead;
            }
            ++k;
        }
    }
    std::cout << "Output statistics" << std::endl;
    std::cout << "Number of reads         :\t" << nRead << std::endl;
    std::cout << "Number of PCR-artifacts :\t" << numPCRArtifacts << std::endl;
    std::cout << "Number of errors        :\t" << nErrors << " (" << (float)nErrors/(float)nGeneratedBases  << "%)"<<std::endl;

    close(rawReads);
    close(preprocessedReads);

}
               
