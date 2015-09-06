#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <string>
#include <algorithm>
#include <chrono>

#include "peak.h"
#include "BamRecordKey.h"

struct Statistics
{
    unsigned totalReads = 0;
    unsigned samePositionReads = 0;
};

template <typename TStream>
void printStatistics(TStream &stream, const Statistics &stats, const bool tabbed = false)
{
    if (tabbed)
    {
        stream << "Total reads" << "\t" << stats.totalReads << std::endl;
        stream << "Same position reads" << "\t" << stats.samePositionReads << std::endl;
    }
    else
    {
        stream << "Total reads                          : " << stats.totalReads << std::endl;
        stream << "Same position reads                  : " << stats.samePositionReads << std::endl;
    }
}

seqan::ArgumentParser buildParser(void)
{
    seqan::ArgumentParser parser;

    setCategory(parser, "MappingAnalyzer");
    setShortDescription(parser, "");
    addUsageLine(parser, " \\fI<READ_FILE1> \\fI<READ_FILE2> \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "");

    addDescription(parser, "");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "mapped reads", true);
    setValidValues(fileArg, seqan::BamFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "SAM or BAM file");

    seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
        "o", "output", "Prefix of the output file.",
        seqan::ArgParseOption::OUTPUT_FILE, "OUTPUT");
    setDefaultValue(outputOpt, "");
    addOption(parser, outputOpt);

    return parser;
}

std::string getFilePath(const std::string& fileName)
{
    std::size_t found = fileName.find_last_of("/\\");
    if (found == std::string::npos || found < 1)
        return std::string();
    return fileName.substr(0, found);
}

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

std::string getFilename(const std::string& fileName)
{
    std::size_t found = fileName.find_last_of("/\\");
    if (found == std::string::npos)
        return std::string();
    return fileName.substr(found + 1);
}



int main(int argc, char const * argv[])
{
    // Additional checks
    seqan::ArgumentParser parser = buildParser();
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if input was successfully parsed.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check if one or two input files (single or paired-end) were given.
    int fileCount = getArgumentValueCount(parser, 0);
    if (fileCount < 2)
    {
        printShortHelp(parser);
        return 1;
    }

    std::string outPrefix;
    seqan::CharString output;
    getOptionValue(output, parser, "output");
    if (output != "")
    {
        outPrefix = seqan::toCString(output);
    }
    else
        outPrefix = getFilePrefix(argv[1]);

    seqan::CharString fileName1;
    getArgumentValue(fileName1, parser, 0, 0);
    seqan::CharString fileName2;
    getArgumentValue(fileName2, parser, 0, 1);

    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn1(seqan::toCString(fileName1));
    seqan::BamFileIn bamFileIn2(seqan::toCString(fileName2));

    typedef std::vector<BamRecordKey<NoBarcode>> ReadsSeen;
    ReadsSeen readsSeen1, readsSeen2;
    Statistics stats;

    std::cout << "reading " << getFilename(seqan::toCString(fileName1)) << " ...";
    auto t1 = std::chrono::steady_clock::now();
    seqan::BamAlignmentRecord record1, record2;

    seqan::BamHeader header;
    readHeader(header, bamFileIn1);
    readHeader(header, bamFileIn2);

    BamRecordKey<NoBarcode> key1(record1);
    BamRecordKey<NoBarcode> key2(record2);

    while (!atEnd(bamFileIn1))
    {
        readRecord(record1, bamFileIn1);
        ++stats.totalReads;

        key1.init(record1);
        while(!atEnd(bamFileIn1) && key2.init(record2) < key1)
            readRecord(record2, bamFileIn2);
        if (key1 == key2)
        {
            // check if readId1 == readId2
            // if same, write read into sameMapped
        }
        else
        {
            // write read into notMappedIn
        }
        readsSeen1.emplace_back(key1);
    }

    auto t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;

    std::fstream fs;
#ifdef _MSC_VER
    fs.open(getFilePrefix(argv[1]) + "MappingAnalyzer_Statistics.txt", std::fstream::out, _SH_DENYNO);
#else
    fs.open(getFilePrefix(argv[1]) + "MappingAnalyzer_Statistics.txt", std::fstream::out);
#endif

    fs.close();
    return 0;
}
