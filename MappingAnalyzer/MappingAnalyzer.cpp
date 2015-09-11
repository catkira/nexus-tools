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
    unsigned readsFile1 = 0;
    unsigned readsFile2 = 0;
    unsigned matchingReads = 0;
};

template <typename TStream>
void printStatistics(TStream &stream, const Statistics &stats, const bool tabbed = false)
{
    if (tabbed)
    {
        stream << "Reads File1" << "\t" << stats.readsFile1 << std::endl;
        stream << "Reads File2" << "\t" << stats.readsFile2 << std::endl;
        stream << "Matching reads" << "\t" << stats.matchingReads << std::endl;
    }
    else
    {
        stream << "Reads File1     : " << stats.readsFile1 << std::endl;
        stream << "Reads File2     : " << stats.readsFile2 << std::endl;
        stream << "Matching reads  : " << stats.matchingReads << std::endl;
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

template <typename TContext>
struct SaveBam
{
    SaveBam(const seqan::BamHeader header, TContext& context, const std::string& filename)
        : bamFileOut(static_cast<TContext>(context))
    {
        if (!open(bamFileOut, (filename + ".bam").c_str()))
        {
            std::cerr << "ERROR: Could not open " << filename << " for writing.\n";
            return;
        }
        writeHeader(bamFileOut, header);
    }
    void write(const seqan::BamAlignmentRecord& record)
    {
        writeRecord(bamFileOut, record);
    }
    void close()
    {
        seqan::close(bamFileOut);
    }
    seqan::BamFileOut bamFileOut;
};


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

std::string getReadId(const std::string id)
{
    auto pos1 = id.find('.');
    auto pos2 = id.find(':');
    if(pos1 == std::string::npos || pos2 == std::string::npos || pos1 >= pos2)
        return std::string("");
    return id.substr(pos1+1, pos2 - pos1 - 1);
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
    //std::vector<std::pair<unsigned int, std::string>> tempIdStorage;
    std::map<std::string, unsigned int> tempIdStorage;

    SaveBam<seqan::BamFileIn> saveBam(header, bamFileIn1, outPrefix + "_consensus_mapped");

    std::cout << std::endl;

    readRecord(record2, bamFileIn2);
    while (!atEnd(bamFileIn1) && !atEnd(bamFileIn2))
    {
        readRecord(record1, bamFileIn1);
        ++stats.readsFile1;

        if ((stats.readsFile1 + stats.readsFile1) % 100000 == 0)
            std::cout << stats.readsFile1 + stats.readsFile1 << " reads processed" << "\r";

        key1.init(record1);     

        while (!atEnd(bamFileIn2) && lessEqualWithoutStrand(key2.init(record2), key1))
        {
            tempIdStorage.emplace(std::make_pair(getReadId(seqan::toCString(record2.qName)), key2.get5EndPosition()));
            readRecord(record2, bamFileIn2);
            ++stats.readsFile2;
        }
        const auto record1ReadId = getReadId(seqan::toCString(record1.qName));
        const auto it = tempIdStorage.find(record1ReadId);
        if (it != tempIdStorage.end() && it->second == static_cast<int>(key1.get5EndPosition()))
        {
            saveBam.write(record1);
            ++stats.matchingReads;
        }
    }
    saveBam.close();

    auto t2 = std::chrono::steady_clock::now();
    std::cout << std::endl;
    std::cout << "elapsed time: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    printStatistics(std::cout, stats);

    std::fstream fs;
#ifdef _MSC_VER
    fs.open(getFilePrefix(argv[1]) + "_MappingAnalyzer_Statistics.txt", std::fstream::out, _SH_DENYNO);
#else
    fs.open(getFilePrefix(argv[1]) + "_MappingAnalyzer_Statistics.txt", std::fstream::out);
#endif

    fs.close();
    return 0;
}
