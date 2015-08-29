#define SEQAN_PROFILE

#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <algorithm>

#include "peak.h"
#include "BamRecordKey.h"

using seqan::_proFloat;
using seqan::sysTime;

struct Statistics
{
    unsigned totalReads = 0;
    unsigned samePositionReads = 0;
};

template <typename TStream>
void printStatistics(TStream &stream, const Statistics &stats, const bool tabbed=false)
{
    if (tabbed)
    {
        stream << "Total reads" << "\t"<< stats.totalReads << std::endl;
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

    setCategory(parser, "Ting");
    setShortDescription(parser, "Peak Caller for Chip-Nexus and Chip-Exo data");
    addUsageLine(parser, " \\fI<READ_FILE>  \\fI<OUTPUT_FILE> \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "");

    addDescription(parser, "");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "mapped reads", true);
    setValidValues(fileArg, seqan::BamFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "SAM or BAM file");

    seqan::ArgParseOption recordWriteBed = seqan::ArgParseOption(
        "b", "bedGraph", "Create a BedGraph file");
    addOption(parser, recordWriteBed);


    seqan::ArgParseOption performPeakCalling = seqan::ArgParseOption(
        "p", "Peak", "Perform peak calling");
    addOption(parser, performPeakCalling);

    seqan::ArgParseOption ratioOpt = seqan::ArgParseOption(
        "t", "tolerance", "Score ratio tolerance between first and second half Window (1.0 := 100%)",
        seqan::ArgParseOption::DOUBLE, "VALUE");
    setDefaultValue(ratioOpt, 2.0);
    setMinValue(ratioOpt, "0.01");
    setMaxValue(ratioOpt, "100");
    addOption(parser, ratioOpt);

    seqan::ArgParseOption halfWindowSizeOpt = seqan::ArgParseOption(
        "w", "size", "Half window size",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(halfWindowSizeOpt, 20);
    setMinValue(halfWindowSizeOpt, "1");
    addOption(parser, halfWindowSizeOpt);

    seqan::ArgParseOption scoreLimitOpt = seqan::ArgParseOption(
        "s", "score", "Score limit",
        seqan::ArgParseOption::DOUBLE, "VALUE");
    setDefaultValue(scoreLimitOpt, 10);
    setMinValue(scoreLimitOpt, "0.0001");
    addOption(parser, scoreLimitOpt);

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

typedef std::map<BamRecordKey<NoBarcode>, std::pair<unsigned int, unsigned int>, CompareBamRecordKey<NoBarcode>> OccurenceMap;

BamRecordKey<NoBarcode> getKey(const OccurenceMap::value_type& val)
{
    return val.first;
}

unsigned getUniqueFrequency(const OccurenceMap::value_type& val)
{
    return val.second.second;
}

bool isReverseStrand(const OccurenceMap::value_type& val)
{
    return (val.first.pos & 0x01) != 0;
}

template <typename TPosition>
TPosition getPosition(const OccurenceMap::value_type& val)
{
    TPosition position;
    position.chromosomeID = static_cast<__int32>(val.first.pos >> 32);
    position.position = static_cast<__int32>(val.first.pos) >> 1;
    return position;    // assume return value optimization
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
    if (fileCount == 0) 
    {
        printShortHelp(parser);
        return 1;
    }

    seqan::CharString fileName2;
    std::string outFilename;
    if (fileCount == 2)
    {
        getArgumentValue(fileName2, parser, 0, 1);
        outFilename = seqan::toCString(fileName2);
    }
    else
        outFilename = getFilePrefix(argv[1]) + std::string("_candidateScores");

    unsigned numRecords;

    seqan::CharString fileName1;
    getArgumentValue(fileName1, parser, 0, 0);

    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn(seqan::toCString(fileName1));

    const bool bedOutputEnabled = seqan::isSet(parser, "b");

    OccurenceMap occurenceMap;
    auto occurenceMapIt = occurenceMap.begin();

    Statistics stats;

    std::cout << "reading file... ";
    SEQAN_PROTIMESTART(loopTime);
    seqan::BamAlignmentRecord record;

    seqan::BamHeader header;
    readHeader(header, bamFileIn);

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        ++stats.totalReads;
        
        const BamRecordKey<NoBarcode> key(record);
        const auto n = ++occurenceMap[key].second;
        if (n > 1)
            ++stats.samePositionReads;
    }

    double loop = SEQAN_PROTIMEDIFF(loopTime);
    std::cout << loop << "s" << std::endl;
    printStatistics(std::cout, stats);

    SEQAN_PROTIMESTART(peakCandidatesTime);

    double scoreLimit = 0.2;
    unsigned int halfWindowWidth = 30;
    double ratioTolerance = 0.2; // allow 20% tolerance in score between first und second halfWindow
    getOptionValue(scoreLimit, parser, "s");
    getOptionValue(halfWindowWidth, parser, "w");
    getOptionValue(ratioTolerance, parser, "t");

    std::cout << std::endl;
    std::cout << "Settings for peak calling" << std::endl;
    std::cout << "half window size: " << halfWindowWidth << std::endl;
    std::cout << "score limit: " << scoreLimit << std::endl;
    std::cout << "ratio tolerance: " << ratioTolerance << std::endl;

    std::cout << "calculating peak candidates...";
    std::vector<PeakCandidate<OccurenceMap>> positionsVector;

    collectForwardCandidates<OccurenceMap>(Range<OccurenceMap>(occurenceMap.begin(), occurenceMap.end()), scoreLimit, halfWindowWidth, ratioTolerance, positionsVector);
    loop = SEQAN_PROTIMEDIFF(peakCandidatesTime);
    std::cout << loop << "s" << std::endl;
    std::cout << "found " << positionsVector.size() << " candidates" << std::endl;

    SaveBed<seqan::BedRecord<seqan::Bed4>> saveBedCandidateScores(outFilename);
    saveBedCandidateScores.writeHeader("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
    forwardCandidatesToBed<OccurenceMap, SaveBed<seqan::BedRecord<seqan::Bed4>>, decltype(bamFileIn.context)>(positionsVector, saveBedCandidateScores, bamFileIn.context);

	return 0;
}
