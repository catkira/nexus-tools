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
    addUsageLine(parser, " \\fI<READ_FILE>  \\fI[OPTIONS]\\fP");
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

    //seqan::ArgParseOption recordWriteBed = seqan::ArgParseOption(
    //    "np", "narrowPeak", "Create a narrow peak file");
    //addOption(parser, recordWriteBed);

    seqan::ArgParseOption ratioOpt = seqan::ArgParseOption(
        "t", "tolerance", "Score ratio tolerance between first and second half Window (1.0 := 100%)",
        seqan::ArgParseOption::DOUBLE, "VALUE");
    setDefaultValue(ratioOpt, 2.0);
    setMinValue(ratioOpt, "0.01");
    setMaxValue(ratioOpt, "100");
    addOption(parser, ratioOpt);

    seqan::ArgParseOption filterChromosomesOpt = seqan::ArgParseOption(
        "fc", "filterChromosomes", "Comma-seperated list of Chromosomes to filter out for calculation of QFragment-Length-Distribution",
        seqan::ArgParseOption::STRING, "LIST");
    setDefaultValue(filterChromosomesOpt, "");
    addOption(parser, filterChromosomesOpt);

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

typedef std::map<BamRecordKey<NoBarcode>, std::pair<unsigned int, unsigned int>> OccurenceMap;

BamRecordKey<NoBarcode> getKey(const OccurenceMap::value_type& val)
{
    return val.first;
}

unsigned getUniqueFrequency(const OccurenceMap::value_type& val)
{
    return val.second.second;
}

template <class T, size_t ROW, size_t COL>
    using Matrix = std::array<std::array<T, COL>, ROW>;

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

    std::string outFilename;
    seqan::CharString output;
    getOptionValue(output, parser, "output");
    if (output != "")
    {
        outFilename = seqan::toCString(output);
    }
    else
        outFilename = getFilePrefix(argv[1]) + std::string("_candidateScores");

    seqan::CharString fileName1;
    getArgumentValue(fileName1, parser, 0, 0);

    seqan::CharString _filterChromosomes;
    seqan::getOptionValue(_filterChromosomes, parser, "fc");
    std::string filterChromosomes = seqan::toCString(_filterChromosomes);

    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn(seqan::toCString(fileName1));

    //const bool narrowPeakEnabled = seqan::isSet(parser, "np");

    OccurenceMap occurenceMap;
    Statistics stats;

    std::cout << "reading file... ";
    auto t1 = std::chrono::steady_clock::now();
    seqan::BamAlignmentRecord record;

    seqan::BamHeader header;
    readHeader(header, bamFileIn);
    const auto chromosomeFilter = calculateChromosomeFilter(filterChromosomes, contigNames(context(bamFileIn)));

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        ++stats.totalReads;
        
        const BamRecordKey<NoBarcode> key(record);
        const auto n = ++occurenceMap[key].second;
        if (n > 1)
            ++stats.samePositionReads;
    }

    auto t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    printStatistics(std::cout, stats);

    
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
    t1 = std::chrono::steady_clock::now();
    
    std::vector<PeakCandidate<OccurenceMap>> peakCandidatesVector;

    const auto range = Range<OccurenceMap>(occurenceMap.begin(), occurenceMap.end());
    auto windowRange = range;
    auto calcScore = [&range, ratioTolerance, halfWindowWidth](const auto _it, auto& _tempSlidingWindowRange)
        {return slidingWindowScore<OccurenceMap>(_it, range, halfWindowWidth, ratioTolerance, _tempSlidingWindowRange);};
    (void)windowRange; // suppress warning

    collectForwardCandidates<OccurenceMap>(range, calcScore, scoreLimit, halfWindowWidth, peakCandidatesVector);
    t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    std::cout << "found " << peakCandidatesVector.size() << " candidates" << std::endl;

    SaveBed<seqan::BedRecord<seqan::Bed4>> saveBedCandidateScores(outFilename);
    saveBedCandidateScores.writeHeader("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
    forwardCandidatesToBed<OccurenceMap, SaveBed<seqan::BedRecord<seqan::Bed4>>, decltype(bamFileIn.context)>(peakCandidatesVector, saveBedCandidateScores, bamFileIn.context);

    //auto filter = [scoreLimit](const PeakCandidate<OccurenceMap>& peakCandidate)
    //    {return peakCandidate.score > 10*scoreLimit ? true : false;};
    //std::map<unsigned int, unsigned int> bindingLengthDistribution;
    //calculateBindingLengthDistribution(peakCandidatesVector, filter, bindingLengthDistribution);
    //std::map<unsigned int, double> bindingCharacteristicsMap;
    const int maxDistance = 1000;
    //calculateScoreDistribution(peakCandidatesVector, calcScore, maxDistance, bindingCharacteristicsMap);
    //auto calcScoreWidth = [&range, ratioTolerance](const auto _it, auto& _tempSlidingWindowRange, auto _halfWindowWidth)
    //{return slidingWindowScore<OccurenceMap>(_it, range, _halfWindowWidth, 0, _tempSlidingWindowRange);};
    //calculateScoreDistribution2(occurenceMap, calcScoreWidth, maxDistance, bindingCharacteristicsMap);

    std::cout << "calculating QFrag-Length-Distribution...";
    t1 = std::chrono::steady_clock::now();
    const auto numChr = seqan::length(contigNames(context(bamFileIn)));
    std::vector<std::vector<unsigned int>> qFragLengthDistribution(maxDistance, std::vector<unsigned int>(numChr));
    calculateQFragLengthDistribution(occurenceMap, qFragLengthDistribution, chromosomeFilter, bamFileIn);
    saveQFragLengthDistribution(getFilePrefix(argv[1]) + "_QFragLengthDistribution.txt", qFragLengthDistribution, bamFileIn);
    t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;

    unsigned int estimatedFragmentLength = 0;
    estimateFragmentLength(qFragLengthDistribution, estimatedFragmentLength);
    std::cout << "estimated fragment length: " << estimatedFragmentLength << std::endl;
    return 0;
}
