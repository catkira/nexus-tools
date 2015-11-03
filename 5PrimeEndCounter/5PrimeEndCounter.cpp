#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <string>
#include <algorithm>
#include <chrono>
#include <cassert>

#include "peak.h"
#include "BamRecordKey.h"

struct Statistics
{
    unsigned totalReads = 0;
    unsigned totalMappedReads = 0;
    unsigned filteredReads = 0;
    unsigned removedReads = 0;
    unsigned totalSamePositionReads = 0;
    unsigned readsAfterFiltering = 0;
    unsigned int couldNotMap = 0;
    unsigned int couldNotMapUniquely = 0;
};

seqan::ArgumentParser buildParser(void)
{
    seqan::ArgumentParser parser;

    setCategory(parser, "5-prime end counter");
    setShortDescription(parser, "Preprocessing Pipeline for Chip-Nexus and Chip-Exo data");
    addUsageLine(parser, " \\fI<READ_FILE1> \\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "");

    addDescription(parser, "");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "positions in bed-format", true);
    setValidValues(fileArg, ".bed");
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "bed file");

    seqan::ArgParseOption readsOpt = seqan::ArgParseOption(
        "i", "input", "Name of the reads file.",
        seqan::ArgParseOption::INPUT_FILE, "INPUT");
    setValidValues(readsOpt, seqan::BamFileIn::getFileExtensions());
    addOption(parser, readsOpt);

    seqan::ArgParseOption radiusOpt = seqan::ArgParseOption(
        "r", "radius", "radius around peaks to scan ",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(radiusOpt, 1000);
    setMinValue(radiusOpt, "1");
    addOption(parser, radiusOpt);

    seqan::ArgParseOption filterChromosomesOpt = seqan::ArgParseOption(
        "fc", "filterChromosomes", "Regular expression to remove chromosomes",
        seqan::ArgParseOption::STRING, "REGEX");
    addOption(parser, filterChromosomesOpt);

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

typedef std::map<BamRecordKey<NoBarcode>, unsigned int> OccurenceMap;

BamRecordKey<NoBarcode> getKey(const OccurenceMap::value_type& val)
{
    return val.first;
}

template <typename TOccurenceMap, typename TChromosomeFilter>
void processBamFile(seqan::BamFileIn& bamFileIn, const TChromosomeFilter& chromosomeFilter, TOccurenceMap &occurenceMap, Statistics& stats)
{
    seqan::BamAlignmentRecord record;
    unsigned tagID = 0;
    std::set<BamRecordKey<WithBarcode>> keySet;

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        if (atEnd(bamFileIn))
            break;
        ++stats.totalReads;
        if (chromosomeFilter.find(record.rID) != chromosomeFilter.end())
        {
            ++stats.filteredReads;
            continue;
        }
        const seqan::BamTagsDict tags(record.tags);
        if (seqan::findTagKey(tagID, tags, seqan::CharString("XM")))
        {
            __int32 tagValue = 0;
            extractTagValue(tagValue, tags, tagID);
            if (tagValue == 0)
                ++stats.couldNotMap;
            else
                ++stats.couldNotMapUniquely;
        }
        if (record.flag != 0x00 && record.flag != 0x10)
            continue;

        ++stats.totalMappedReads;
        const BamRecordKey<NoBarcode> pos(record);
        OccurenceMap::mapped_type &mapItem = occurenceMap[pos];
        ++mapItem;
    }
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
    if (fileCount < 1) {
        printShortHelp(parser);
        return 1;
    }

    unsigned int radius = 1;
    getOptionValue(radius, parser, "r");

    seqan::CharString readsFileName;
    getOptionValue(readsFileName, parser, "i");

    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn(seqan::toCString(readsFileName));
    
    std::string outFilename;
    outFilename = getFilePrefix(argv[1]) + std::string("_filtered");

    seqan::CharString _filterChromosomes;
    seqan::getOptionValue(_filterChromosomes, parser, "fc");
    std::string filterChromosomes = seqan::toCString(_filterChromosomes);

    OccurenceMap occurenceMap;
    Statistics stats;

    std::cout << "read bam file... ";
    auto t1 = std::chrono::steady_clock::now();
    seqan::BamAlignmentRecord record;
    seqan::BamHeader header;
    readHeader(header, bamFileIn);
    const auto chromosomeFilterSet = calculateChromosomeFilter(filterChromosomes, contigNames(context(bamFileIn)));
    const auto chromosomes = contigNames(context(bamFileIn));
    processBamFile(bamFileIn, chromosomeFilterSet, occurenceMap, stats);
    auto t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;

    std::vector<std::pair<unsigned int, unsigned int>> hits(radius * 2 + 1);

    t1 = std::chrono::steady_clock::now();
    std::cout << "calculating 5'-ends around peaks... ";

    for (unsigned int fileIndex = 0;fileIndex < static_cast<unsigned int>(fileCount); ++fileIndex)
    {
        seqan::CharString fileName_;
        getArgumentValue(fileName_, parser, fileIndex, 0);
        const std::string fileName = seqan::toCString(fileName_);

        std::ifstream infile(fileName);
        std::string chromosome, dummy;
        unsigned int start, end;
        while (infile >> chromosome >> start >> end >> dummy)
        {
            int rID = -1;
            for (unsigned int i = 0;i < length(chromosomes);++i)
                if (chromosomes[i] == chromosome)
                {
                    rID = i;
                    break;
                }
            if (rID == -1)
            {
                std::cout << "invalid chromosome name: " << chromosome << " in file " << fileName << std::endl;
                return -1;
            }
            seqan::BamAlignmentRecord record;
            record.beginPos = std::max<int>(start - radius, 0);
            record.rID = rID;
            record.flag = 0;
            unsigned int index = 0;
            if (start < radius)
                index += radius - start;
            while (record.beginPos <= static_cast<__int32>(start + radius))
            {
                BamRecordKey<NoBarcode> pos(record);
                hits[index].first += occurenceMap[pos];
                pos.init(pos.getRID(), pos.get5EndPosition(), true);
                hits[index].second += occurenceMap[pos];
                ++record.beginPos;
                ++index;
            }
        }

        std::fstream fs;
        const std::string outFilename = getFilePrefix(fileName) + std::string("_5PrimeEnds.tab");
        std::cout << "writing " << outFilename << std::endl;
#ifdef _MSC_VER
        fs.open(outFilename, std::fstream::out, _SH_DENYNO);
#else
        fs.open(outFilename, std::fstream::out);
#endif
        int i = - static_cast<int>(radius);
        for (const auto& hit : hits)
            fs << i++ << "\t" << hit.first << "\t" << hit.second << std::endl;
        fs.close();
    }
    t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;


	return 0;
}
