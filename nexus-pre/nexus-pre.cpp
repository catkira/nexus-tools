#define SEQAN_PROFILE

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <algorithm>

using namespace seqan;


struct Statistics
{
    unsigned totalReads = 0;
    unsigned totalMappedReads = 0;
    unsigned removedReads = 0;
    unsigned totalSamePositionReads = 0;
    unsigned readsAfterFiltering = 0;

};

void printStatistics(const Statistics &stats, const bool clusterFiltering)
{
    std::cout << "Total reads                          : " << stats.totalReads << std::endl;
    std::cout << "Total mapped reads                   : " << stats.totalMappedReads << std::endl;
    std::cout << "Total duplet reads                   : " << stats.totalSamePositionReads << std::endl;
    std::cout << "After barcode filtering              : " << stats.totalMappedReads - stats.removedReads << " (-" << stats.removedReads << ")" <<std::endl;
    if (clusterFiltering)
        std::cout << "After cluster filtering              : " << stats.readsAfterFiltering << " (-" << stats.totalMappedReads - stats.removedReads - stats.readsAfterFiltering << ")" << std::endl;
}

// Todo: remove bad copypaste code
void printStatistics(std::fstream &fs, const Statistics &stats, const bool clusterFiltering)
{
    fs << "Total reads                          : " << stats.totalReads << std::endl;
    fs << "Total mapped reads                   : " << stats.totalMappedReads << std::endl;
    fs << "Total duplet reads                   : " << stats.totalSamePositionReads << std::endl;
    fs << "After barcode filtering              : " << stats.totalMappedReads - stats.removedReads << " (-" << stats.removedReads << ")" << std::endl;
    if (clusterFiltering)
        fs << "After cluster filtering              : " << stats.readsAfterFiltering << " (-" << stats.totalMappedReads - stats.removedReads - stats.readsAfterFiltering << ")" << std::endl;
}

seqan::ArgumentParser buildParser(void)
{
    seqan::ArgumentParser parser;

    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Filtering Toolkit of seqan_flexbar.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of SeqAn-Flexbar (a reimplementation of"
        " the original flexbar[1]) and can be used to filter reads and apply "
        "sequence independent trimming options");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
        "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
        "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "mapped reads", true);
    setValidValues(fileArg, seqan::BamFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "SAM or BAM file");

    seqan::ArgParseOption recordFilterCluster = seqan::ArgParseOption(
        "f", "cluster size", "Minimum number of mapped reads at the same position",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(recordFilterCluster, 0);
    setMinValue(recordFilterCluster, "0");
    addOption(parser, recordFilterCluster);

    seqan::ArgParseOption recordSortBam = seqan::ArgParseOption(
        "s", "sort", "Sort BAM file");
    addOption(parser, recordSortBam);

    seqan::ArgParseOption recordWriteBam = seqan::ArgParseOption(
        "b", "BAI file", "Create a BAI file (activates the sort option)",
        seqan::ArgParseArgument::INPUT_FILE, "BAM/BAI filename");
    addOption(parser, recordWriteBam);

    seqan::ArgParseOption recordOpt = seqan::ArgParseOption(
        "r", "records", "Number of records to be read in one run.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(recordOpt, 10000);
    setMinValue(recordOpt, "10");
    addOption(parser, recordOpt);

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

class WithBarcode {};
class NoBarcode {};

bool isRev(const BamAlignmentRecord &record)
{
    return (record.flag & 0x10) != 0;
};

template <typename THasBarcode>
struct CompareBamRecordKey
{
    template <typename TBamRecordKey>
    bool operator()(const TBamRecordKey &lhs, const TBamRecordKey& rhs) const
    {
        return lhs.pos < rhs.pos;
    }
};

template <>
struct CompareBamRecordKey<WithBarcode>
{
    template <typename TBamRecordKey>
    bool operator()(const TBamRecordKey &lhs, const TBamRecordKey& rhs) const
    {
        if (lhs.pos != rhs.pos)
            return lhs.pos < rhs.pos;
        if (lhs.barcode.empty() == false && rhs.barcode.empty() == false
            && lhs.barcode != rhs.barcode)
            return lhs.barcode < rhs.barcode;
        return false;
    }
};

template <typename TContext>
struct SaveBam
{
    template <typename TContext>
    SaveBam(const BamHeader header, const TContext& context, const std::string& filename)
        : bamFileOut((BamFileOut::TDependentContext)context)
    {
        if (!open(bamFileOut, filename.c_str()))
        {
            std::cerr << "ERROR: Could not open " << filename << " for writing.\n";
            return;
        }
        writeHeader(bamFileOut, header);
    }
    void write(const BamAlignmentRecord& record)
    {
        writeRecord(bamFileOut, record);
    }
    void close()
    {
        seqan::close(bamFileOut);
    }
    BamFileOut bamFileOut;
};

template <typename THasBarcode>
struct BamRecordKey
{
    BamRecordKey(const uint64_t pos) : pos(pos){};
    BamRecordKey(const BamAlignmentRecord &record)
    {
        pos = (uint64_t)record.rID << 32 | (record.beginPos + (isRev(record) == true ? length(record.seq) : 0)) << 1 | (uint64_t)isRev(record);
    };
    typedef CompareBamRecordKey<WithBarcode> TCompareBamRecordKey;
    uint64_t pos;
};

template <>
struct BamRecordKey<WithBarcode> : BamRecordKey<NoBarcode>
{
    BamRecordKey(const BamAlignmentRecord &record) : BamRecordKey<NoBarcode>(record)
    {
        const std::string idString = toCString(record.qName);
        size_t const posStart = idString.find("TL:") + 3;
        if (posStart == std::string::npos)
            return;
        size_t posEnd = idString.find(':', posStart);
        if (posEnd == std::string::npos)
            posEnd = idString.length();
        barcode = idString.substr(posStart, posEnd - posStart);
    }
    typedef CompareBamRecordKey<NoBarcode> TCompareBamRecordKey;
    std::string barcode;
};

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
    if (!(fileCount == 1 || fileCount == 2)) {
        printShortHelp(parser);
        return 1;
    }

    unsigned numRecords;
    getOptionValue(numRecords, parser, "r");

    seqan::CharString fileName1, fileName2;
    getArgumentValue(fileName1, parser, 0, 0);

    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn(seqan::toCString(fileName1));
    
    std::string outFilename;
    if (fileCount == 2)
    {
        getArgumentValue(fileName2, parser, 0, 1);
        outFilename = seqan::toCString(fileName2);
    }
    else
        outFilename = getFilePrefix(argv[1]) + std::string("_filtered.bam");

    const bool filter = seqan::isSet(parser, "f");
    const bool sort = seqan::isSet(parser, "s");

    //BamAlignmentRecord record;
    // Copy records.
    //std::set<std::string> keySet;
    std::set<BamRecordKey<WithBarcode>, CompareBamRecordKey<WithBarcode>> keySet;
    std::map<BamRecordKey<NoBarcode>, unsigned, CompareBamRecordKey<NoBarcode>> occurenceMapUnique;
    std::map<BamRecordKey<NoBarcode>, unsigned, CompareBamRecordKey<NoBarcode>> occurenceMap;
    seqan::String<BamAlignmentRecord> records;

    Statistics stats;

    std::cout << "sorting reads... ";
    SEQAN_PROTIMESTART(loopTime);
    BamAlignmentRecord record;
    bool writeOutputStage1 = !filter;

    BamHeader header;
    readHeader(header, bamFileIn);
    
    SaveBam<BamFileOut::TDependentContext> saveBam(header, bamFileIn.context, outFilename);

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        ++stats.totalReads;
        if (record.flag == 0x00 || record.flag == 0x10)
        {
            ++stats.totalMappedReads;
            const BamRecordKey<WithBarcode> key(record);
            const BamRecordKey<NoBarcode> pos(record);
            if (keySet.find(key) != keySet.end())
            {
                ++stats.removedReads;
            }
            else
            {
                saveBam.write(record);
                keySet.insert(key);
                ++occurenceMapUnique[pos];
            }
            // every read is stored in occurenceMap
            ++occurenceMap[pos];
        }
    }
    saveBam.close();
    clear(keySet);

    double loop = SEQAN_PROTIMEDIFF(loopTime);
    std::cout << loop << "s" << std::endl;

    const std::string outFilename2 = getFilePrefix(argv[1]) + std::string("_filtered2.bam");
    if (filter)
    {
        SEQAN_PROTIMESTART(loopTime);
        std::cout << "filtering reads... ";
        unsigned clusterSize = 0;
        getOptionValue(clusterSize, parser, "f");
        // Open output file, BamFileOut accepts also an ostream and a format tag.
        //while (!atEnd(bamFileIn2))
        unsigned ex = 0;
        BamFileIn bamFileIn2(seqan::toCString(outFilename));
        readHeader(header, bamFileIn2);
        SaveBam<BamFileOut::TDependentContext> saveBam2(header, bamFileIn.context, outFilename2);
        while (!atEnd(bamFileIn2))
        {
            readRecord(record, bamFileIn2);
            const auto occurenceMapUniqueIt =  occurenceMapUnique.find(BamRecordKey<NoBarcode>(record));
            if (occurenceMapUniqueIt->second >= clusterSize)
            {
                //sortedBamVector2.emplace_back(sortedBamVector[i]);
                saveBam2.write(record);
                ++stats.readsAfterFiltering;
                //++keyMapIt;
            }
        }
        saveBam2.close();

        loop = SEQAN_PROTIMEDIFF(loopTime);
        std::cout << loop << "s" << std::endl;

    }

    // occurenceMap = number of mappings for each location in genome
    // duplicationRate[x] = number of locations with x mappings
    SEQAN_PROTIMESTART(finalProcessing);
    std::cout << "calculating unique/non unique duplication Rate... ";
    std::vector<unsigned> duplicationRateUnique;
    std::for_each(occurenceMapUnique.begin(), occurenceMapUnique.end(), [&](auto &it)
    {
        if (duplicationRateUnique.size() < it.second)
            duplicationRateUnique.resize(it.second);
        ++duplicationRateUnique[it.second - 1];
    });

    std::vector<unsigned> duplicationRate;
    std::for_each(occurenceMap.begin(), occurenceMap.end(), [&](auto &it)
    {
        if (duplicationRate.size() < it.second)
            duplicationRate.resize(it.second);
        ++duplicationRate[it.second - 1];
    });
    clear(occurenceMap);
    clear(occurenceMapUnique);


    std::fstream fs,fs2,fs3;
    fs.open(getFilePrefix(argv[1]) + "_duplication_rate_positions.txt", std::fstream::out, _SH_DENYNO);
    fs << "rate" << "\t" << "unique" << "\t" << "non unique" << std::endl;
    fs2.open(getFilePrefix(argv[1]) + "_duplication_rate_reads.txt", std::fstream::out, _SH_DENYNO);
    fs2 << "rate" << "\t" << "unique" << "\t" << "non unique" << std::endl;
    unsigned maxLen = size(duplicationRateUnique) > size(duplicationRate) ? size(duplicationRateUnique) : size(duplicationRate);
    duplicationRateUnique.resize(maxLen);
    std::vector<unsigned>::iterator it = duplicationRateUnique.begin();
    for (unsigned i = 0; i < maxLen;++i)
    {
        if (i > 0)
            stats.totalSamePositionReads += (duplicationRate[i] * (i));
        //std::cout << i+1 << " " << duplicationRateUnique[i] << " " << duplicationRate[i] << std::endl;
        fs << i + 1 << "\t" << duplicationRateUnique[i] << "\t" << duplicationRate[i] << std::endl;
        fs2 << i + 1 << "\t" << duplicationRateUnique[i]*(i+1) << "\t" << duplicationRate[i]*(i+1) <<std::endl;
    }
    loop = SEQAN_PROTIMEDIFF(finalProcessing);
    std::cout << loop << "s" << std::endl;
    clear(duplicationRateUnique);
    clear(duplicationRate);

    if (sort)
    {
        SEQAN_PROTIMESTART(sortTime);
        std::cout << "sorting BAM file... ";
        unsigned clusterSize = 0;
        getOptionValue(clusterSize, parser, "f");
        // Open output file, BamFileOut accepts also an ostream and a format tag.
        //while (!atEnd(bamFileIn2))
        std::string inFilename2;
        if (filter)
            inFilename2 = outFilename2;
        else
            inFilename2 = outFilename;

        BamFileIn bamFileIn3(seqan::toCString(inFilename2));
        readHeader(header, bamFileIn3);
        std::vector<BamAlignmentRecord> records;

        while (!atEnd(bamFileIn3))
        {
            readRecord(record, bamFileIn3);
            records.emplace_back(std::move(record));
        }
        close(bamFileIn3);
        std::sort(records.begin(), records.end(), [](auto &lhs, auto& rhs) {return CompareBamRecordKey<NoBarcode>()(BamRecordKey<NoBarcode>(lhs), BamRecordKey<NoBarcode>(rhs));});

        const std::string outFilename2 = getFilePrefix(inFilename2) + std::string("_sorted.bam");
        SaveBam<BamFileOut::TDependentContext> saveBam3(header, bamFileIn.context, outFilename2);
        for (auto record : records)
            saveBam3.write(record);
        saveBam3.close();


        loop = SEQAN_PROTIMEDIFF(sortTime);
        std::cout << loop << "s" << std::endl;
    }

    printStatistics(stats, seqan::isSet(parser, "f"));
    fs3.open(getFilePrefix(argv[1]) + "_statistics.txt", std::fstream::out, _SH_DENYNO);
    printStatistics(fs3, stats, seqan::isSet(parser, "f"));

	return 0;
}
