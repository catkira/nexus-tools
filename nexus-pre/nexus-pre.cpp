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

    seqan::ArgParseOption recordWriteBed = seqan::ArgParseOption(
        "b", "bedGraph", "Create a BedGraph file");
    addOption(parser, recordWriteBed);

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

template <typename TBedRecord>
struct SaveBed
{
    using BedRecord = TBedRecord;
    SaveBed(const std::string& filename) : bedFileOut()
    {
        if (!open(bedFileOut, (filename + ".bed").c_str()))
        {
            std::cerr << "ERROR: Could not open " << filename << " for writing.\n";
            return;
        }
    }
    void write(TBedRecord& record)
    {
        writeRecord(bedFileOut, record);
    }
    void writeHeader(const seqan::CharString& header)
    {
        seqan::write(bedFileOut.iter, header);
    }
    void close()
    {
        seqan::close(bedFileOut);
    }
    seqan::BedFileOut bedFileOut;
};


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


typedef std::map<BamRecordKey<NoBarcode>, unsigned, CompareBamRecordKey<NoBarcode>> OccurenceMapUnique;
typedef std::map<BamRecordKey<NoBarcode>, unsigned, CompareBamRecordKey<NoBarcode>> OccurenceMap;

BamRecordKey<NoBarcode> getKey(const OccurenceMapUnique::value_type& val)
{
    return val.first;
}

unsigned getFrequency(const OccurenceMapUnique::value_type& val)
{
    return val.second;
}

bool isReverseStrand(const OccurenceMapUnique::value_type& val)
{
    return (val.first.pos & 0x01) != 0;
}

template <typename TPosition>
TPosition getPosition(const OccurenceMapUnique::value_type& val)
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
    if (!(fileCount == 1 || fileCount == 2)) {
        printShortHelp(parser);
        return 1;
    }

    unsigned numRecords;
    getOptionValue(numRecords, parser, "r");

    seqan::CharString fileName1, fileName2;
    getArgumentValue(fileName1, parser, 0, 0);

    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn(seqan::toCString(fileName1));
    
    std::string outFilename;
    if (fileCount == 2)
    {
        getArgumentValue(fileName2, parser, 0, 1);
        outFilename = seqan::toCString(fileName2);
    }
    else
        outFilename = getFilePrefix(argv[1]) + std::string("_filtered");

    const bool filter = seqan::isSet(parser, "f");
    const bool sort = seqan::isSet(parser, "s");
    const bool bedOutputEnabled = seqan::isSet(parser, "b");

    typedef std::map<BamRecordKey<NoBarcode>, unsigned, CompareBamRecordKey<NoBarcode>> OccurenceMapUnique;
    typedef std::map<BamRecordKey<NoBarcode>, unsigned, CompareBamRecordKey<NoBarcode>> OccurenceMap;
    std::set<BamRecordKey<WithBarcode>, CompareBamRecordKey<WithBarcode>> keySet;
    OccurenceMapUnique occurenceMapUnique;
    OccurenceMap occurenceMap;

    Statistics stats;

    std::cout << "barcode filtering... ";
    SEQAN_PROTIMESTART(loopTime);
    seqan::BamAlignmentRecord record;

    seqan::BamHeader header;
    readHeader(header, bamFileIn);
    
    SaveBam<seqan::BamFileIn> saveBam(header, bamFileIn, outFilename);

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        ++stats.totalReads;
        if (record.flag == 0x00 || record.flag == 0x10)
        {
            ++stats.totalMappedReads;
            const BamRecordKey<WithBarcode> key(record);
            const BamRecordKey<NoBarcode> pos(record);

            const auto insertResult = keySet.insert(std::move(key));

            if (!insertResult.second)  // element was not inserted because it existed already
            {
                ++stats.removedReads;
            }
            else
            {
                saveBam.write(record);
                ++occurenceMapUnique[pos];
            }
            // every read is stored in occurenceMap
            ++occurenceMap[pos];
        }
    }
    saveBam.close();
    seqan::clear(keySet);

    double loop = SEQAN_PROTIMEDIFF(loopTime);
    std::cout << loop << "s" << std::endl;

    const std::string outFilename2 = getFilePrefix(argv[1]) + std::string("_filtered2");
    if (filter)
    {
        SEQAN_PROTIMESTART(loopTime);
        std::cout << "filtering reads... ";
        unsigned clusterSize = 0;
        getOptionValue(clusterSize, parser, "f");
        // Open output file, BamFileOut accepts also an ostream and a format tag.
        //while (!atEnd(bamFileIn2))
        seqan::BamFileIn bamFileIn2(seqan::toCString(outFilename + ".bam"));
        //clear(header);
        readHeader(header, bamFileIn2);
        SaveBam<seqan::BamFileIn> saveBam2(header, bamFileIn2, outFilename2);
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

    SaveBed<seqan::BedRecord<seqan::Bed4>> saveBedForwardStrand(outFilename + "_forward");
    SaveBed<seqan::BedRecord<seqan::Bed4>> saveBedReverseStrand(outFilename + "_reverse");
    if (bedOutputEnabled)
    {
        saveBedForwardStrand.writeHeader("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
        saveBedReverseStrand.writeHeader("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
    }

    seqan::BedRecord<seqan::Bed4> bedRecord;

    std::vector<unsigned> duplicationRateUnique;
    std::for_each(occurenceMapUnique.begin(), occurenceMapUnique.end(), [&](const OccurenceMapUnique::value_type& val)
    {
        if (bedOutputEnabled)
        {
            bedRecord.rID = static_cast<int32_t>(val.first.pos >> 32);
            bedRecord.ref = contigNames(bamFileIn.context)[bedRecord.rID];
            if (val.first.pos & 0x01)    // reverse strand
            {
                // I think this -1 is not neccessary, but its here to reproduce the data from the CHipNexus paper exactly
                bedRecord.beginPos = (static_cast<int32_t>(val.first.pos) >> 1) - 1; 
                bedRecord.endPos = bedRecord.beginPos + 1;
                bedRecord.name = std::to_string(-static_cast<int32_t>(val.second)); // abuse name as val parameter in BedGraph
                saveBedReverseStrand.write(bedRecord);
            }
            else    // forward strand
            {
                bedRecord.beginPos = (static_cast<int32_t>(val.first.pos) >> 1);
                bedRecord.endPos = bedRecord.beginPos + 1;
                bedRecord.name = std::to_string(val.second); // abuse name as val parameter in BedGraph
                saveBedForwardStrand.write(bedRecord);
            }
        }
        if (duplicationRateUnique.size() < val.second)
            duplicationRateUnique.resize(val.second);
        ++duplicationRateUnique[val.second - 1];
    });
    saveBedForwardStrand.close();
    saveBedReverseStrand.close();

    std::vector<unsigned> duplicationRate;
    std::for_each(occurenceMap.begin(), occurenceMap.end(), [&](const OccurenceMapUnique::value_type& val)
    {
        if (duplicationRate.size() < val.second)
            duplicationRate.resize(val.second);
        ++duplicationRate[val.second - 1];
    });

    std::fstream fs,fs2,fs3;
#ifdef _MSC_VER
    fs.open(getFilePrefix(argv[1]) + "_duplication_rate_positions.txt", std::fstream::out, _SH_DENYNO);
    fs2.open(getFilePrefix(argv[1]) + "_duplication_rate_reads.txt", std::fstream::out, _SH_DENYNO);
#else
    fs.open(getFilePrefix(argv[1]) + "_duplication_rate_positions.txt", std::fstream::out);
    fs2.open(getFilePrefix(argv[1]) + "_duplication_rate_reads.txt", std::fstream::out);
#endif
    fs << "rate" << "\t" << "unique" << "\t" << "non unique" << std::endl;
    fs2 << "rate" << "\t" << "unique" << "\t" << "non unique" << std::endl;
    unsigned maxLen = duplicationRateUnique.size() > duplicationRate.size() ? duplicationRateUnique.size() : duplicationRate.size();
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
    duplicationRateUnique.clear();
    duplicationRate.clear();

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

        seqan::BamFileIn bamFileIn3(seqan::toCString(inFilename2));
        readHeader(header, bamFileIn3);
        std::vector<seqan::BamAlignmentRecord> records;

        while (!atEnd(bamFileIn3))
        {
            readRecord(record, bamFileIn3);
            records.emplace_back(std::move(record));
        }
        close(bamFileIn3);
        std::sort(records.begin(), records.end(), [](const seqan::BamAlignmentRecord& lhs, const seqan::BamAlignmentRecord& rhs) {return CompareBamRecordKey<NoBarcode>()(BamRecordKey<NoBarcode>(lhs), BamRecordKey<NoBarcode>(rhs));});

        const std::string outFilename2 = getFilePrefix(inFilename2) + std::string("_sorted.bam");
        SaveBam<seqan::BamFileIn> saveBam3(header, bamFileIn, outFilename2);
        for (auto record : records)
            saveBam3.write(record);
        saveBam3.close();


        loop = SEQAN_PROTIMEDIFF(sortTime);
        std::cout << loop << "s" << std::endl;
    }

    printStatistics(stats, seqan::isSet(parser, "f"));
#ifdef _MSV_VER
    fs3.open(getFilePrefix(argv[1]) + "_statistics.txt", std::fstream::out, _SH_DENYNO);
#else
    fs3.open(getFilePrefix(argv[1]) + "_statistics.txt", std::fstream::out);
#endif
    printStatistics(fs3, stats, seqan::isSet(parser, "f"));

    SEQAN_PROTIMESTART(peakCandidatesTime);
    std::cout << "calculating peak candidates...";
    std::vector<PeakCandidate<OccurenceMapUnique>> positionsVector;
    const int scoreLimit = 40;
    collectForwardCandidates<OccurenceMapUnique>(Range<OccurenceMapUnique>(occurenceMapUnique.begin(), occurenceMapUnique.end()), scoreLimit, 50, positionsVector);
    loop = SEQAN_PROTIMEDIFF(peakCandidatesTime);
    std::cout << loop << "s" << std::endl;
    std::cout << "found " << positionsVector.size() << " candidates" << std::endl;

    SaveBed<seqan::BedRecord<seqan::Bed4>> saveBedCandidateScores(outFilename + "_candidateScores");
    saveBedCandidateScores.writeHeader("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
    forwardCandidatesToBed<OccurenceMapUnique, SaveBed<seqan::BedRecord<seqan::Bed4>>, decltype(bamFileIn.context)>(positionsVector, saveBedCandidateScores, bamFileIn.context);

	return 0;
}
