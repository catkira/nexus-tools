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
#include <boost/program_options.hpp>
#include <algorithm>
#include <chrono>

#include "peak.h"
#include "BamRecordKey.h"

namespace po = boost::program_options;

struct Statistics
{
    unsigned totalReads = 0;
    unsigned totalMappedReads = 0;
    unsigned removedReads = 0;
    unsigned totalSamePositionReads = 0;
    unsigned readsAfterFiltering = 0;
    unsigned int couldNotMap = 0;
    unsigned int couldNotMapUniquely = 0;
};

template <typename TStream>
void printStatistics(TStream &stream, const Statistics &stats, const bool clusterFiltering, const bool tabbed = false)
{
    if (tabbed)
    {
        stream << "Total reads" << "\t" << stats.totalReads << std::endl;
        stream << "Mapped reads" << "\t" << stats.totalMappedReads << std::endl;
        stream << "Non mappable reads" << "\t" << stats.couldNotMap << std::endl;
        stream << "Non uniquely mappable reads" << "\t" << stats.couldNotMapUniquely << std::endl;
        stream << "After barcode filtering" << "\t" << stats.totalMappedReads - stats.removedReads << "\t" << " (-" << stats.removedReads << ")" << std::endl;
        stream << "Total duplet reads" << "\t" << stats.totalSamePositionReads << std::endl;
        if (clusterFiltering)
            stream << "After cluster filtering" << "\t" << stats.readsAfterFiltering << "\t" << " (-" << stats.totalMappedReads - stats.removedReads - stats.readsAfterFiltering << ")" << std::endl;
    }
    else
    {
        stream << "Total reads                          : " << stats.totalReads << std::endl;
        stream << "Mapped reads                         : " << stats.totalMappedReads << std::endl;
        stream << "Non mappable reads                   : " << stats.couldNotMap << std::endl;
        stream << "Non uniquely mappable reads          : " << stats.couldNotMapUniquely << std::endl;
        stream << "After barcode filtering              : " << stats.totalMappedReads - stats.removedReads << " (-" << stats.removedReads << ")" << std::endl;
        stream << "Total duplet reads                   : " << stats.totalSamePositionReads << std::endl;
        if (clusterFiltering)
            stream << "After cluster filtering              : " << stats.readsAfterFiltering << " (-" << stats.totalMappedReads - stats.removedReads - stats.readsAfterFiltering << ")" << std::endl;
    }
}

po::options_description buildParser(void)
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("compression", po::value<int>(), "set compression level")
        ("input-file", po::value< std::vector<std::string> >(), "input file")
        ;
    return desc;
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
    auto desc = buildParser();
    
    po::positional_options_description p;
    p.add("input-file", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
        options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    // Check if one or two input files (single or paired-end) were given.
    if (vm.count("input-file") < 1)
    {
        std::cout << desc << "\n";
        return 1;
    }

    if(vm["input-file"].as< std::vector<std::string> >().size() < 2)
    {
        std::cout << desc << "\n";
        return 1;
    }

    const auto fileName1  = vm["input-file"].as< std::vector<std::string> >()[0];
    const auto fileName2 = vm["input-file"].as< std::vector<std::string> >()[1];

    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn(fileName1.c_str());

    std::string outFilename;

    std::set<BamRecordKey<WithBarcode>, CompareBamRecordKey<WithBarcode>> keySet;
    OccurenceMap occurenceMap;

    Statistics stats;

    std::cout << "barcode filtering... ";
    auto t1 = std::chrono::steady_clock::now();
    seqan::BamAlignmentRecord record;

    seqan::BamHeader header;
    readHeader(header, bamFileIn);

    SaveBam<seqan::BamFileIn> saveBam(header, bamFileIn, outFilename);

    unsigned tagID = 0;

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        ++stats.totalReads;
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
        const BamRecordKey<WithBarcode> key(record);
        const BamRecordKey<NoBarcode> pos(record);
        const auto insertResult = keySet.insert(std::move(key));
        OccurenceMap::mapped_type &mapItem = occurenceMap[pos];
        if (!insertResult.second)  // element was not inserted because it existed already
            ++stats.removedReads;
        else
        {
            saveBam.write(record);
            ++mapItem.second; // unique hits
        }
        // non unique hits
        ++mapItem.first;
    }
    saveBam.close();
    seqan::clear(keySet);

    auto t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;

    const std::string outFilename2 = getFilePrefix(argv[1]) + std::string("_filtered2");

    // occurenceMap = number of mappings for each location in genome
    // duplicationRate[x] = number of locations with x mappings
    t1 = std::chrono::steady_clock::now();
    std::cout << "calculating unique/non unique duplication Rate... ";


    seqan::BedRecord<seqan::Bed4> bedRecord;

    std::vector<unsigned> duplicationRateUnique;
    std::vector<unsigned> duplicationRate;
    std::for_each(occurenceMap.begin(), occurenceMap.end(), [&](const OccurenceMap::value_type& val)
    {
        if (duplicationRateUnique.size() < val.second.second)
            duplicationRateUnique.resize(val.second.second);
        ++duplicationRateUnique[val.second.second - 1];
        if (duplicationRate.size() < val.second.first)
            duplicationRate.resize(val.second.first);
        ++duplicationRate[val.second.first - 1];
    });

    std::fstream fs, fs2, fs3;
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
        fs2 << i + 1 << "\t" << duplicationRateUnique[i] * (i + 1) << "\t" << duplicationRate[i] * (i + 1) << std::endl;
    }
    t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    duplicationRateUnique.clear();
    duplicationRate.clear();
    return 0;
}
