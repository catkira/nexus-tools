#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <string>

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
        "C.  FLEXBAR�Flexible Barcode and Adapter Processing for "
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

    seqan::ArgParseOption recordWriteBam = seqan::ArgParseOption(
        "b", "BAM file", "Create a BAM + BAI file",
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

    unsigned records;
    getOptionValue(records, parser, "r");


    seqan::CharString fileName1, fileName2;
    getArgumentValue(fileName1, parser, 0, 0);

    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn(seqan::toCString(fileName1));

    // Open output file, BamFileOut accepts also an ostream and a format tag.
    BamFileOut bamFileOut(bamFileIn);
    
    std::string outFilename;
    if (fileCount == 2)
    {
        getArgumentValue(fileName2, parser, 0, 1);
        outFilename = seqan::toCString(fileName2);
    }
    else
        outFilename = getFilePrefix(argv[1]) + std::string("_filtered.sam");
    if (!open(bamFileOut, outFilename.c_str()))
    {
        std::cerr << "ERROR: Could not open " << outFilename << " for writing.\n";
        return 1;
    }

    // Copy header.
    BamHeader header;
    readHeader(header, bamFileIn);
    writeHeader(bamFileOut, header);

    BamAlignmentRecord record;
    // Copy records.
    std::set<std::string> keySet;
    std::map<std::string, unsigned> occurenceMapUnique;
    std::map<std::string, unsigned> occurenceMap;

    Statistics stats;

    std::cout << "sorting reads..." << std::endl;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        ++stats.totalReads;
        if (record.flag == 0x00 || record.flag == 0x10)
        {
            ++stats.totalMappedReads;
            std::string idString = toCString(record.qName);
            size_t const posStart = idString.find("TL:") + 3;
            if (posStart == std::string::npos)
                continue;
            size_t posEnd = idString.find(':', posStart);
            if (posEnd == std::string::npos)
                posEnd = idString.length();
            std::string const barcode = idString.substr(posStart, posEnd - posStart);
            // pos = chromosome + strand + position
            std::string const strand = record.flag == 0x00 ? "+" : "-";
            const size_t offset = record.flag == 0x00 ? 0 : seqan::length(record.seq);
            std::string const pos = std::to_string(record.rID) + strand + ":" + std::to_string(record.beginPos+offset);
            std::string const key = barcode + ":" + pos;
            if (keySet.find(key) != keySet.end())
            {
                ++stats.removedReads;
                //std::cout << "found" << std::endl;
            }
            else
            {
                keySet.emplace(std::move(key));
                writeRecord(bamFileOut, record);
                // only count each key once for the occurenceMapUnique 
                ++occurenceMapUnique[pos];
            }
            // every read is stored in occurenceMap
            ++occurenceMap[pos];
        }
    }
    close(bamFileIn);
    close(bamFileOut);

    if (seqan::isSet(parser, "f"))
    {
        BamFileIn bamFileIn2(seqan::toCString(outFilename));
        unsigned clusterSize = 0;
        getOptionValue(clusterSize, parser, "f");
        std::cout << "filtering reads..." << std::endl;
        // Open output file, BamFileOut accepts also an ostream and a format tag.
        setPosition(bamFileIn2, 0);
        BamFileOut bamFileOut2(bamFileIn2);
        const std::string outFilename2 = getFilePrefix(argv[1]) + std::string("_filtered2.sam");
        if (!open(bamFileOut2, outFilename2.c_str()))
        {
            std::cerr << "ERROR: Could not open " << outFilename2 << " for writing.\n";
            return 1;
        }
        // Copy header.
        BamHeader header;
        readHeader(header, bamFileIn2);
        writeHeader(bamFileOut2, header);
        while (!atEnd(bamFileIn2))
        {
            readRecord(record, bamFileIn2);
            if (record.flag != 0x00 && record.flag != 0x10)
                continue;

            std::string idString = toCString(record.qName);
            size_t const posStart = idString.find("TL:") + 3;
            if (posStart == std::string::npos)
                continue;
            size_t posEnd = idString.find(':', posStart);
            if (posEnd == std::string::npos)
                posEnd = idString.length();
            std::string const barcode = idString.substr(posStart, posEnd - posStart);
            // pos = chromosome + strand + position
            std::string const strand = record.flag == 0x00 ? "+" : "-";
            const size_t offset = record.flag == 0x00 ? 0 : seqan::length(record.seq);
            std::string const pos = std::to_string(record.rID) + strand + ":" + std::to_string(record.beginPos + offset);
            const auto& it = occurenceMapUnique.find(pos);
            if (it != occurenceMapUnique.end() && it->second >= clusterSize)
            {
                writeRecord(bamFileOut2, record);
                ++stats.readsAfterFiltering;
            }
        }
        close(bamFileOut2);
    }
    std::fstream fs;
    fs.open(getFilePrefix(argv[1]) + "_statistics.txt", std::fstream::out, _SH_DENYNO);

    // occurenceMap = number of mappings for each location in genome
    // duplicationRate[x] = number of locations with x mappings
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

    std::cout << "calculating unique/non unique duplication Rate..." << std::endl;
    fs << "Duplication Rate unique/non unique" << std::endl;
    unsigned maxLen = size(duplicationRateUnique) > size(duplicationRate) ? size(duplicationRateUnique) : size(duplicationRate);
    duplicationRateUnique.resize(maxLen);
    std::vector<unsigned>::iterator it = duplicationRateUnique.begin();
    for (unsigned i = 0; i < maxLen;++i)
    {
        if (i > 0)
            stats.totalSamePositionReads += (duplicationRate[i] * (i));
        //std::cout << i+1 << " " << duplicationRateUnique[i] << " " << duplicationRate[i] << std::endl;
        fs << i+1 << " " << duplicationRateUnique[i] << " " << duplicationRate[i] << " " 
            << duplicationRateUnique[i]*(i+1) << " " << duplicationRate[i]*(i+1) <<std::endl;
    }
    printStatistics(stats, seqan::isSet(parser, "f"));
    printStatistics(fs, stats, seqan::isSet(parser, "f"));

   
	return 0;
}
