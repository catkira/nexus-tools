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
#include <boost/unordered_map.hpp>
#include <algorithm>
#include <chrono>

#include "peak.h"
#include "BamRecordKey.h"
#include <seqan/bam_io.h>

struct Hit
{
    /** \brief The position of a hit corresponds to an 5' end of a fragment. */
    int pos;
    /** \brief A hit is either on the forward strand (T) or the reverse strand (F). */
    bool strand;
    /** \brief True if hit belongs to a qfrag. */
    bool is_q_hit;
    Hit() :
        is_q_hit(false)
    {}
};

struct Chromosome
{
    /** \brief The name of the chromosome e.g. chr1 */
    seqan::CharString name;
    /** \brief The length of the chromosome */
    int len;
    /** \brief Length of the first read in treatment input data */
    int read_len_chip;
    /** \brief The total number of hits on the chromosome for ChIP sample*/
    int hit_num_chip;
    /** \brief The total number of hits on the chromosome for control sample */
    int hit_num_ctrl;
    /** \brief The total number of redundant ChIP hits that were removed*/
    int redundant_hits_removed_num_chip;
    /** \brief The total number of redundant control hits that were removed*/
    int redundant_hits_removed_num_ctrl;
    /** \brief The total number of redundant ChIP hits that were removed*/
    int redundant_hits_num_chip;
    /** \brief The total number of redundant control hits that were removed*/
    int redundant_hits_num_ctrl;
    /** \brief The number of forward strand hits on the chromosome ChIP sample */
    int f_hit_num_chip;
    /** \brief The number of forward strand hits on the chromosome control sample */
    int f_hit_num_ctrl;
    /** \brief The number of reverse strand hits on the chromosome ChIP sample */
    int r_hit_num_chip;
    /** \brief The number of reverse strand hits on the chromosome control sample */
    int r_hit_num_ctrl;
    /** \brief A vector containing all ChIP hits on the chromosome */
    std::vector<Hit> CHIP_HITS;
    /** \brief A vector containing all control hits on the chromosome */
    std::vector<Hit> CTRL_HITS;
    /** \brief The total observed number of qfrags on the chromosome for ChIP sample */
    int qfrag_num_chip;
    /** \brief A vector containing all summits on the chromosome */
    //std::vector<Summit> SUMMITS;
    /** \brief A map containing the hamming distances for all shift sizes */
    std::map<int, int> HAMMING_DISTANCES;
    /** \brief The total number of summits on the chromosome */
    int sum_num;
    /** \brief Relative strand correlation RSC */
    double rsc;

    Chromosome() :
        read_len_chip(-1), hit_num_chip(0), hit_num_ctrl(0), redundant_hits_removed_num_chip(0), redundant_hits_removed_num_ctrl(0), f_hit_num_chip(0), f_hit_num_ctrl(0), r_hit_num_chip(0), r_hit_num_ctrl(0), qfrag_num_chip(0), sum_num(0), rsc(0)
    {}
};

bool compareHitsByPos(const Hit& a, const Hit& b)
{
    return a.pos < b.pos;
}

int RemoveRedundantHits(std::vector<Hit> &HITS, bool keep)
{
    std::vector<Hit> HITS_RMDUP;

    // if the chromosome object has no hits
    if (HITS.size() == 0)
    {
        HITS = HITS_RMDUP;
        return 0;
    }

    int seen_strands;
    int cur_pos = HITS[0].pos;
    HITS_RMDUP.push_back(HITS[0]);
    if (HITS[0].strand == true) { seen_strands = 0; }
    else { seen_strands = 1; }

    for (unsigned int j = 1;j<HITS.size();j++)
    {
        if (cur_pos != HITS[j].pos)
        {
            cur_pos = HITS[j].pos;
            if (HITS[j].strand == true) { seen_strands = 0; }
            else { seen_strands = 1; }
            HITS_RMDUP.push_back(HITS[j]);
        }
        else
        {
            // only f strand has been seen AND hit is on r strand
            if (seen_strands == 0 && !HITS[j].strand)
            {
                HITS_RMDUP.push_back(HITS[j]);
                seen_strands = 2;
            }
            // only r strand has been seen AND hit is on f strand
            else if (seen_strands == 1 && HITS[j].strand)
            {
                HITS_RMDUP.push_back(HITS[j]);
                seen_strands = 2;
            }
        }
    }
    int redundant_hits = HITS.size() - HITS_RMDUP.size();
    if (!keep)
    {
        HITS = HITS_RMDUP;
    }
    return redundant_hits;
}


int ReadAlignmentFile(std::vector<Chromosome> &chromosome, int &chr_num, seqan::CharString chip_sample, seqan::CharString control_sample, bool keep_dup, int thread_num)
{
    // Open input stream, BamStream can read SAM and BAM files
    seqan::BamFileIn bamStreamInChIP;
    if (!open(bamStreamInChIP, toCString(chip_sample)))
    {
        std::cerr << "ERROR: Could not open " << chip_sample << " !\n";
        return 1;
    }

    // create chromosome objects
    seqan::BamHeader header;
    seqan::readHeader(header, bamStreamInChIP);

    const auto& bamContext = seqan::context(bamStreamInChIP);
    chr_num = seqan::length(contigNames(bamContext));

    // maps chip-chromosome name to chip-rID
    boost::unordered_map <std::string, int> chr_map_aux;
    // maps control-chromosome rID to chip-chromosome rID
    boost::unordered_map <int, int> chr_map;

    for (int i = 0;i<chr_num;i++)
    {
        Chromosome c;
        c.name = contigNames(bamContext)[i];
        c.len = contigLengths(bamContext)[i];
        c.hit_num_chip = 0;
        chromosome.push_back(c);
        chr_map_aux[toCString(c.name)] = i;
    }

    // read chip hits to chromosome objects
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStreamInChIP))
    {
        seqan::readRecord(record, bamStreamInChIP);

        // get the read length from the first alignment record
        if (chromosome[record.rID].read_len_chip == -1)
        {
            chromosome[record.rID].read_len_chip = length(record.seq);
        }

        Hit hit;
        if (record.flag == 16)
        {
            hit.pos = record.beginPos + length(record.seq) - 1;
            hit.strand = 1;
            chromosome[record.rID].CHIP_HITS.push_back(hit);
            chromosome[record.rID].hit_num_chip++;
            chromosome[record.rID].r_hit_num_chip++;
        }
        if (record.flag == 0)
        {
            hit.pos = record.beginPos;
            hit.strand = 0;
            chromosome[record.rID].CHIP_HITS.push_back(hit);
            chromosome[record.rID].hit_num_chip++;
            chromosome[record.rID].f_hit_num_chip++;
        }
    }



    // sort hits according to starting position
    SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
        for (int i = 0;i<chr_num;i++)
        {
            sort(chromosome[i].CHIP_HITS.begin(), chromosome[i].CHIP_HITS.end(), compareHitsByPos);
            sort(chromosome[i].CTRL_HITS.begin(), chromosome[i].CTRL_HITS.end(), compareHitsByPos);
        }

    if (keep_dup)
    {
        SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
            for (int i = 0;i<chr_num;i++)
            {
                chromosome[i].redundant_hits_num_chip = RemoveRedundantHits(chromosome[i].CHIP_HITS, true);

                if (control_sample != "None")
                {
                    if (chromosome[i].CTRL_HITS.size() == 0) { continue; }
                    chromosome[i].redundant_hits_num_ctrl = RemoveRedundantHits(chromosome[i].CTRL_HITS, true);
                }
            }
        return 0;
    }

    // if not keep-dup: remove redundant hits
    // --------------------------------------

    SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
        for (int i = 0;i<chr_num;i++)
        {
            chromosome[i].redundant_hits_removed_num_chip = RemoveRedundantHits(chromosome[i].CHIP_HITS, false);
            chromosome[i].redundant_hits_num_chip = chromosome[i].redundant_hits_removed_num_chip;

            // re-determine f_hit_num_chip
            chromosome[i].f_hit_num_chip = 0;
            chromosome[i].r_hit_num_chip = 0;
            for (unsigned int j = 0;j<chromosome[i].CHIP_HITS.size();j++)
            {
                if (chromosome[i].CHIP_HITS[j].strand)
                {
                    chromosome[i].f_hit_num_chip++;
                }
                else
                {
                    chromosome[i].r_hit_num_chip++;
                }
            }
            chromosome[i].hit_num_chip = chromosome[i].f_hit_num_chip + chromosome[i].r_hit_num_chip;
        }

    if (control_sample != "None")
    {
        SEQAN_OMP_PRAGMA(parallel for num_threads(thread_num))
            for (int i = 0;i<chr_num;i++)
            {
                // if there are hits for chip but not for control
                if (chromosome[i].CTRL_HITS.size() == 0) { continue; }

                chromosome[i].redundant_hits_removed_num_ctrl = RemoveRedundantHits(chromosome[i].CTRL_HITS, false);
                chromosome[i].redundant_hits_num_ctrl = chromosome[i].redundant_hits_removed_num_ctrl;

                // re-determine f_hit_num_chip
                chromosome[i].f_hit_num_ctrl = 0;
                chromosome[i].r_hit_num_ctrl = 0;
                for (unsigned int j = 0;j<chromosome[i].CTRL_HITS.size();j++)
                {
                    if (chromosome[i].CTRL_HITS[j].strand)
                    {
                        chromosome[i].f_hit_num_ctrl++;
                    }
                    else
                    {
                        chromosome[i].r_hit_num_ctrl++;
                    }
                }
                chromosome[i].hit_num_ctrl = chromosome[i].f_hit_num_ctrl + chromosome[i].r_hit_num_ctrl;
            }
    }

    return 0;
}

int writeQFragsToBED(Chromosome &chromosome, std::string out_prefix, int min, int max, int &numQFrags)
{
    // fiddle seqan::CharString to std::string
    char* foo = toCString(chromosome.name);
    std::string chr_name = foo;

    int c_hit = 0;
    int pos = 0;

    while (pos < chromosome.len)
    {
        // at position pos are one or more hits
        while (c_hit < chromosome.hit_num_chip && pos == chromosome.CHIP_HITS[c_hit].pos)
        {
            // current hit is on the forward strand
            if (chromosome.CHIP_HITS[c_hit].strand == 0)
            {
                // look for qfrags
                int j = c_hit + 1; // start examination with the next hit

                                   // second hit for qfrag is at most max bases apart
                while (j < chromosome.hit_num_chip && chromosome.CHIP_HITS[j].pos <= (chromosome.CHIP_HITS[c_hit].pos + max))
                {
                    // second hit of the qfrag is at least min bases apart
                    // and on the reverse strand
                    if (chromosome.CHIP_HITS[j].strand == 1 && chromosome.CHIP_HITS[c_hit].pos + min <= chromosome.CHIP_HITS[j].pos)
                    {
                        // write qfrags
                        //out << chromosome.name << "\t" << chromosome.CHIP_HITS[c_hit].pos << "\t" << chromosome.CHIP_HITS[j].pos + 1 << "\n";
                        ++numQFrags;
                    }
                    j++;
                }
            }
            c_hit++;
        }
        pos++;
    }
    return 0;
}


int main(int argc, char const * argv[])
{
    std::vector<Chromosome> chromosomes;
    int num;
    int numQFrags = 0;
    ReadAlignmentFile(chromosomes, num, "P:\\data\\preprocessed\\SRR1175698\\SRR1175698_filtered_sorted.bam", "None", true, 1);
    std::ofstream out;
    out.open("P:\\qfrags.txt");
    for (const auto chromosome : chromosomes)
        out << chromosome.name << "\t";

    for (int n = 1;n < 100;n++)
    {
        numQFrags = 0;
        //for(auto chromosome : chromosomes)
        std::cout << "numQFrags: " << n << ": ";
        for (int i = 0; i < chromosomes.size(); ++i)
        {
            numQFrags = 0;
            writeQFragsToBED(chromosomes[i], "P:\\qfrags", n, n, numQFrags);
            std::cout << numQFrags << "\t";
            out << numQFrags << "\t";
        }
        std::cout << std::endl;
        out << std::endl;

    }
    out.close();
    // Additional checks
//    auto desc = buildParser();
//    
//    po::positional_options_description p;
//    p.add("input-file", -1);
//    po::variables_map vm;
//    po::store(po::command_line_parser(argc, argv).
//        options(desc).positional(p).run(), vm);
//    po::notify(vm);
//
//    if (vm.count("help")) {
//        std::cout << desc << "\n";
//        return 1;
//    }
//
//    // Check if one or two input files (single or paired-end) were given.
//    if (vm.count("input-file") < 1)
//    {
//        std::cout << desc << "\n";
//        return 1;
//    }
//
//    if(vm["input-file"].as< std::vector<std::string> >().size() < 2)
//    {
//        std::cout << desc << "\n";
//        return 1;
//    }
//
//    const auto fileName1  = vm["input-file"].as< std::vector<std::string> >()[0];
//    const auto fileName2 = vm["input-file"].as< std::vector<std::string> >()[1];
//
//    // Open input file, BamFileIn can read SAM and BAM files.
//    seqan::BamFileIn bamFileIn(fileName1.c_str());
//
//    std::string outFilename;
//
//    std::set<BamRecordKey<WithBarcode>, CompareBamRecordKey<WithBarcode>> keySet;
//    OccurenceMap occurenceMap;
//
//    Statistics stats;
//
//    std::cout << "barcode filtering... ";
//    auto t1 = std::chrono::steady_clock::now();
//    seqan::BamAlignmentRecord record;
//
//    seqan::BamHeader header;
//    readHeader(header, bamFileIn);
//
//    SaveBam<seqan::BamFileIn> saveBam(header, bamFileIn, outFilename);
//
//    unsigned tagID = 0;
//
//    while (!atEnd(bamFileIn))
//    {
//        readRecord(record, bamFileIn);
//        ++stats.totalReads;
//        const seqan::BamTagsDict tags(record.tags);
//        if (seqan::findTagKey(tagID, tags, seqan::CharString("XM")))
//        {
//            __int32 tagValue = 0;
//            extractTagValue(tagValue, tags, tagID);
//            if (tagValue == 0)
//                ++stats.couldNotMap;
//            else
//                ++stats.couldNotMapUniquely;
//        }
//        if (record.flag != 0x00 && record.flag != 0x10)
//            continue;
//
//        ++stats.totalMappedReads;
//        const BamRecordKey<WithBarcode> key(record);
//        const BamRecordKey<NoBarcode> pos(record);
//        const auto insertResult = keySet.insert(std::move(key));
//        OccurenceMap::mapped_type &mapItem = occurenceMap[pos];
//        if (!insertResult.second)  // element was not inserted because it existed already
//            ++stats.removedReads;
//        else
//        {
//            saveBam.write(record);
//            ++mapItem.second; // unique hits
//        }
//        // non unique hits
//        ++mapItem.first;
//    }
//    saveBam.close();
//    seqan::clear(keySet);
//
//    auto t2 = std::chrono::steady_clock::now();
//    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
//
//    const std::string outFilename2 = getFilePrefix(argv[1]) + std::string("_filtered2");
//
//    // occurenceMap = number of mappings for each location in genome
//    // duplicationRate[x] = number of locations with x mappings
//    t1 = std::chrono::steady_clock::now();
//    std::cout << "calculating unique/non unique duplication Rate... ";
//
//
//    seqan::BedRecord<seqan::Bed4> bedRecord;
//
//    std::vector<unsigned> duplicationRateUnique;
//    std::vector<unsigned> duplicationRate;
//    std::for_each(occurenceMap.begin(), occurenceMap.end(), [&](const OccurenceMap::value_type& val)
//    {
//        if (duplicationRateUnique.size() < val.second.second)
//            duplicationRateUnique.resize(val.second.second);
//        ++duplicationRateUnique[val.second.second - 1];
//        if (duplicationRate.size() < val.second.first)
//            duplicationRate.resize(val.second.first);
//        ++duplicationRate[val.second.first - 1];
//    });
//
//    std::fstream fs, fs2, fs3;
//#ifdef _MSC_VER
//    fs.open(getFilePrefix(argv[1]) + "_duplication_rate_positions.txt", std::fstream::out, _SH_DENYNO);
//    fs2.open(getFilePrefix(argv[1]) + "_duplication_rate_reads.txt", std::fstream::out, _SH_DENYNO);
//#else
//    fs.open(getFilePrefix(argv[1]) + "_duplication_rate_positions.txt", std::fstream::out);
//    fs2.open(getFilePrefix(argv[1]) + "_duplication_rate_reads.txt", std::fstream::out);
//#endif
//    fs << "rate" << "\t" << "unique" << "\t" << "non unique" << std::endl;
//    fs2 << "rate" << "\t" << "unique" << "\t" << "non unique" << std::endl;
//    unsigned maxLen = duplicationRateUnique.size() > duplicationRate.size() ? duplicationRateUnique.size() : duplicationRate.size();
//    duplicationRateUnique.resize(maxLen);
//    std::vector<unsigned>::iterator it = duplicationRateUnique.begin();
//    for (unsigned i = 0; i < maxLen;++i)
//    {
//        if (i > 0)
//            stats.totalSamePositionReads += (duplicationRate[i] * (i));
//        //std::cout << i+1 << " " << duplicationRateUnique[i] << " " << duplicationRate[i] << std::endl;
//        fs << i + 1 << "\t" << duplicationRateUnique[i] << "\t" << duplicationRate[i] << std::endl;
//        fs2 << i + 1 << "\t" << duplicationRateUnique[i] * (i + 1) << "\t" << duplicationRate[i] * (i + 1) << std::endl;
//    }
//    t2 = std::chrono::steady_clock::now();
//    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
//    duplicationRateUnique.clear();
//    duplicationRate.clear();
    return 0;
}
