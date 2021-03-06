// ==========================================================================
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================

#ifndef GENERALSTATS_H
#define GENERALSTATS_H

template <typename TLen>
struct AdapterTrimmingStats
{
    using LenType = TLen;

    std::vector<std::vector<TLen>> removedLength;
    std::vector<TLen> numRemoved;
    TLen overlapSum;
    TLen minOverlap, maxOverlap;

    AdapterTrimmingStats() : overlapSum(0),
        minOverlap(std::numeric_limits<TLen>::max()), maxOverlap(0) {};

    AdapterTrimmingStats& operator+= (AdapterTrimmingStats const& rhs)
    {
        overlapSum += rhs.overlapSum;
        minOverlap = minOverlap < rhs.minOverlap ? minOverlap : rhs.minOverlap;
        maxOverlap = maxOverlap < rhs.maxOverlap ? rhs.maxOverlap : maxOverlap;
        {
            const auto len = rhs.removedLength.size();
            if (removedLength.size() < len)
                removedLength.resize(std::max(removedLength.size(), len));
            for (TLen i = 0;i < len;++i)
            {
                const auto len2 = rhs.removedLength[i].size();
                if (removedLength[i].size() < len2)
                    removedLength[i].resize(len2);
                for (TLen k = 0;k < len2;++k)
                    removedLength[i][k] += rhs.removedLength[i][k];
            }
        }

        {
            const auto len = rhs.numRemoved.size();
            if (numRemoved.size() < len)
                numRemoved.resize(len);
            for (TLen i = 0;i < len;++i)
                numRemoved[i] += rhs.numRemoved[i];
        }
        return *this;
    }
    void clear()
    {
        overlapSum = 0;
        minOverlap = std::numeric_limits<TLen>::max();
        maxOverlap = 0;
        auto it = numRemoved.begin();
        while (it != numRemoved.end())
        {
            *it = 0;
            ++it;
        }
    }
};

struct GeneralStats
{
    unsigned removedN;       //Number of deleted sequences due to N's
    unsigned removedDemultiplex;
    unsigned removedQuality;
    unsigned long uncalledBases;//Number of uncalled bases (evtl. Masked) in surviving sequences
    unsigned removedShort;  //Number of deleted sequences due to shortness.
    unsigned int readCount;
    float processTime;
    float readTime;
    float writeTime;
    std::vector<unsigned int> matchedBarcodeReads;

    using TAdapterTrimmingStats = AdapterTrimmingStats<unsigned int>;
    TAdapterTrimmingStats adapterTrimmingStats;

    void clear()
    {
        removedN = removedDemultiplex = removedQuality = uncalledBases = removedShort = readCount = 0;
        processTime = readTime = writeTime = 0;
        auto it = matchedBarcodeReads.begin();
        while (it != matchedBarcodeReads.end())
        {
            *it = 0;
            ++it;
        }
        adapterTrimmingStats.clear();
    };

    GeneralStats(): removedN(0), removedDemultiplex(0), removedQuality(0), uncalledBases(0), removedShort(0), readCount(0), processTime(0), readTime(0), writeTime(0) {};
    GeneralStats(unsigned int N, unsigned int numAdapters) : GeneralStats() 
    { 
        matchedBarcodeReads.resize(N); 
        adapterTrimmingStats.numRemoved.resize(numAdapters);
    };
    GeneralStats(const GeneralStats& rhs) = default;
    GeneralStats(GeneralStats&& rhs) = default;

    GeneralStats& operator=(const GeneralStats& rhs) = default;
    GeneralStats& operator=(GeneralStats&& rhs) = default;


    GeneralStats& operator+=(const GeneralStats& rhs)
    {
        removedN += rhs.removedN;
        removedDemultiplex += rhs.removedDemultiplex;
        removedQuality += rhs.removedQuality;
        uncalledBases += rhs.uncalledBases;
        removedShort += rhs.removedShort;
        readCount += rhs.readCount;
        processTime += rhs.processTime;
        readTime += rhs.readTime;
        writeTime += rhs.writeTime;
        if (matchedBarcodeReads.size() != rhs.matchedBarcodeReads.size())
            matchedBarcodeReads.resize(rhs.matchedBarcodeReads.size());
        matchedBarcodeReads = matchedBarcodeReads + rhs.matchedBarcodeReads;
        adapterTrimmingStats += rhs.adapterTrimmingStats;
        return *this;
    }
};

#endif
