#ifndef PEAK_H_
#define PEAK_H_

#include <algorithm>
#include <vector>

struct SingleStrandPosition
{
    const int INVALID_ID = -1;
    const int INVALID_POSITION = -1;

    SingleStrandPosition() :
        chromosomeID(INVALID_ID), position(INVALID_POSITION) {};

    int chromosomeID;
    int position;
};

struct DoubleStrandPosition : SingleStrandPosition
{
    DoubleStrandPosition() :
        SingleStrandPosition(), reverseStrand(false) {};

    bool reverseStrand;
};

template <typename TEdgeDistribution>
unsigned slidingWindowScore(typename TEdgeDistribution::const_iterator startIt, typename TEdgeDistribution::const_iterator endIt,
    const unsigned frequencyLimit, const unsigned widthLimit)
{
    TEdgeDistribution::key_type key = getKey(*startIt);
    while(getKey(*startIt).)
    return 0;
}

template <typename TEdgeDistribution, typename TPosition>
void collectForwardCandidates(typename TEdgeDistribution::const_iterator startIt, typename TEdgeDistribution::const_iterator endIt, 
    const unsigned frequencyLimit, const unsigned widthLimit, std::vector<TPosition>& candidatePositions)
{
    for (TEdgeDistribution::const_iterator it = startIt; it != endIt; ++it)
    {
        if (!isReverseStrand(*it) && getFrequency(*it) > frequencyLimit)
        {
            TEdgeDistribution::const_iterator findResult = std::find_if(next(it, 1), next(it, widthLimit),
                [frequencyLimit](auto element)->bool {return isReverseStrand(element) && (getFrequency(element) > frequencyLimit);});
            if (findResult != next(it, widthLimit))
                candidatePositions.push_back(getPosition<TPosition>(*it));
        }
    }
}


#endif  // #ifndef PEAK_H_