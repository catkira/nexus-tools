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
int slidingWindowScore(typename const TEdgeDistribution::const_iterator centerIt, typename TEdgeDistribution::const_iterator startIt, typename TEdgeDistribution::const_iterator endIt,
    const unsigned widthLimit)
{
    TEdgeDistribution::const_iterator runningIt = centerIt;
    TEdgeDistribution::const_iterator prevRunningIt = runningIt;
    int score = 0;
    int distance = 0;
    if (centerIt == startIt)
        return 0;
    bool ok = false;
    while (runningIt != startIt && (ok = calculateDistance(getKey(*runningIt), getKey(*centerIt), distance)) && distance <= widthLimit)
    {
        if (isReverseStrand(*runningIt))
            score -= getFrequency(*runningIt);
        else
            score += getFrequency(*runningIt);
        prevRunningIt = runningIt--;
    }
    if (!ok) // score across chromosomes is not supported, return 0
        return 0;
    runningIt = std::next(centerIt, 1);
    if (runningIt == endIt) // if centerIt is last before chromosome end, return 0 
        return 0;
    prevRunningIt = centerIt;
    while (runningIt != endIt && (ok = calculateDistance(getKey(*centerIt), getKey(*runningIt), distance)) && distance <= widthLimit)
    {
        if (isReverseStrand(*runningIt))
            score += getFrequency(*runningIt);
        else
            score -= getFrequency(*runningIt);
        prevRunningIt = runningIt++;
    }
    if (!ok) // score across chromosomes is not supported
        return 0;
    return score;    // assume return value optimization
}

template <typename TEdgeDistribution, typename TPosition>
void collectForwardCandidates(typename TEdgeDistribution::const_iterator startIt, typename TEdgeDistribution::const_iterator endIt, 
    const int scoreLimit, const unsigned widthLimit, std::vector<TPosition>& candidatePositions)
{
    for (TEdgeDistribution::const_iterator it = startIt; it != endIt; ++it)
    {
        if (slidingWindowScore<TEdgeDistribution>(it, startIt, endIt, widthLimit) >= scoreLimit)
            candidatePositions.push_back(getPosition<TPosition>(*it));
    }
}


#endif  // #ifndef PEAK_H_