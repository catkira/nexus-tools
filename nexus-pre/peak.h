#ifndef PEAK_H_
#define PEAK_H_

#include <algorithm>
#include <vector>

template<typename TContainer>
using Range = std::pair<typename TContainer::const_iterator, typename TContainer::const_iterator>;

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
int slidingWindowScore(typename const TEdgeDistribution::const_iterator centerIt, Range<TEdgeDistribution> range,
    const unsigned widthLimit, Range<TEdgeDistribution>& windowRange)
{
    TEdgeDistribution::const_iterator runningIt = centerIt;
    windowRange.first = centerIt;
    int score = 0;
    int distance = 0;
    if (centerIt == range.first)
        return 0;
    bool ok = false;
    while (runningIt != range.first && (ok = calculateDistance(getKey(*runningIt), getKey(*centerIt), distance)) && distance <= widthLimit)
    {
        if (isReverseStrand(*runningIt))
            score -= getFrequency(*runningIt);
        else
            score += getFrequency(*runningIt);
        windowRange.first = runningIt--;
    }
    if (!ok) // score across chromosomes is not supported, return 0
        return 0;
    runningIt = std::next(centerIt, 1);
    if (runningIt == range.second) // if centerIt is last before chromosome end, return 0 
        return 0;
    windowRange.second = centerIt;
    while (runningIt != range.second && (ok = calculateDistance(getKey(*centerIt), getKey(*runningIt), distance)) && distance <= widthLimit)
    {
        if (isReverseStrand(*runningIt))
            score += getFrequency(*runningIt);
        else
            score -= getFrequency(*runningIt);
        windowRange.second = runningIt++;
    }
    if (!ok) // score across chromosomes is not supported
        return 0;
    return score;    // assume return value optimization
}

template <typename TEdgeDistribution>
void collectForwardCandidates(Range<TEdgeDistribution> range,
    const int scoreLimit, const unsigned widthLimit, std::vector<Range<TEdgeDistribution>>& candidatePositions)
{
    int score = 0;
    int lastScore = 0;
    Range<TEdgeDistribution> slidingWindowRate, tempSlidingWindowRange;
    for (TEdgeDistribution::const_iterator it = range.first; it != range.second; ++it)
    {
        score = slidingWindowScore<TEdgeDistribution>(it, range, widthLimit, tempSlidingWindowRange);
        if (score >= scoreLimit && score > lastScore)
        {
            lastScore = score;
            slidingWindowRate = tempSlidingWindowRange;
            continue;
        }
        if (lastScore > 0)
        {
            candidatePositions.push_back(slidingWindowRate);
            lastScore = 0;
        }
    }
}


#endif  // #ifndef PEAK_H_