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
struct PeakCandidate
{
    Range<TEdgeDistribution> range;
    typename TEdgeDistribution::const_iterator centerIt;
    int score;
};

template <typename TEdgeDistribution>
int slidingWindowScore(typename const TEdgeDistribution::const_iterator centerIt, Range<TEdgeDistribution> range,
    const unsigned widthLimit, const int halfScoreLimit, Range<TEdgeDistribution>& windowRange)
{
    TEdgeDistribution::const_iterator runningIt = centerIt;
    windowRange.first = centerIt;
    int score = 0;
    int half1score = 0;
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
    if (score < halfScoreLimit) // prevent matches where there is only a rising in the 2nd half
        return 0;
    half1score = score;
    windowRange.second = centerIt;
    while (runningIt != range.second && (ok = calculateDistance(getKey(*centerIt), getKey(*runningIt), distance)) && distance <= widthLimit)
    {
        if (isReverseStrand(*runningIt))
            score += getFrequency(*runningIt);
        else
            score -= getFrequency(*runningIt);
        windowRange.second = runningIt++;
    }
    if (score - half1score < halfScoreLimit)
        return 0;
    double ratio = static_cast<double>(score - half1score) / static_cast<double>(half1score);
    if (ratio < 0.3 || ratio > 0.7)
        return 0;
    if (!ok) // score across chromosomes is not supported
        return 0;
    return score;    // assume return value optimization
}

template <typename TEdgeDistribution>
void collectForwardCandidates(Range<TEdgeDistribution> range,
    const int scoreLimit, const unsigned widthLimit, std::vector<PeakCandidate<TEdgeDistribution>>& candidatePositions)
{
    int score = 0;
    int lastScore = 0;
    PeakCandidate<TEdgeDistribution> peakCandidate;
    Range<TEdgeDistribution> slidingWindowRange, tempSlidingWindowRange;
    for (TEdgeDistribution::const_iterator it = range.first; it != range.second; ++it)
    {
        score = slidingWindowScore<TEdgeDistribution>(it, range, widthLimit, scoreLimit/2, tempSlidingWindowRange);
        if (score >= scoreLimit && score > lastScore)
        {
            lastScore = score;
            slidingWindowRange = tempSlidingWindowRange;
            continue;
        }
        if (lastScore > 0)
        {
            peakCandidate.centerIt = it;
            peakCandidate.range = range;
            peakCandidate.score = score;
            candidatePositions.push_back(peakCandidate);
            it = slidingWindowRange.second;
            lastScore = 0;
        }
    }
}

template <typename TEdgeDistribution, typename TWriter, typename TContext>
void forwardCandidatesToBed(const std::vector<PeakCandidate<TEdgeDistribution>>& candidatePositions, TWriter& writer, TContext& context)
{
    TWriter::BedRecord bedRecord;
    for (const auto& element : candidatePositions)
    {
        bedRecord.rID = static_cast<int32_t>(element.centerIt->first.pos >> 32);
        bedRecord.ref = contigNames(context)[bedRecord.rID];    // ADL
        bedRecord.beginPos = (static_cast<int32_t>(element.centerIt->first.pos) >> 1) - 1;
        bedRecord.endPos = bedRecord.beginPos + 1;
        bedRecord.name = std::to_string(element.score); // abuse name as score parameter in BedGraph

        writer.write(bedRecord);
    }
}


#endif  // #ifndef PEAK_H_