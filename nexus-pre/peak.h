#ifndef PEAK_H_
#define PEAK_H_

#include <algorithm>
#include <vector>
#include <functional>

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
    PeakCandidate() : score(0) {};
    void clear()
    {
        score = 0;
    }
    Range<TEdgeDistribution> range;
    typename TEdgeDistribution::const_iterator centerIt;
    double score;
};

template <typename TEdgeDistribution>
double slidingWindowScore(const typename TEdgeDistribution::const_iterator centerIt, Range<TEdgeDistribution> range,
    const unsigned halfWindowWidth, const double scoreLimit, const double ratioTolerance, Range<TEdgeDistribution>& windowRange)
{
    typename TEdgeDistribution::const_iterator runningIt = centerIt;
    windowRange.first = centerIt;
    double score = 0;
    double half1score = 0;
    int distance = 0;
    if (centerIt == range.first)
        return 0;
    bool ok = false;
    while (runningIt != range.first && (ok = calculateDistance(getKey(*runningIt), getKey(*centerIt), distance)) && distance <= halfWindowWidth)
    {
        if (isReverseStrand(*runningIt))
            score -= getUniqueFrequency(*runningIt);
        else
            score += getUniqueFrequency(*runningIt);
        windowRange.first = runningIt--;
    }
    if (!ok) // score across chromosomes is not supported, return 0
        return 0;
    runningIt = std::next(centerIt, 1);
    if (runningIt == range.second) // if centerIt is last before chromosome end, return 0 
        return 0;
    half1score = score;
    windowRange.second = centerIt;
    while (runningIt != range.second && (ok = calculateDistance(getKey(*centerIt), getKey(*runningIt), distance)) && distance <= halfWindowWidth)
    {
        if (isReverseStrand(*runningIt))
            score += getUniqueFrequency(*runningIt);
        else
            score -= getUniqueFrequency(*runningIt);
        windowRange.second = runningIt++;
    }
    double ratio = static_cast<double>(score - half1score) / static_cast<double>(half1score);
    if (ratio < (1 - ratioTolerance) || ratio >(1 + ratioTolerance))
    {
        //if (getPosition(getKey(*centerIt)) > 2456589 && getPosition(getKey(*centerIt)) < 2456631)
        //    std::cout << "ratio fail: " << ratio <<std::endl;
        return 0;
    }
    if (!ok) // score across chromosomes is not supported
        return 0;
    return score;    // assume return value optimization
}

template <typename TEdgeDistribution, typename TLambda>
void plateauAdjustment(PeakCandidate<TEdgeDistribution>& peakCandidate, TLambda& calcScore)
{
    unsigned int plateauCount = 0;
    Range<TEdgeDistribution> tempSlidingWindowRange;
    auto prevIt = peakCandidate.centerIt;
    for (typename TEdgeDistribution::const_iterator it = std::next(peakCandidate.centerIt,1); it != peakCandidate.range.second; prevIt = it++)
    {
        //if (getPosition(getKey(*prevIt)) == getPosition(getKey(*it)))
        //    continue;
        if (calcScore(it, tempSlidingWindowRange) > (peakCandidate.score * 0.9))
            ++plateauCount;
        else
            break;
    }
    if (plateauCount > 0)
        peakCandidate.centerIt = std::next(peakCandidate.centerIt, plateauCount / 2);
}

template <typename TEdgeDistribution>
void collectForwardCandidates(const Range<TEdgeDistribution> range,
    const double scoreLimit, const unsigned halfWindowWidth, const double ratioTolerance, typename std::vector<PeakCandidate<TEdgeDistribution>>& candidatePositions)
{
    double tempScore = 0;
    int checkAhead = 0;
    int plateauCount = 0;
    auto calcScore = [&range, halfWindowWidth, scoreLimit, ratioTolerance](const auto _it, auto& _tempSlidingWindowRange) 
        {return slidingWindowScore<TEdgeDistribution>(_it, range, halfWindowWidth, scoreLimit, ratioTolerance, _tempSlidingWindowRange);};
    PeakCandidate<TEdgeDistribution> peakCandidate;
    Range<TEdgeDistribution> tempSlidingWindowRange;
    auto prevIt = range.first;
    for (typename TEdgeDistribution::const_iterator it = range.first; it != range.second; prevIt=it++)
    {
        //if (prevIt != it && getPosition(getKey(*prevIt)) == getPosition(getKey(*it)))   // check if the next object is on the same position, but only another strand
        //    continue;
        tempScore = calcScore(it, tempSlidingWindowRange);
        if (tempScore >= scoreLimit && peakCandidate.score == 0)    // scan until first match
        {
            checkAhead = halfWindowWidth;    // start checkAhead
            peakCandidate.score = tempScore;
            peakCandidate.range = tempSlidingWindowRange;
            peakCandidate.centerIt = it;
        }
        if(checkAhead > 0)  // after first match, checkAhead for better peak candidates
        {
            if (tempScore > peakCandidate.score)  // is the current candidate better
            {
                peakCandidate.score = tempScore;
                peakCandidate.range = tempSlidingWindowRange;
                peakCandidate.centerIt = it;
            }
            --checkAhead;
            continue;
        }
        if (peakCandidate.score > 0)    // checking finished
        {
            if (getPosition(getKey(*peakCandidate.centerIt)) > 2456589 && getPosition(getKey(*peakCandidate.centerIt)) < 2456631)
            {
                std::cout << "\npeak score: " << peakCandidate.score << std::endl;
                std::cout << "pos before plateau adjustment: " << getPosition(getKey(*peakCandidate.centerIt)) << std::endl;
            }
            plateauAdjustment(peakCandidate, calcScore);
            if (getPosition(getKey(*peakCandidate.centerIt)) > 2456589 && getPosition(getKey(*peakCandidate.centerIt)) < 2456631)
                std::cout << "pos after plateau adjustment: " << getPosition(getKey(*peakCandidate.centerIt)) << std::endl;
            candidatePositions.push_back(peakCandidate);
            it = peakCandidate.range.second;
            tempScore = 0;
            checkAhead = 0;
            peakCandidate.clear();
        }
    }
}

template <typename TEdgeDistribution, typename TWriter, typename TContext>
void forwardCandidatesToBed(const typename std::vector<PeakCandidate<TEdgeDistribution>>& candidatePositions, TWriter& writer, TContext& context)
{
    typename TWriter::BedRecord bedRecord;
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
