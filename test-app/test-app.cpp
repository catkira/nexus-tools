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


auto length(const std::string& str)
{
    return str.length();
}

void resize(std::string& str, const std::string::size_type size)
{
    str.resize(size);
}

template <typename TContainer>
void runBenchmark()
{
    TContainer initSequence = "CAT";
    auto sequence = initSequence;

    auto t1 = std::chrono::steady_clock::now();
    for (unsigned int n = 0;n < 100000000;++n)
        sequence += initSequence;
    auto t2 = std::chrono::steady_clock::now();
    std::cout << "build string: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    auto sequence2 = sequence;

    t1 = std::chrono::steady_clock::now();
    //sequence.erase(std::remove(begin(sequence), end(sequence), 'C'), sequence.end());
    //sequence.erase(std::remove(begin(sequence), end(sequence), 'C'), sequence.end());
    t2 = std::chrono::steady_clock::now();

    auto t1loop = std::chrono::steady_clock::now();

    const auto len = length(sequence2);
    unsigned int keep = 0;
    for (unsigned int n = 0;n < len;n++)
    {
        if (sequence2[n] != 'C')
            sequence2[keep++] = sequence2[n];
    }
    resize(sequence2, keep);

    auto t2loop = std::chrono::steady_clock::now();


    std::cout << "std algorithm" 
        << "\n\ttime = " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" 
        << "\n\tstringLength = " << length(sequence) << std::endl;
    //for (unsigned int n = 0;n < 10;++n)
    //    std::cout << sequence[n];
    //std::cout << std::endl;

    std::cout << "for loop"
        <<"\n\ttime=" << std::chrono::duration_cast<std::chrono::duration<float>>(t2loop - t1loop).count() << "s"
        <<"\n\tstringLength=" << length(sequence) << std::endl;
    //for (unsigned int n = 0;n < 10;++n)
    //    std::cout << sequence2[n];
    //std::cout << std::endl;
}

int main(int argc, char const * argv[])
{
    seqan::CharString str1, str2;
    str1 = "abcd";
    std::cout << "str1: " << str1 << "  str2:" << str2 << std::endl;
    str2 = std::move(str1);
    std::cout << "str1: " << str1 << "  str2:" << str2 << std::endl;


    std::vector<seqan::CharString> stringVector = { "hallo", "welt","ni","hao","hello","world"};

    for (auto& element : stringVector)
        std::cout << element << " ";
    std::cout << std::endl;

    unsigned int keep = 0;
    for (auto& element : stringVector)
        if (element == "ni" || element == "hao")
            stringVector[keep++] = std::move(element);
    resize(stringVector, keep);

    for (auto& element : stringVector)
        std::cout << element << " ";
    std::cout << std::endl;


    return 0;
    std::cout << "std::string\n";
    runBenchmark<std::string>();
    std::cout << "\n\nseqan::CharString\n";
    runBenchmark<seqan::CharString>();
    std::cout << "std::string\n";
    runBenchmark<std::string>();
    return 0;
}
