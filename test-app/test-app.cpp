#include <string>
#include <algorithm>
#include <chrono>
#include <future>
#include <iostream>
#include <utility>

using namespace std;

template <class F, class... Ts>
void for_each_argument(F f, Ts&&... a) {
    auto C = std::make_tuple(f(std::forward<Ts>(a))...);
    //(void)std::initializer_list<int>{(f(std::forward<Ts>(a)), 0)...};
}

template<typename Trem, typename TremVal, typename... TContainer>
unsigned int _eraseSeqs(const bool multiThreaded, const Trem& rem, const TremVal remVal, TContainer&&... container)
{
    const auto numRemoveElements = std::count(begin(rem), end(rem), remVal);
    auto eraseElementsWrapper = [&rem, numRemoveElements, remVal, multiThreaded](auto& seq)  // erase Elements using the remove erase idiom
    {
        auto eraseElements = [rem, numRemoveElements, remVal](auto* seq)  // erase Elements using the remove erase idiom
        {
            auto t1 = std::chrono::steady_clock::now();
            const auto beginAddr = (void*)&*begin(*seq);
            std::remove_if(begin(*seq), end(*seq),
                [&rem, &beginAddr, remVal](const auto& element) {
                return rem[&element - beginAddr] == remVal;});
            seq->resize(seq->size() - numRemoveElements);
            auto t2 = std::chrono::steady_clock::now();
            std::cout << "thread id: "<<std::this_thread::get_id() << " time: "<< std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s"  <<std::endl;
        };
        if(multiThreaded)
            return std::async(std::launch::async, eraseElements, &seq);
        eraseElements(&seq);
        return std::future<void>();
    };
    for_each_argument(eraseElementsWrapper, std::forward<TContainer>(container)...);
    return rem.size() - numRemoveElements;
}

template<typename T, T...Is>
struct Index {};

template<typename Trem, typename TremVal, typename TTuple, typename std::size_t... I>
unsigned int _eraseSeqs(const bool multiThreaded, const Trem& rem, const TremVal remVal, TTuple&& tuple, Index<size_t,I...>)
{
    return _eraseSeqs(multiThreaded, rem, remVal, std::get<I>(std::forward<TTuple>(tuple))...);
}

template<typename T, T N, T D>
struct frac
{
    using type = T;
    void print()
    {
        std::cout << nom << "/" << den << std::endl;
    }

    static const int nom = N;
    static const int den = D;
};

template <typename T, typename U>
struct issame
{
};

template <typename T>
struct issame<T, T>
{
    using type = T;
};

template<typename F1, typename F2, typename = issame < F1::type, F2::type >>
struct mult
{
    using res = frac<typename F1::type,F1::nom*F2::nom, F1::den*F2::den>;
};

template<typename F1, typename F2>
struct add
{
    using res = frac<typename F1::type,F1::nom*F2::den + F1::den*F2::nom, F1::den*F2::den>;
};

template<typename T, int N, T...Is>
struct _build_index : _build_index<T, N - 1, N - 1, Is...>
{};

template<typename T, T...Is>
struct _build_index<T, 0, Is...>
{
    using index = Index<T, Is...>;
};

template<typename T, T N>
using build_index = typename _build_index<T, N>::index;


int main(int argc, char const * argv[])
{

    mult<frac<int,2,2>,frac<int,1,2>>::res().print();
    add<frac<int,2, 2>, frac<int,1, 2>>::res().print();
    //return 0;

    std::vector < std::string > stringVector(10000000);
    std::vector<bool> rem(stringVector.size());
    for (unsigned int n = 0;n < stringVector.size();++n)
    {
        stringVector[n] = (n % 2 == 0) ? "even" : "odd";
        rem[n] = n % 2;
    }


    auto test = build_index<size_t,5>{};
    auto test2 = make_integer_sequence<size_t,5>{};

    auto stringVectorTuple = std::make_tuple(stringVector, stringVector, stringVector, stringVector);
    auto index = build_index <size_t, std::tuple_size<decltype(stringVectorTuple)>::value >{};
    for (unsigned int n = 0;n < 10;++n)
        std::cout << std::get<0>(stringVectorTuple)[n] << " ";
    std::cout << std::endl;

    std::cout << "multithreaded" << std::endl;
    auto t1 = std::chrono::steady_clock::now();
    _eraseSeqs(true, rem, true, stringVectorTuple, index);
    auto t2 = std::chrono::steady_clock::now();
    std::cout << "erase time: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;


    stringVectorTuple = std::make_tuple(stringVector, stringVector, stringVector, stringVector);
    std::cout << "singlethreaded" << std::endl;
    t1 = std::chrono::steady_clock::now();
    _eraseSeqs(false, rem, true, stringVectorTuple, index);
    t2 = std::chrono::steady_clock::now();
    std::cout << "erase time: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;

    stringVectorTuple = std::make_tuple(stringVector, stringVector, stringVector, stringVector);
    std::cout << "multithreaded" << std::endl;
    t1 = std::chrono::steady_clock::now();
    _eraseSeqs(true, rem, true, stringVectorTuple, index);
    t2 = std::chrono::steady_clock::now();
    std::cout << "erase time: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;



    for (unsigned int n = 0;n < 10;++n)
        std::cout << std::get<0>(stringVectorTuple)[n] << " ";
    std::cout << std::endl;
    //auto f = std::bind(_eraseSeqs, rem, true);
    //call(f, a);

    return 0;
}
