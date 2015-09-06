#ifndef BAM_RECORD_KEY_H_
#define BAM_RECORD_KEY_H_

class WithBarcode {};
class NoBarcode {};

template<typename A, typename B>
using disable_if_same_or_derived =
typename std::enable_if<
    !std::is_base_of<A, typename
    std::remove_reference<B>::type
    >::value
>::type;

//template <typename THasBarcode>
//struct CompareBamRecordKey
//{
//    template <typename TBamRecordKey>
//    bool operator()(const TBamRecordKey &lhs, const TBamRecordKey& rhs) const
//    {
//        return lhs.pos < rhs.pos;
//    }
//};
//
//template <>
//struct CompareBamRecordKey<WithBarcode>
//{
//    template <typename TBamRecordKey>
//    bool operator()(const TBamRecordKey &lhs, const TBamRecordKey& rhs) const
//    {
//        if (lhs.pos != rhs.pos)
//            return lhs.pos < rhs.pos;
//        if (lhs.barcode.empty() == false && rhs.barcode.empty() == false
//            && lhs.barcode != rhs.barcode)
//            return lhs.barcode < rhs.barcode;
//        return false;
//    }
//};

bool isRev(const seqan::BamAlignmentRecord &record)
{
    return (record.flag & 0x10) != 0;
}

template <typename THasBarcode = NoBarcode>
struct BamRecordKey
{
    BamRecordKey init(const seqan::BamAlignmentRecord &record)
    {
        pos = static_cast<uint64_t>(record.rID) << 32 |
            static_cast<uint64_t>((record.beginPos + static_cast<uint64_t>((isRev(record) == true ? length(record.seq) : 0)))) << 1 |
            static_cast<uint64_t>(isRev(record));
        return *this;
    }
    BamRecordKey(const uint64_t pos) : pos(pos) {};

    template <typename TRecord, typename X =
        disable_if_same_or_derived<BamRecordKey, TRecord >>
    BamRecordKey(TRecord&& record)
    {
        init(std::forward<TRecord>(record));
    }

    __int32 getPosition() const
    {
        return static_cast<__int32>(pos) >> 1;
    }

    __int32 getRID() const
    {
        return static_cast<__int32>(pos >> 32);
    }

    bool isReverseStrand() const
    {
        return (pos & 0x01) != 0;
    }
private:
    template <typename THasBarcode>
    friend struct CompareBamRecordKey;
    bool friend operator<(const BamRecordKey<NoBarcode>& lhs, const BamRecordKey<NoBarcode>& rhs);
    uint64_t pos;
};

bool operator<(const BamRecordKey<NoBarcode>& lhs, const BamRecordKey<NoBarcode>& rhs)
{
    return lhs.pos < rhs.pos;
}

template <>
struct BamRecordKey<WithBarcode> : BamRecordKey<NoBarcode>
{
    template <typename TRecord>
    BamRecordKey(TRecord&& record) : BamRecordKey<NoBarcode>(record)
    {
        const std::string idString = toCString(record.qName);
        size_t const posStart = idString.find("TL:") + 3;
        if (posStart == std::string::npos)
            return;
        size_t posEnd = idString.find(':', posStart);
        if (posEnd == std::string::npos)
            posEnd = idString.length();
        barcode = idString.substr(posStart, posEnd - posStart);
    }
private:
    std::string barcode;
    bool friend operator<(const BamRecordKey<WithBarcode>& lhs, const BamRecordKey<WithBarcode>& rhs);
};

bool operator<(const BamRecordKey<WithBarcode>& lhs, const BamRecordKey<WithBarcode>& rhs)
{
    if (operator<(static_cast<BamRecordKey<NoBarcode>>(lhs), static_cast<BamRecordKey<NoBarcode>>(rhs)))
        return true;
    if (lhs.barcode.empty() == false && rhs.barcode.empty() == false
        && lhs.barcode != rhs.barcode)
        return lhs.barcode < rhs.barcode;
    return false;
}

// returns false if keys are from different chromosomes
template <typename THasBarcode>
bool calculateDistance(const BamRecordKey<THasBarcode>& key1, const BamRecordKey<THasBarcode>& key2, int& distance)
{
    if (key1.getRID() != key2.getRID())
    {
        //distance = 0;
        return false;
    }
    distance = key2.getPosition() - key1.getPosition();
    return true;
}

#endif
