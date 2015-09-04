#ifndef BAM_RECORD_KEY_H_
#define BAM_RECORD_KEY_H_

class WithBarcode {};
class NoBarcode {};

template <typename THasBarcode>
struct CompareBamRecordKey
{
    template <typename TBamRecordKey>
    bool operator()(const TBamRecordKey &lhs, const TBamRecordKey& rhs) const
    {
        return lhs.pos < rhs.pos;
    }
};

template <>
struct CompareBamRecordKey<WithBarcode>
{
    template <typename TBamRecordKey>
    bool operator()(const TBamRecordKey &lhs, const TBamRecordKey& rhs) const
    {
        if (lhs.pos != rhs.pos)
            return lhs.pos < rhs.pos;
        if (lhs.barcode.empty() == false && rhs.barcode.empty() == false
            && lhs.barcode != rhs.barcode)
            return lhs.barcode < rhs.barcode;
        return false;
    }
};

bool isRev(const seqan::BamAlignmentRecord &record)
{
    return (record.flag & 0x10) != 0;
}

template <typename THasBarcode = NoBarcode>
struct BamRecordKey
{
    BamRecordKey(const uint64_t pos) : pos(pos) {};
    BamRecordKey(const seqan::BamAlignmentRecord &record)
    {
        pos = static_cast<uint64_t>(record.rID) << 32 |
            static_cast<uint64_t>((record.beginPos + static_cast<uint64_t>((isRev(record) == true ? length(record.seq) : 0)))) << 1 |
            static_cast<uint64_t>(isRev(record));
    };
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
    typedef CompareBamRecordKey<WithBarcode> TCompareBamRecordKey;
private:
    template <typename THasBarcode>
    friend struct CompareBamRecordKey;
    uint64_t pos;
};


template <>
struct BamRecordKey<WithBarcode> : BamRecordKey<NoBarcode>
{
    BamRecordKey(const seqan::BamAlignmentRecord &record) : BamRecordKey<NoBarcode>(record)
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
    typedef CompareBamRecordKey<NoBarcode> TCompareBamRecordKey;
    std::string barcode;
};

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
