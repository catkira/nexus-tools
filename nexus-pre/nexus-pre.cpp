#include <seqan/bam_io.h>
#include <string>

using namespace seqan;

int main(int argc, char const * argv[])
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " IN.sam\n";
        return 1;
    }
	// Open input file, BamFileIn can read SAM and BAM files.
	BamFileIn bamFileIn(argv[1]);

	// Open output file, BamFileOut accepts also an ostream and a format tag.
	BamFileOut bamFileOut(bamFileIn);
    std::string outFilename = std::string(argv[1]) + std::string("_filtered.sam");
    if (!open(bamFileOut, outFilename.c_str()))
    {
        std::cerr << "ERROR: Could not open " << outFilename << " for writing.\n";
        return 1;
    }

	// Copy header.
	BamHeader header;
	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

    BamAlignmentRecord record;
    // Copy records.
    std::set<std::string> keySet;

	while (!atEnd(bamFileIn))
	{
        readRecord(record, bamFileIn);
        if (record.flag == 0x00 || record.flag == 0x10)
        {
            std::string idString = toCString(record.qName);
            size_t const posStart = idString.find("TL:") + 3;
            if (posStart == std::string::npos)
                continue;
            size_t posEnd = idString.find(':', posStart);
            if (posEnd == std::string::npos)
                posEnd = idString.length();
            std::string const barcode = idString.substr(posStart, posEnd - posStart);
            std::string const key = barcode + ":" + std::to_string(record.rID) + ":" + std::to_string(record.beginPos);
            //std::cout << "barcode: " << barcode << "  key: " << key << std::endl;
            if (keySet.find(key) != keySet.end())
            {
                //std::cout << "found" << std::endl;
            }
            else
            {
                keySet.insert(key);
                writeRecord(bamFileOut, record);
            }
        }
	}

	return 0;
}
