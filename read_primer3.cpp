#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "read_primer3.hpp"
#include "libprimervg.hpp"

/*
Terrible read in - but gets the job done

Populates the vgprimers object using the primer3
output file.
*/
int read_primer_entry(std::ifstream* file, vgprimers* primers)
{
    std::string line;
    bool entry = false;
    /* Readlines until file or current entry ends */
    while (std::getline(*file, line) && line != "=")
    {
        /* lines in TAG=DATUM format, entry ends with '=' */
        int split;
        if ((split=line.find('=')) == -1)
        {
            /* Not inherently breaking, but print warning */
            std::cout << "Illegal line found in primer3 output" << std::endl;
            continue;
        }
        /* Isolate the tag and data of the line */
        entry = true; /* At least one set of primers was generated */
        std::string tag = line.substr(0, split);
        std::string datum = line.substr(split+1, line.size()-(split+1));

        if ("SEQUENCE_ID" == tag)
        {
            primers->name = datum;
        }
        else if ("SEQUENCE_TEMPLATE" == tag)
        {
            primers->sequence = datum;
        }
        else if ("PRIMER_LEFT_" == tag.substr(0, 12))
        {
            /* Get pair number (0-based) */
            int end_num = tag.find_last_of("0123456789");
            int npair = atoi(tag.substr(12, end_num-11).c_str());
            /* Ensure enough primer pairs are allocated */
            while (npair >= primers->primer_pairs.size())
                primers->primer_pairs.push_back(primerpair());
            /* Save sequences of each primer for output */
            if ("_SEQUENCE" == tag.substr(end_num+1, tag.size()-end_num))
            {
                primers->primer_pairs[npair].left.sequence = datum;
            }
            /* Save coordinates of primer */
            else if (tag.length() == end_num+1)
            {
                int separator = datum.find_first_of(',');
                int l_pos = atoi( datum.substr(0, separator).c_str() );
                int len = atoi( datum.substr(separator+1, datum.length()-(separator)).c_str() );
                primers->primer_pairs[npair].left.position = l_pos;
                primers->primer_pairs[npair].left.length = len;
            }
        }
        else if ("PRIMER_RIGHT_" == tag.substr(0, 13))
        {
            /* Get pair number (0-based) */
            int end_num = tag.find_last_of("0123456789");
            int npair = atoi(tag.substr(13, end_num-12).c_str());
            /* Ensure enough primer pairs are allocated */
            while (npair >= primers->primer_pairs.size())
                primers->primer_pairs.push_back(primerpair());
            /* Save sequences of each primer for output */
            if ("_SEQUENCE" == tag.substr(end_num+1, tag.size()-end_num))
            {
                primers->primer_pairs[npair].right.sequence = datum;
            }
            /* Save coordinates of primer */
            else if (tag.length() == end_num+1)
            {
                int separator = datum.find_first_of(',');
                int l_pos = atoi( datum.substr(0, separator).c_str() );
                int len = atoi( datum.substr(separator+1, datum.length()-(separator)).c_str() );
                primers->primer_pairs[npair].right.position = l_pos;
                primers->primer_pairs[npair].right.length = len;
            }
        }
    }
    /* Sanity checks for seq, name, and primers */
    if (primers->name.empty() || !primers->primer_pairs.size() || primers->sequence.empty())
    {
        std::cout << "Bad entry read, incomplete data" << std::endl;
        return false;
    }
    /* Successful entry found? */
    return entry;
}