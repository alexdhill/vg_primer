#include <string>
#include <iostream>
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <gbwtgraph/gbz.h>
#include <bdsg/snarl_distance_index.hpp>

#include "libprimervg.hpp"
#include "read_primer3.hpp"
#include "vgfilters.hpp"

void print_vgprimers(vgprimers *primers);
void print_help();

int main(int argc, char* argv[]) // ./vg_primers .p3 .gam .dist
{
    /* Too simple for real option handling */
    if (argc != 4)
    {
        print_help();
        return (-1);
    }

    /* Open primer3 output */
    std::ifstream infile(argv[1]);
    if (infile.fail())
    {
        std::cout << "Failed to load P3 output" << std::endl;
        return (-1);
    }
    /* Init our primers (to be filled) */
    vgprimers* primers;

    /* Right now assumes single entry - easy to exand to higher throughput */
    primers = new vgprimers();
    if (!read_primer_entry(&infile, primers))
    {
        std::cout << "Failed to read P3 output" << std::endl;
        return(-1);
    }

    /* Read in aligned sequence, graph and index */
    vg::Path seq_path;
    std::vector<vg::Path> gam_out;

    /* Prep gam file to read */
    std::string gam_path(argv[2]);
    auto gam_file = std::ifstream(gam_path);
    if (gam_file.fail())
    {
        std::cout << "Failed to load GAM output" << std::endl;
        return (-1);
    }

    /* Read file and save paths */
    auto read_path = [&](vg::Alignment& alignment) {
        gam_out.push_back(alignment.path());
    };
    vg::io::for_each_parallel<vg::Alignment>(gam_file, read_path);
    /* Keep the best alignment */
    seq_path = gam_out.front();

    /* Prep distance index to read */
    std::string dist_path(argv[3]);
    bdsg::SnarlDistanceIndex dist_index;
    dist_index.deserialize(dist_path);

    /* Apply vg 'filters' */
    for (int i = 0; i < primers->primer_pairs.size(); i++)
        filter_primers(&seq_path, &(primers->primer_pairs[i]), &dist_index);

    /* Print in a meh format */
    print_vgprimers(primers);

    /* Leak free (I hope) */
    delete primers;
    primers = nullptr;

    return 0;
}

/* Pretty straight forward */
void print_help()
{
    std::cout << "FORMAT: vg_primers <primer_list.txt> <alignment.gam> <distance_index.dist>" << std::endl;
}

/* Print results */
void print_vgprimers(vgprimers *primers)
{
    std::cout << primers->name << std::endl;
    std::cout << primers->sequence << std::endl;
    for (int i = 0; i < primers->primer_pairs.size(); ++i)
    {
        std::cout << "=" << std::endl;
        primerpair pair = primers->primer_pairs[i];
        /* Print primers */
        std::cout << "Left primer " << pair.left.position << "," << pair.left.length << ": " << pair.left.sequence << std::endl;
        std::cout << "Right primer " << pair.right.position << "," << pair.right.length << ": " << pair.right.sequence << std::endl;
        /* Print snarls */
        std::cout << "Snarled primers:";
        if (pair.left.snarled) std::cout << " L";
        if (pair.right.snarled) std::cout << " R";
        if (!(pair.right.snarled || pair.left.snarled)) std::cout << " none";
        std::cout << std::endl;
        /* Print distances */
        std::cout << "Est max distance: " << pair.max_dist << std::endl;
        std::cout << "Est min distance: " << pair.min_dist << std::endl;
    }
}