#ifndef LIBPRIMERVG_H
#define LIBPRIMERVG_H

#include <string>
#include <vector>

typedef struct primer
{
    /* The primer sequence */
    std::string sequence;
    /* Position of the primer relative to the given sequence strand */
    unsigned int position;
	/* Length of the primer */
	int length;
    /* Primer is in a snarl */
    bool snarled;

} primer;

typedef struct primerpair
{
    
    primer left;
    primer right;
    unsigned int min_dist;
    unsigned int max_dist;
} primerpair;

typedef struct vgprimers
{
    std::string name;
    std::string sequence;
    std::vector<primerpair> primer_pairs;
} vgprimers;

#endif