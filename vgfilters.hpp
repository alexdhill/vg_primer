#ifndef VGFILTERS_H
#define VGFILTERS_H

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>
#include <bdsg/snarl_distance_index.hpp>
#include <handlegraph/snarl_decomposition.hpp>
#include <vg/vg.pb.h>

#include "libprimervg.hpp"

void filter_primers(vg::Path* seq_alignment, primerpair* primer_pair, bdsg::SnarlDistanceIndex *dist_index);

#endif