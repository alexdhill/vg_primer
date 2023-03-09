#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>
#include <bdsg/snarl_distance_index.hpp>
#include <handlegraph/snarl_decomposition.hpp>
#include <vg/vg.pb.h>

#include "vgfilters.hpp"
#include "libprimervg.hpp"


void filter_primers(vg::Path* seq_alignment, primerpair* primer_pair, bdsg::SnarlDistanceIndex *dist_index)
{
/*
	Find representative nodes for each primer using the sequence path in the graph.
	For each node, determine snarl involvement, and calculate the distances between
	them.
*/
	/*
		Find the positions correlating to the left and right primers
	*/

	/* Def node ranges to search for snarls */
	unsigned int left_node1;
	unsigned int left_node2;
	unsigned int right_node1;
	unsigned int right_node2;

	/* Prefix and suffix needed as the distance calculations */
	int prefix = 0; /* On the node that overlaps the left side of the left primer, the number of bases before the primer */
	int suffix = 0; /* On the node that overlaps the right side of the right primer, the number of bases after the primer */

	/* Iterators */
	unsigned int pos = 0; /* Current pos in sequence */
	unsigned int curr_node = 0; /* Current node in path */

	/* Find node overlapping start of left primer and find overlap size (prefix) */
	while (pos < (primer_pair->left.position) && (curr_node < seq_alignment->mapping_size()))
	{
		/* Find num bases aligned to this node */
		auto edits = seq_alignment->mapping(curr_node).edit();
		auto n_edits = edits.size();
		unsigned int matched_bases = 0;
		for (int i = 0; i < n_edits; i++)
		{
			matched_bases += edits[i].from_length();
		}
		/* Move forward */
		pos += matched_bases;
		++curr_node;
	}
	/* Reached end of sequence */
	if (curr_node >= seq_alignment->mapping_size())
	{
		std::cout << "Looking for " << primer_pair->left.position << std::endl;
		std::cout << pos << " in " << curr_node << "/" << seq_alignment->mapping_size() << std::endl;
		std::cout << "Ran off end of mapping trying to locate left primer" << std::endl;
		exit(-2);
	}
	/* Calculate prefix */
	for (int i = 0; i < seq_alignment->mapping(curr_node-1).edit_size(); ++i)
	{
		prefix += seq_alignment->mapping(curr_node-1).edit(i).from_length();
	}
	prefix -= (pos - primer_pair->left.position);
	left_node1 = curr_node-1;

	/* Find right node of left primer */
	while (pos < (primer_pair->left.position + primer_pair->left.length) && (curr_node < seq_alignment->mapping_size()))
	{
		/* Find num bases aligned to this node across edits */
		auto edits = seq_alignment->mapping(curr_node).edit();
		auto n_edits = edits.size();
		unsigned int matched_bases = 0;
		for (int i = 0; i < n_edits; i++)
		{
			matched_bases += edits[i].from_length();
		}
		pos += matched_bases;
		++curr_node;
	}
	left_node2 = curr_node-1;

	if (pos < (primer_pair->left.position + primer_pair->left.length)) {
		std::cout << "Ran off end of mapping trying to locate left primer" << std::endl;
		exit(-2);
	}

	/* Find node overlapping end of R primer and find overlap size (suffix) */
	while (pos < (primer_pair->right.position-primer_pair->right.length) && (curr_node < seq_alignment->mapping_size()))
	{
		/* Find num bases aligned to this node across edits */
		auto edits = seq_alignment->mapping(curr_node).edit();
		auto n_edits = edits.size();
		unsigned long int matched_bases = 0;
		for (int i = 0; i < n_edits; i++)
		{
			matched_bases += edits[i].from_length();
		}
		pos += matched_bases;
		++curr_node;
	}
	right_node1 = curr_node-1;

	if (pos < (primer_pair->right.position - primer_pair->right.length)) {
		std::cout << "Ran off end of mapping trying to locate right primer" << std::endl;
		exit(-2);
	}

	/* Find node overlapping end of R primer and find overlap size (suffix) */
	while (pos < primer_pair->right.position && (curr_node < seq_alignment->mapping_size()))
	{
		/* Find num bases aligned to this node across edits */
		auto edits = seq_alignment->mapping(curr_node).edit();
		auto n_edits = edits.size();
		unsigned long int matched_bases = 0;
		for (int i = 0; i < n_edits; i++)
		{
			matched_bases += edits[i].from_length();
		}
		pos += matched_bases;
		++curr_node;
	}

	if (pos < primer_pair->right.position) {
		std::cout << "Ran off end of mapping trying to locate right primer" << std::endl;
		exit(-2);
	}
	right_node2 = curr_node-1;

	for (int i = 0; i < seq_alignment->mapping(curr_node-1).edit_size(); i++)
	{
		suffix += seq_alignment->mapping(curr_node-1).edit(i).from_length();
	}
	suffix -= (pos - primer_pair->right.position);

	/*
		Find distances between nodes in distance index.
	*/
	vg::Position left_node = seq_alignment->mapping(left_node1).position();
	vg::Position right_node = seq_alignment->mapping(right_node2).position();

	primer_pair->min_dist = dist_index->minimum_distance(
		left_node.node_id(), false, left_node.offset()+prefix-1,
		right_node.node_id(), false, right_node.offset()+suffix, true
	);

	primer_pair->max_dist = dist_index->maximum_distance(
		left_node.node_id(), false, left_node.offset()+prefix-1,
		right_node.node_id(), false, right_node.offset()+suffix, true
	);

	/*
		Snarl gathering
	*/
	/* Check left primer for snarl */
	primer_pair->left.snarled = true;
	for (int i = left_node1; i <= left_node2; ++i)
	{
		auto curr = dist_index->get_node_net_handle(seq_alignment->mapping(i).position().node_id());
		auto curr_parent = dist_index->get_parent(curr);
		if ((dist_index->get_handle_type(curr_parent) == 0) ||
			(dist_index->get_handle_type(dist_index->get_parent(curr_parent)) == 0))
		{
			primer_pair->left.snarled = false;
			break;
		}
	}

	/* Check right primer for snarl */
	primer_pair->right.snarled = true;
	for (int i = right_node1; i <= right_node2; ++i)
	{
		auto curr = dist_index->get_node_net_handle(seq_alignment->mapping(i).position().node_id());
		auto curr_parent = dist_index->get_parent(curr);
		if ((dist_index->get_handle_type(curr_parent)) ||
			(dist_index->get_handle_type(dist_index->get_parent(curr_parent))))
		{
			primer_pair->right.snarled = false;
			break;
		}
	}
}


	/* Attempt 2 */

	// /* Def node ranges to search for snarls */
	// unsigned int left_node1 = -1;
	// unsigned int left_node2 = -1;
	// unsigned int right_node1 = -1;
	// unsigned int right_node2 = -1;

	// /* Prefix and suffix needed as the distance calculations */
	// unsigned int prefix = 0; /* On the node that overlaps the left side of the left primer, the number of bases before the primer */
	// unsigned int suffix = 0; /* On the node that overlaps the right side of the right primer, the number of bases after the primer */

	// /* Iterate through sequence until the ends of all primers are found */
	// unsigned int curr_pos = 0;
	// unsigned int curr_node = 0;
	// while (curr_node < seq_alignment->mapping_size())
	// {
	// 	auto curr_mapping = seq_alignment->mapping(curr_node);
	// 	unsigned int curr_pos_node = 0;
	// 	unsigned int curr_edit = 0;
	// 	while (curr_edit < curr_mapping.edit_size())
	// 	{
	// 		unsigned int curr_edit_size = curr_mapping.edit(curr_edit).from_length();
	// 		unsigned int curr_pos_edit = 0;
	// 		while (curr_pos_edit < curr_edit_size)
	// 		{
	// 			if (curr_pos == primer_pair->left.position)
	// 			{
	// 				left_node1 = curr_node;
	// 				prefix = curr_pos_edit;
	// 			}
	// 			else if (curr_pos == (primer_pair->left.position + primer_pair->left.length))
	// 			{
	// 				left_node2 = curr_node;
	// 			}
	// 			else if (curr_pos == (primer_pair->right.position - primer_pair->right.length))
	// 			{
	// 				right_node1 = curr_node;
	// 			} else if (curr_pos == primer_pair->right.position)
	// 			{
	// 				right_node2 = curr_node;
	// 				suffix = curr_pos_node;
	// 				break;
	// 			}
	// 			/* Increment my many counters */
	// 			++curr_pos_edit;
	// 			++curr_pos_node;
	// 			++curr_pos;
	// 		}
	// 		if (right_node2 != -1) break;
	// 		++curr_edit;
	// 	}
	// 	if (right_node2 != -1) break;
	// 	++curr_node;
	// }

	// /* Calculate distances */
	// vg::Position left_node = seq_alignment->mapping(left_node1).position();
	// vg::Position right_node = seq_alignment->mapping(right_node2).position();

	// primer_pair->min_dist = dist_index->minimum_distance(
	// 	left_node.node_id(), false, left_node.offset()+prefix-1,
	// 	right_node.node_id(), false, right_node.offset()+suffix, true
	// );

	// primer_pair->max_dist = dist_index->maximum_distance(
	// 	left_node.node_id(), false, left_node.offset()+prefix-1,
	// 	right_node.node_id(), false, right_node.offset()+suffix, true
	// );

	// /* Find snarls | ROOT_HANDLE=0, NODE_HANDLE, SNARL_HANDLE, CHAIN_HANDLE, SENTINEL_HANDLE */
	// primer_pair->left.snarled = false;
	// for (int i = left_node1; i <= left_node2; ++i)
	// {
	// 	/* Get node and parent */
	// 	auto curr = dist_index->get_node_net_handle(seq_alignment->mapping(i).position().node_id());
	// 	auto curr_parent = dist_index->get_parent(curr);
	// 	/* If parent or grandparent are root - it is not in a snarl */
	// 	if (!((dist_index->get_handle_type(curr_parent) == 0) ||
	// 		(dist_index->get_handle_type(dist_index->get_parent(curr_parent)) == 0)))
	// 	{
	// 		primer_pair->left.snarled = true;
	// 		break;
	// 	}
	// }

	// primer_pair->right.snarled = false;
	// for (int i = right_node1; i <= right_node2; ++i)
	// {
	// 	/* Get node and parent */
	// 	auto curr = dist_index->get_node_net_handle(seq_alignment->mapping(i).position().node_id());
	// 	auto curr_parent = dist_index->get_parent(curr);
	// 	/* If parent or grandparent are root - it is not a snarl */
	// 	if (!((dist_index->get_handle_type(curr_parent)) ||
	// 		(dist_index->get_handle_type(dist_index->get_parent(curr_parent)))))
	// 	{
	// 		primer_pair->right.snarled = true;
	// 		break;
	// 	}
	// }