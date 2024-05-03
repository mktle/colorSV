#ifndef TOPOLOGY_SEARCH_H
#define TOPOLOGY_SEARCH_H

#include <string>
#include <unordered_map>
#include <unordered_set>

#include "argument_parser.h"

namespace topology_search{
	bool check_args(ArgumentParser& user_args);
	bool index_link_file(ArgumentParser& user_args, std::unordered_map<int, std::streampos>& index_table);
	int utg_to_int(std::string& utg_id);
	bool get_split_alignments(ArgumentParser& user_args, std::unordered_set<std::string>& candidates);
	bool run_topology_search(ArgumentParser& user_args, std::unordered_map<int, std::streampos>& index_table, std::unordered_set<std::string>& candidates, std::unordered_set<std::string>& result);
	bool get_neighbors(std::string& target_utg, std::ifstream& link_file, std::unordered_set<std::string>& neighbor_list, int bin_size, std::unordered_map<int, std::streampos>& index_table);

	bool write_final_paf(ArgumentParser& args, std::unordered_set<std::string>& sv_set);
}

#endif

