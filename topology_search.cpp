#include "argument_parser.h"
#include "topology_search.h"

#include <assert.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>

/* Checks that the user input all required flags */
bool topology_search::check_args(ArgumentParser& user_args){
    // check required flags
    std::list<std::string> required {"-o", "--r_graph"};
    if (!user_args.check_required_flags(required)){
        return false;
    }

    if (user_args.args.count("-t") == 0){
        user_args.args.insert({"-t", "3"});
    }

    if (user_args.args.count("-k") == 0){
        user_args.args.insert({"-k", "10"});
    }

    if (user_args.args.count("--index_bin_size") == 0){
        user_args.args.insert({"--index_bin_size", "100"});
    }
    return true;
}

// Index link file so the search is a bit faster
bool topology_search::index_link_file(ArgumentParser& user_args, std::unordered_map<int, std::streampos>& index_table){

    std::ifstream link_file(user_args.args["--r_graph"]);
    int bin_size {std::stoi(user_args.args["--index_bin_size"])};

    char line_type;
    std::string line_info;

    // jump to the area of the file with the 
    bool reached_links {false};
    while(!reached_links){
        char next = link_file.peek();
        if (next != 'L'){
            link_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }else{
            reached_links = true;
        }
    }

    // record the file position of every kth unitig, where k is the --index-bin_size argument
    while (link_file >> line_type){
        link_file >> line_info;
        int id {topology_search::utg_to_int(line_info)};

        // check if this the kth unitig; if so, record the location
        if (id % bin_size == 1){
            index_table[id] = link_file.tellg();

            link_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // skip past the lines with the same unitig so we don't record the same unitig again
            while (link_file >> line_type){
                link_file >> line_info;
                int next_id {topology_search::utg_to_int(line_info)};
                if (next_id == id){
                    // still on the same unitig, so skip to next line
                    link_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }else{
                    // on a new unitig, so can go back to original loop of checking unitig IDs
                    break;
                }
            }
        }
        // this will technically mean we skip over the next line after a recorded unitig
        // in practice, shouldn't matter because the recorded unitigs should not be consecutive
        link_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    link_file.close();
    return true;
}

// converts unitig ID to integer representation
// e.g., utg000016l becomes 16
int topology_search::utg_to_int(std::string& utg_id){
    // trim the "utg" at the beginning and the "l" at the end
    return std::stoi(utg_id.substr(3, utg_id.length() - 4));
}

bool topology_search::get_split_alignments(ArgumentParser& user_args, std::unordered_set<std::string>& candidates){
    // TODO: refactor alignment type parsing
    std::ifstream in_file(user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs_mapq_filtered.paf");
    std::ofstream out_file(user_args.args["-o"] + "/intermediate_output/sv_calls_utg_ids.paf");

    std::string prev_unitig;
    std::string unitig_id;
    std::string line_info;

    int p_count {0};

    // get the first unitig id and alignment type
    in_file >> prev_unitig;

    // skip to the tp tag
    // since alignment is run internally, in_file format will always be the same
    for (int i{0}; i < 16; i++){
        in_file >> line_info;
    }

    // double check we're on the right tag
    if (line_info.substr(0, 2) != "tp"){
        std::cout << "[topology_search::get_split_alignments][ERROR] unexpected formatting when parsing tumor-only unitig alignment file\n";
        return false;
    }
    if (line_info.back() == 'P'){
        p_count += 1;
    }
    in_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // iterate over lines in the paf file
    while(in_file >> unitig_id){
        if (unitig_id != prev_unitig){
            // new block of unitigs, so check if previous unitig had a split alignment
            if (p_count > 1){
                candidates.insert(prev_unitig);
            }

            // reset for next utg block
            prev_unitig = unitig_id;
            p_count = 0;
        }
        // check alignment type
        for (int i{0}; i < 16; i++){
            in_file >> line_info;
        }

        // double check we're on the right tag
        if (line_info.substr(0, 2) != "tp"){
            std::cout << "[topology_search::get_split_alignments][ERROR] unexpected formatting when parsing tumor-only unitig alignment file\n";
            return false;
        }
        if (line_info.back() == 'P'){
            p_count += 1;
        }
        in_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return true;
}

bool topology_search::run_topology_search(ArgumentParser& user_args, std::unordered_map<int, std::streampos>& index_table, std::unordered_set<std::string>& candidates, std::unordered_set<std::string>& result){

    int max_steps {std::stoi(user_args.args["-k"])};
    int bin_size {std::stoi(user_args.args["--index_bin_size"])};
    std::ifstream link_file(user_args.args["--r_graph"]);

    // run topology search on every candidate by iterating through the set
    std::unordered_set<std::string>::iterator itr;
    int cand_idx = 0;
    for (itr = candidates.begin(); itr != candidates.end(); itr++){

        std::string target_utg {*itr};

        // track the target unitig's neighbors, since we want to check if they can reach each other without the target
        std::unordered_set<std::string> target_neighbors;
        if (!get_neighbors(target_utg, link_file, target_neighbors, bin_size, index_table)){
            // this unitig is not in link file, so we will mark as false positive 
            continue;
        }

        // mark all candidates as "seen" because we only want to see if paths between neighbors exist without any candidate unitigs
        std::unordered_set<std::string> seen_nodes;
        seen_nodes.insert(candidates.begin(), candidates.end());

        // run BFS for k steps from one neighbor and see if we can get to all of the other neighbors
        std::string start_node {*target_neighbors.begin()};
        target_neighbors.erase(start_node);
        seen_nodes.insert(start_node);

        std::queue<std::string> to_traverse;
        // start traveling at a random neighbor
        to_traverse.push(start_node);
        // use '*' to separate search layers (so we can keep track of distance)
        to_traverse.push("*");

        int steps_taken {0};
        int num_neighbors_seen {0};
        bool to_remove {false};
        while (steps_taken <= max_steps){
            // get next node to explore
            std::string node {to_traverse.front()};
            to_traverse.pop();

            if (node == "*"){
                // we've reached the end of one search layer
                steps_taken++;
                to_traverse.push("*");
            }else{
                // check if this node is one of the target's neighbors
                if (target_neighbors.count(node)){
                    num_neighbors_seen += 1;
                    if (num_neighbors_seen == (int)target_neighbors.size()){
                        // successfully found a local path without candidate utgs
                        // so mark as a false positive
                        to_remove = true;
                        break;
                    }
                }

                // get current node's neighbors
                std::unordered_set<std::string> curr_neighbors;
                if (!get_neighbors(node, link_file, curr_neighbors, bin_size, index_table)){
                    std::cout << "[topology_search::run_topology_search][ERROR] error parsing link file when getting neighbors for " << node << '\n';
                    return false;
                }
            
                // add current node's neighbors to queue if they haven't already been explored
                std::unordered_set<std::string>::iterator itr2;
                for (itr2 = curr_neighbors.begin(); itr2 != curr_neighbors.end(); itr2++){
                    if (!seen_nodes.count(*itr2)){
                        to_traverse.push(*itr2);
                        seen_nodes.insert(*itr);
                    }
                }
            }
        }
        if (!to_remove){
            result.insert(target_utg);
            std::cout << target_utg << '\n';
        }
        cand_idx += 1;
        if (cand_idx % 50 == 0){
            std::cout << "[topology_search::run_topology_search] " << cand_idx << '/' << candidates.size() << " candidates checked\n";
        }
    }
    return true;
}

bool topology_search::get_neighbors(std::string& target_utg, std::ifstream& link_file, std::unordered_set<std::string>& neighbor_list, int bin_size, std::unordered_map<int, std::streampos>& index_table){
    std::string link_info;
    int utg_int_id {utg_to_int(target_utg)};

    // round current utg ID down to the beginning of the current bin to get lookup index
    int link_index {(utg_int_id/bin_size)*bin_size + 1};

    // jump to unitig location in link file
    link_file.seekg(index_table[link_index], std::ios::beg);

    // special case for the first unitig in a bin because of how the file was indexed
    // link table will point to position AFTER the first unitig name
    // i.e., if your first line in the bin is:
    //           L  utg1    +   utg901 ........
    // then the first character read from the line will be '+'
    // so this statement checks if our current unitig is the beginning of the bin
    if (utg_int_id % bin_size == 1){
        link_file >> link_info;
        link_file >> link_info;
        neighbor_list.insert(link_info);
    }
    link_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // sometimes unitigs are not reported in link file
    // need to track this, since candidate unitigs not in link files should be marked as a false positive
    bool found_utg {false};
    while(link_file >> link_info){
        // next character is unitig ID
        link_file >> link_info;
        if (found_utg && link_info != target_utg){
            // if we've already seen the target unitig, but the next unitig we're checking doesn't match, then that means we've added all neighbors
            break;
        }
        if (link_info == target_utg){
            found_utg = true;
            // found the right unitig, so add this neighbor
            link_file >> link_info;
            link_file >> link_info;
            neighbor_list.insert(link_info);
        }
        link_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // reset failbit to prevent bug after reaching EOF
    if (link_file.fail()){
        auto state = std::cin.rdstate();
        state &= ~std::ios_base::failbit;
        link_file.clear(state);
    }
    return found_utg;
}
