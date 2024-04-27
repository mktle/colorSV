#include "argument_parser.h"
#include "topology_search.h"

#include <assert.h>
#include <fstream>
#include <iostream>
#include <limits>
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
        user_args.args.insert({"--index_bin_size", "10000"});
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
