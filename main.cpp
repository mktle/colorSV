#include "argument_parser.h"
#include "preprocess.h"
#include "topology_search.h"

#include <cstring>
#include <iostream>
#include <iterator>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <unordered_map>

int main(int argc, char* argv[]){
    ArgumentParser input(argc, argv);
    std::cout << '\n';
    if (input.args["command"] == "--help"){
        print_help();
        return 0;
    }

    if (input.args["command"] == "preprocess"){
        // check that the user input all required flags
        if (!preprocess::check_args(input)){
            return 1;
        }

        // perform file and directory setup
        if (!preprocess::file_setup(input)){
            return 1;
        }

        std::cout << "[preprocess] filtering r_utg unitigs\n";

        if (!preprocess::filter_unitigs(input, true)){
            return 1;
        }

        std::cout << "[preprocess] performing unitig alignment\n";

        if(!preprocess::align_unitigs(input)){
            return 1;
        }
    }else if (input.args["command"] == "call"){
        if (!topology_search::check_args(input)){
            return 1;
        }

        std::unordered_map<int, std::streampos> link_index;
        if (!topology_search::index_link_file(input, link_index)){
            return 1;
        }

        std::unordered_set<std::string> candidate_utgs;
        if (!topology_search::get_split_alignments(input, candidate_utgs)){
            return 1;
        }

        std::cout << "[call] number of candidate unitigs before topology search: " << candidate_utgs.size() << '\n';

        std::unordered_set<std::string> final_svs;
        if (!topology_search::run_topology_search(input, link_index, candidate_utgs, final_svs)){
            return 1;
        }

        std::cout << "[call] number of unitigs that pass all filters: " << final_svs.size() << '\n';

    }else{
        std::cout << "Undefined command\n";
    }
    std::ofstream cmd_file(input.args["-o"] + "/command.txt", std::ios_base::app);

    std::cout << "*************************\n";
    cmd_file << "*************************\n";
    std::map<std::string, std::string> :: iterator it;
    for (it=input.args.begin();it !=input.args.end();++it){
        std::cout << it->first << ' ' <<it->second << '\n';
        cmd_file << it->first << ' ' <<it->second << '\n';
    }
    cmd_file.close();
    return 0;
}
