#include "argument_parser.h"
#include "preprocess.h"
#include "topology_search.h"

#include <cstring>
#include <iostream>
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

        if (!topology_search::run_topology_search(input)){
            return 1;
        }
    }else{
        std::cout << "Undefined command\n";
    }
    std::ofstream cmd_file(input.args["-o"] + "/command.txt");

    std::cout << "*************************\n";
    std::map<std::string, std::string> :: iterator it;
    for (it=input.args.begin();it !=input.args.end();++it){
        std::cout << it->first << ' ' <<it->second << '\n';
        cmd_file << it->first << ' ' <<it->second << '\n';
    }
    cmd_file.close();
    return 0;
}
