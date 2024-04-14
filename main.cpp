#include "argument_parser.h"
#include "preprocess.h"
#include "snv.h"
#include "sv.h"

#include <cstring>
#include <iostream>
#include <map>
#include <fstream>
#include <stdlib.h>

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

        std::cout << "[preprocess] filtering p_utg unitigs\n";

        if (!preprocess::filter_unitigs(input, false)){
            return 1;
        }

        std::cout << "[preprocess] filtering out specified genomic regions\n";

        if(!preprocess::filter_regions(input)){
            return 1;
        }
    } else if (input.args["command"] == "snv"){
        if (!snv::check_args(input)){
            return 1;
        }
        if (!snv::run_pileup(input)){
            return 1;
        }
        if (!snv::call_snvs(input)){
            return 1;
        }
    } else if (input.args["command"] == "sv"){
        if (!sv::check_args(input)){
            return 1;
        }

        if (!sv::call_cis_chrom_svs(input)){
            return 1;
        }
    } else {
        std::cout << "Undefined command\n";
    }
    std::ofstream cmd_file(input.args["-o"] + "/command.txt");

    std::cout << "*************************\n";
    std::map<std::string, std::string> :: iterator it;
    for(it=input.args.begin();it !=input.args.end();++it){
        std::cout << it->first << ' ' <<it->second << '\n';
        cmd_file << it->first << ' ' <<it->second << '\n';
    }
    cmd_file.close();
    return 0;
}
