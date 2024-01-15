#include "argument_parser.h"
#include "preprocess.h"
#include "snv.h"

#include <cstring>
#include <iostream>
#include <map>
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

        if (!preprocess::filter_unitigs(input)){
            return 1;
        }

        if(!preprocess::filter_regions(input)){
            return 1;
        }
    } else if (input.args["command"] == "snv"){
        if (!snv::check_args(input)){
            return 1;
        }
        std::cout << "Performing SNV calling!\n";
    } else if (input.args["command"] == "translocation"){
        std::cout << "Performing translocation calling!\n";
    } else {
        std::cout << "Undefined command\n";
    }
    std::cout << "*************************\n";
    std::map<std::string, std::string> :: iterator it;
    for(it=input.args.begin();it !=input.args.end();++it){
        std::cout << it->first << ' ' <<it->second << '\n';
    }
    return 0;
}
