#include "argument_parser.h"

#include <cstring>
#include <iostream>
#include <stdexcept>

ArgumentParser::ArgumentParser(int &argc, char** argv){
    if(argc == 1){
        this->args.insert({"command", "--help"});
    }
    // first argument should indicate valid command; otherwise throw error
    else if(std::strcmp(*(argv + 1), "preprocess") and std::strcmp(*(argv + 1), "snv") and std::strcmp(*(argv + 1), "translocation") and std::strcmp(*(argv + 1), "--help")){
        throw std::invalid_argument("Command not found, see sv-caller --help for valid commands");
    }else {
        // add the rest of the command-line options to map of args
        this->args.insert({"command", *(argv+1)});
        int i{2};
        int prev_opt_index{-1};
        while (*(argv + i) != '\0'){
            // check if this is an option by looking for '-' prefix
            if (*(argv + i)[0] == '-'){
                // if this is not the first option, add previous option to map
                if (prev_opt_index != -1){
                    // check whether they specified arguments after flag
                    if (prev_opt_index == i - 1){
                        throw std::invalid_argument("Must specify option after flags");
                    }
                    // if valid, add options to argument map
                    std::string opts{*(argv + prev_opt_index + 1)};
                    for (int k{2}; prev_opt_index + k < i; k++){
                        opts += " ";
                        opts += *(argv + prev_opt_index + k);
                    }
                    this->args.insert({*(argv + prev_opt_index), opts});
                }
                prev_opt_index = i;
            }
            i++;
        }
        // push final option
        // check whether they specified arguments after flag
        if (prev_opt_index == i - 1){
            throw std::invalid_argument("Must specify option after flags");
        }
        // if valid, add options to argument map
        std::string opts{*(argv + prev_opt_index + 1)};
        for (int k{2}; prev_opt_index + k < i; k++){
            opts += " ";
            opts += *(argv + prev_opt_index + k);
        }
        this->args.insert({*(argv + prev_opt_index), opts});
    }
}

void print_help(){
    std::cout << "Usage: sv-caller <mutation type> <-o OUTPUT PATH> [options]\n";
    std::cout << "Commands and options:\n";
    std::cout << "  preprocess\n";
    std::cout << "     --graph      STR     path to assembly graph file\n";
    std::cout << "     --normal-ids STR     normal sample identifiers\n";
    std::cout << "     --tumor-ids  STR     tumor sample identifiers\n";
    std::cout << "  snv\n";
    std::cout << "     -f           FLOAT   minimum fraction of tumor reads supporting alt allele (default 0.1)\n";
    std::cout << "     -m           INT     minimum number of supporting tumor reads (default 10)\n";
    std::cout << "  translocation\n";
    std::cout << "     -k           INT     maximum number of steps in topology search\n";
}
