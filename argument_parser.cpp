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
        for(int i = 2; argc > 3 && *(argv + i) != '\0'; i += 2){
            this->args.insert({*(argv+i), *(argv + i + 1)});
        }
    }
}

void print_help(){
    std::cout << "Usage: sv-caller <mutation type> <-o OUTPUT PATH> [options]\n";
    std::cout << "Commands and options:\n";
    std::cout << "  -- preprocess\n";
    std::cout << "  -- snv\n";
    std::cout << "     -f FLOAT     minimum fraction of tumor reads supporting alt allele (default 0.1)\n";
    std::cout << "     -m INT       minimum number of supporting tumor reads (default 10)\n";
    std::cout << "  -- translocation\n";
    std::cout << "     -k INT       maximum number of steps in topology search\n";
}
