#include "argument_parser.h"

#include <cstring>
#include <iostream>
#include <stdexcept>

ArgumentParser::ArgumentParser(int &argc, char** argv){
    // print list of options when there's no argument
    if(argc == 1){
        std::cout << "Usage: sv-caller <mutation type> <-o OUTPUT PATH> [options]\n";
        std::cout << "Commands and options:\n";
        std::cout << "  -- preprocess\n";
        std::cout << "  -- snv\n";
        std::cout << "     -f FLOAT     minimum fraction of tumor reads supporting alt allele (default 0.1)\n";
        std::cout << "     -m INT       minimum number of supporting tumor reads (default 10)\n";
        std::cout << "  -- translocation\n";
        std::cout << "     -k INT       maximum number of steps in topology search\n";
        this->args.insert({"command", "help"});
    }
    // first argument should indicate snv or translocation; otherwise throw error
    else if(std::strcmp(*(argv + 1), "preprocess") and std::strcmp(*(argv + 1), "snv") and std::strcmp(*(argv + 1), "translocation")){
        throw std::invalid_argument("Command not found");
    }else {
        this->args.insert({"command", *(argv+1)});

        for(int i = 2; argc > 2 && *(argv + i) != '\0'; i += 2){
            this->args.insert({*(argv+i), *(argv + i + 1)});
        }
    }
}
    
void echo(int &argc, char **argv){
    std::cout << argc << '\n';
    std::cout << *argv << '\n';
    if (argc > 2){
        std::cout << *(argv + 4) << '\n';
        std::cout << *(argv + 5) << '\n';
    }
}

