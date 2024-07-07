#include "argument_parser.h"

#include <cstring>
#include <list>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>

ArgumentParser::ArgumentParser(int &argc, char** argv){
    if(argc == 1){
        this->args.insert({"command", "--help"});
    }
    // first argument should indicate valid command; otherwise throw error
    else if(std::strcmp(*(argv + 1), "preprocess") && std::strcmp(*(argv + 1), "call") && std::strcmp(*(argv + 1), "--help") && std::strcmp(*(argv + 1), "sv")){
        throw std::invalid_argument("Command not found, see colorSV --help for valid commands");
    }else {
        std::string executable {*(argv)};
        this->args.insert({"exe_path", executable.substr(0, executable.find_last_of('/') + 1)});

        // add the rest of the command-line options to map of args
        this->args.insert({"command", *(argv+1)});
        int i{2};
        int prev_opt_index{-1};
        while (*(argv + i) != 0){
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
                        opts += ",";
                        opts += *(argv + prev_opt_index + k);
                    }
                    this->args.insert({*(argv + prev_opt_index), opts});
                }
                prev_opt_index = i;
            }
            i++;
        }
        if (argc > 2){
            // push final flag
            // check whether they specified arguments after flag
            if (prev_opt_index == i - 1){
                throw std::invalid_argument("Must specify option after flags");
            }
            std::string opts{*(argv + prev_opt_index + 1)};
            for (int k{2}; prev_opt_index + k < i; k++){
                opts += ",";
                opts += *(argv + prev_opt_index + k);
            }
            this->args.insert({*(argv + prev_opt_index), opts});
        }
    }
}

bool ArgumentParser::check_required_flags(std::list<std::string>& r_flags){
    for (auto it = r_flags.begin(); it != r_flags.end(); it++){
        if (this->args.count(*it) == 0){
            std::cout << "[ArgumentParser::check_required_flags][ERROR] missing required flag " << *it << '\n';
            // required flag not in args list
            return false;
        }
    }
    return true;
}

bool ArgumentParser::check_file(std::string opt, std::string ext){
    struct stat buffer;   
    if (stat(this->args[opt].c_str(), &buffer) != 0){
        std::cout << "[ArgumentParser::check_file][ERROR] file specified in " << opt << " is not readable\n";
        return false;
    }

    size_t del_pos {this->args[opt].find_last_of('.')};
    if (del_pos == std::string::npos || this->args[opt].substr(del_pos) != ext){
        std::cout << "[ArgumentParser::check_file][ERROR] file specified in " << opt << " must be type " << ext << '\n';
        return false;
    }

    return true;
}

void print_help(){
    std::cout << "Usage: colorSV <command> -o <output path> [options]\n";
    std::cout << "Commands and options:\n";
    std::cout << "  * preprocess\n";
    std::cout << "     <required flags>\n";
    std::cout << "          --graph             STR     path to assembly graph file\n";
    std::cout << "          --reference         STR     path to reference genome file\n";
    std::cout << "          --tumor-ids         STR     tumor sample identifiers\n";
    std::cout << "          --read-sep          STR     delimiter in read names (e.g., / or .)\n";
    std::cout << "     [optional flags] \n";
    std::cout << "          --min-reads         INT     minimum number of reads when identifying tumor-only unitigs [2]\n";
    std::cout << "          --min-mapq          INT     minimum MAPQ required for tumor-only unitig alignments [10]\n";
    std::cout << "          -t                  INT     number of threads during alignment [3]\n";
    std::cout << "  * call\n";
    std::cout << "     <required flags>\n";
    std::cout << "          --graph             STR     path to assembly graph file\n";
    std::cout << "          --filter            STR     path to BED file with regions to ignore (e.g., centromeres)\n";
    std::cout << "     [optional flags] \n";
    std::cout << "          -k                  INT     maximum number of steps in topology search [10]\n";
    std::cout << "          -q                  INT     minimum MAPQ of alignments when extracting breakpoints [15]\n";
    std::cout << "          -Q                  INT     minimum MAPQ for alignment ends when extracting breakpoints [15]\n";
    std::cout << "          --index-bin-size    INT     unitigs per bin when indexing assembly graph file [100]\n";
}
