#include "argument_parser.h"

#include <boost/filesystem.hpp>

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
            // required flag not in args list
            return false;
        }
    }
    return true;
}

bool ArgumentParser::check_file(std::string opt, std::string ext){
    struct stat buffer;   
    if (stat(this->args[opt].c_str(), &buffer) != 0){
        std::cout << "[ERROR] file specified in " << opt << " is not readable\n";
        return false;
    }

    if (boost::filesystem::extension(this->args[opt]) != ext){
        std::cout << "[ERROR] file specified in " << opt << " must be type " << ext << '\n';
        return false;
    }
    return true;
}

void print_help(){
    std::cout << "Usage: sv-caller <mutation type> <-o OUTPUT PATH> [other flags]\n";
    std::cout << "Commands and options:\n";
    std::cout << "  * preprocess\n";
    std::cout << "     <required flags>\n";
    std::cout << "          --r_graph           STR     path to r_utg assembly graph file\n";
    std::cout << "          --u_graph           STR     path to u_utg assembly graph file\n";
    std::cout << "          --reference     STR     path to reference genome file\n";
    std::cout << "          -t              INT     number of threads\n";
    std::cout << "          --tumor-ids     STR     tumor sample identifiers\n";
    std::cout << "     [optional flags] \n";
    std::cout << "          --min-reads     INT     minimum number of reads when identifying tumor-only unitigs (default 2)\n";
    std::cout << "          --filter        STR     path to files with regions to filter (e.g., centromeres)\n";
    std::cout << "  * snv\n";
    std::cout << "     [optional flags] \n";
    std::cout << "          --frac          FLOAT   minimum fraction of tumor reads supporting alt allele (default 0.1)\n";
    std::cout << "          --normal-reads  STR     path to .bam file with normal reads\n";
    std::cout << "          --pileup        STR     path to existing pileup file\n";
    std::cout << "          --reference     STR     path to reference genome file\n";
    std::cout << "          --tumor-reads   STR     path to .bam file with tumor reads\n";
    std::cout << "          -m              INT     minimum number of supporting tumor reads (default 10)\n";
    std::cout << "  * translocation\n";
    std::cout << "     <required flags>\n";
    std::cout << "          --graph         STR     path to assembly graph file\n";
    std::cout << "     [optional flags] \n";
    std::cout << "          -k              INT     maximum number of steps in topology search (default 10)\n";
}
