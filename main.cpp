#include "argument_parser.h"

#include <cstring>
#include <iostream>
#include <stdlib.h>

int main(int argc, char* argv[]){
    ArgumentParser input(argc, argv);
    if (input.args["command"] == "--help"){
        print_help();
        return 0;
    }
    if (input.args.count("-o") == 0){
        std::cout << "[ERROR] must specify output directory with -o\n";
        return 1;
    } else if (input.args["command"] == "preprocess"){
        std::cout << "Performing preprocessing!\n";
        std::cout << "Output path: " << input.args["-o"] << '\n';
    } else if (input.args["command"] == "snv"){
        std::cout << "Performing SNV calling!\n";
    } else if (input.args["command"] == "translocation"){
        std::cout << "Performing translocation calling!\n";
    } else {
        std::cout << "Undefined command\n";
    }
    //echo(argc, argv);
    //system("pwd");
    return 0;
}
