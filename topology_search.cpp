#include "argument_parser.h"
#include "topology_search.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

/* Checks that the user input all required flags */
bool topology_search::check_args(ArgumentParser& user_args){
    // check required flags
    std::list<std::string> required {"-o", "--r_graph"};
    if (!user_args.check_required_flags(required)){
        return false;
    }

    if (user_args.args.count("-t") == 0){
        user_args.args.insert({"-t", "3"});
    }

    if (user_args.args.count("-k") == 0){
        user_args.args.insert({"-k", "10"});
    }

    if (user_args.args.count("--index_bin_size") == 0){
        user_args.args.insert({"--index_bin_size", "10000"});
    }
    return true;
}

// Index link file so the search is a bit faster
bool topology_search::index_link_file(ArgumentParser& user_args, std::unordered_map<int, std::streampos>& index_table){

    std::ifstream link_file(user_args.args["--r_graph"]);

    std::streampos stream_loc;
    char line_type;
    std::string line_info;

    // jump to the area of the file with the 
    bool reached_links {false};
    while(!reached_links){
        char next = link_file.peek();
        if (next != 'L'){
            link_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }else{
            reached_links = true;
        }
    }

    // record the start of the links location
    stream_loc = link_file.tellg();
    index_table[1] = stream_loc;
    link_file.close();

    std::ifstream test_stream(user_args.args["--r_graph"]);
    test_stream.seekg(index_table[1]);

    for (int i {0}; i < 5; i++){
        test_stream >> line_info;
        std::cout << line_info << ' ';
    }
    std::cout << '\n';
    link_file.close();

    return true;
}
