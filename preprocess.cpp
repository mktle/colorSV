#include "argument_parser.h"
#include "preprocess.h"

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

bool preprocess::file_setup(ArgumentParser& user_args){
    // create output directory if it doesn't already exist
    std::string cmd {"mkdir -p " + user_args.args["-o"] + "/intermediate_output/"};
    system(cmd.c_str());

    return true;
}

/* Checks that the user input all required flags */
bool preprocess::check_args(ArgumentParser& user_args){
    // check required flags
    std::list<std::string> preprocess_required {"-o", "--r_graph", "--p_graph", "--tumor-ids", "--reference", "-t"};
    if (!user_args.check_required_flags(preprocess_required)){
        return false;
    }

    // add default flag values if user did not specify values
    if (user_args.args.count("--min-reads") == 0){
        user_args.args.insert({"--min-reads", "2"});
    }

    if (user_args.args.count("-t") == 0){
        user_args.args.insert({"-t", "3"});
    }

    if (!user_args.check_file("--r_graph", ".gfa")){
        return false;
    }

    if (!user_args.check_file("--p_graph", ".gfa")){
        return false;
    }

    if (!user_args.check_file("--reference", ".fa")){
        return false;
    }

    return true;
}

/* Identifies tumor-only unitigs from given .gfa file */
bool preprocess::filter_unitigs(ArgumentParser& user_args, bool is_r_utg){
    std::cout << "[preprocess::filter_unitigs] classifying unitigs in .gfa file\n";

    std::ifstream gfa_file;

    std::ofstream out_tumor_utg;
    std::ofstream out_tumor_fa;

    // open input and output files
    if(is_r_utg){
        gfa_file.open(user_args.args["--r_graph"]);
        out_tumor_utg.open(user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.txt");
        out_tumor_fa.open(user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.fa");
    }else{
        gfa_file.open(user_args.args["--p_graph"]);

        out_tumor_utg.open(user_args.args["-o"] + "/intermediate_output/p_utg_tumor_only_unitigs.txt");
        out_tumor_fa.open(user_args.args["-o"] + "/intermediate_output/p_utg_tumor_only_unitigs.fa");
    }

    // parse tumor sample IDs into list
    std::vector<std::string> tumor_ids;
    std::string s = user_args.args["--tumor-ids"];
    std::string::const_iterator start = s.begin();
    std::string::const_iterator end = s.end();
    std::string::const_iterator next = std::find(start, end, ',');
    while (next != end) {
        tumor_ids.push_back(std::string(start, next));
        start = next + 1;
        next = std::find(start, end, ',');
    }
    tumor_ids.push_back(std::string(start, next));

    int read_thresh{ std::stoi(user_args.args["--min-reads"]) };
    if (!is_r_utg){
        // since p_utg is being used to remove false positives, allow all tumor-only unitigs regardless of number of reads
        read_thresh = 1;
    }

    // variables for parsing .gfa file
    std::string line_type;
    std::string unitig; 
    std::string ignore;
    std::string segment;
    std::string read_id;
    int num_tumor_reads{0};
    int num_healthy_reads{0};

    // double check we're starting with a segment line
    gfa_file >> line_type;
    if (line_type.c_str()[0] != 'S'){
        std::cout << "[preprocess::filter_unitigs][ERROR] .gfa file must begin with S line\n";
        return false;
    }
    gfa_file >> unitig >> segment;

    gfa_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');


    while(gfa_file >> line_type){
        if (line_type.c_str()[0] == 'S'){
            // end of the node, so reset for the next node
            // first, write unitig info if this is a tumor-only node
            if(num_tumor_reads >= read_thresh && num_healthy_reads == 0){
                out_tumor_utg << unitig << ' ' << num_tumor_reads << '\n';
                out_tumor_fa << '>' << unitig << '\n' << segment << '\n';
            }

            gfa_file >> unitig >> segment;

            num_tumor_reads = 0;
            num_healthy_reads = 0;
        }else if (line_type.c_str()[0] == 'A'){
            gfa_file >> ignore >> ignore >> ignore >> read_id;
            // check if this read comes from a non-tumor sample
            if (std::find(std::begin(tumor_ids), std::end(tumor_ids), read_id.substr(0, read_id.find('/'))) != std::end(tumor_ids)){
                num_tumor_reads++;
            }else{
                num_healthy_reads++;
            }
        }
        gfa_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if(num_tumor_reads >= read_thresh && num_healthy_reads == 0){
        out_tumor_utg << unitig << ' ' << num_tumor_reads << '\n';
        out_tumor_fa << '>' << unitig << '\n' << segment << '\n';
    }

    // close input and output files
    gfa_file.close();
    out_tumor_utg.close();
    out_tumor_fa.close();

    return true;
}

bool preprocess::align_unitigs(ArgumentParser& user_args){
    // r_utg tumor-only unitig alignment to reference
    std::cout << "[preprocess::align_unitigs] aligning r_utg tumor unitigs to reference genome\n\n";
    std::string cmd{"./minimap2 -cx lr:hq -t" + user_args.args["-t"] + " --ds " + user_args.args["--reference"] + " " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.fa > " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.paf"};
    system(cmd.c_str());

    return true;
}
