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
    std::list<std::string> preprocess_required {"-o", "--graph", "--tumor-ids", "--reference", "-t", "--read-sep"};
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

    if (user_args.args.count("--min-mapq") == 0){
        user_args.args.insert({"--min-mapq", "10"});
    }

    if (!user_args.check_file("--graph", ".gfa")){
        std::cout << "[preprocess::check_args][ERROR] could not open graph file: " << user_args.args["--graph"] << '\n';
        return false;
    }

    if (!user_args.check_file("--reference", ".fa")){
        std::cout << "[preprocess::check_args][ERROR] could not open reference file: " << user_args.args["--reference"] << '\n';
        return false;
    }

    return true;
}

/* Identifies tumor-only unitigs from given .gfa file */
bool preprocess::filter_unitigs(ArgumentParser& user_args){
    std::ifstream gfa_file(user_args.args["--graph"]);
    std::ofstream out_tumor_utg_all(user_args.args["-o"] + "/intermediate_output/all_tumor_only_unitigs.txt");
    std::ofstream out_tumor_utg_thresh(user_args.args["-o"] + "/intermediate_output/thresh_tumor_only_unitigs.txt");
    std::ofstream out_tumor_fa(user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.fa");

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

    int read_thresh {std::stoi(user_args.args["--min-reads"])};
    char read_delim {user_args.args["--read-sep"][0]};

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
            if(num_tumor_reads >= 0 && num_healthy_reads == 0){
                out_tumor_utg_all << unitig << '\n';
                if (num_tumor_reads >= read_thresh){
                    out_tumor_utg_thresh << unitig << '\n';
                    out_tumor_fa << '>' << unitig << '\n' << segment << '\n';
                }
            }

            gfa_file >> unitig >> segment;

            num_tumor_reads = 0;
            num_healthy_reads = 0;
        }else if (line_type.c_str()[0] == 'A'){
            gfa_file >> ignore >> ignore >> ignore >> read_id;
            // check if this read comes from a non-tumor sample
            if (std::find(std::begin(tumor_ids), std::end(tumor_ids), read_id.substr(0, read_id.find(read_delim))) != std::end(tumor_ids)){
                num_tumor_reads++;
            }else{
                num_healthy_reads++;
            }
        }
        gfa_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if(num_tumor_reads >= 0 && num_healthy_reads == 0){
        out_tumor_utg_all << unitig << '\n';
        if (num_tumor_reads >= read_thresh){
            out_tumor_utg_thresh << unitig << '\n';
            out_tumor_fa << '>' << unitig << '\n' << segment << '\n';
        }
    }

    // close input and output files
    gfa_file.close();
    out_tumor_utg_all.close();
    out_tumor_utg_thresh.close();
    out_tumor_fa.close();

    return true;
}

bool preprocess::align_unitigs(ArgumentParser& user_args){
    // r_utg tumor-only unitig alignment to reference

    std::string cmd{"./minimap2 -cx lr:hq -t" + user_args.args["-t"] + " --ds " + user_args.args["--reference"] + " " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.fa > " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.paf"};
    system(cmd.c_str());

    // filter to only keep alignments with minimum MAPQ score
    cmd = "awk '$12 >= " + user_args.args["--min-mapq"] + "' " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.paf > " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs_mapq_filtered.paf";
    system(cmd.c_str());

    return true;
}
