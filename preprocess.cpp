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
    std::list<std::string> preprocess_required {"-o", "--graph", "--tumor-ids", "--reference", "-t"};
    if (!user_args.check_required_flags(preprocess_required)){
        std::cout << "[ERROR] preprocess command missing required flags. See sv-caller --help\n";
        return false;
    }

    // add default flag values if user did not specify values
    if (user_args.args.count("--min-reads") == 0){
        user_args.args.insert({"--min-reads", "2"});
    }

    if (user_args.args.count("-t") == 0){
        user_args.args.insert({"-t", "3"});
    }

    if (!user_args.check_file("--graph", ".gfa")){
        return false;
    }

    if (!user_args.check_file("--reference", ".fa")){
        return false;
    }

    return true;
}

/* Identifies tumor-only unitigs from given .gfa file */
bool preprocess::filter_unitigs(ArgumentParser& user_args){
    std::cout << "[PREPROCESS] classifying unitigs in .gfa file\n";
    // open input and output files
    std::ifstream gfa_file(user_args.args["--graph"]);
    std::ofstream out_tumor_utg(user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.txt");
    std::ofstream out_tumor_fa(user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.fa");

    std::ofstream out_healthy_utg(user_args.args["-o"] + "/intermediate_output/healthy_only_unitigs.txt");
    std::ofstream out_healthy_fa(user_args.args["-o"] + "/intermediate_output/healthy_only_unitigs.fa");

    std::ofstream out_mixed_utg(user_args.args["-o"] + "/intermediate_output/mixed_only_unitigs.txt");
    std::ofstream out_mixed_fa(user_args.args["-o"] + "/intermediate_output/mixed_only_unitigs.fa");

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

    int read_thresh{stoi(user_args.args["--min-reads"])};
    std::string line_type;
    std::string unitig; 
    std::string ignore;
    std::string segment;
    std::string read_id;

    // double check we're starting with a segment line
    gfa_file >> line_type;
    if (line_type.c_str()[0] != 'S'){
        std::cout << "[ERROR] .gfa file must begin with S line\n";
        return false;
    }
    gfa_file >> unitig >> segment;

    gfa_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    int num_reads{0};
    bool tumor{false};
    bool normal{false};

    while(gfa_file >> line_type){
        if (line_type.c_str()[0] == 'S'){
            // end of the node, so reset for the next node
            // first, write unitig info if this is a tumor-only node
            if (tumor && !normal && num_reads >= read_thresh){
                out_tumor_utg << unitig << '\n';
                out_tumor_fa << '>' << unitig << '\n' << segment << '\n';
            }else if (!tumor && normal){
                out_healthy_utg << unitig << '\n';
                out_healthy_fa << '>' << unitig << '\n' << segment << '\n';
            }else if (tumor && normal){
                out_mixed_utg << unitig << '\n';
                out_mixed_fa << '>' << unitig << '\n' << segment << '\n';
            }else{
                assert(tumor && num_reads < read_thresh);
            }
            gfa_file >> unitig >> segment;

            tumor = false;
            normal = false;
            num_reads = 0;
        }else if (line_type.c_str()[0] == 'A'){
            num_reads++;
            gfa_file >> ignore >> ignore >> ignore >> read_id;
            // check if this read comes from a non-tumor sample
            if (std::find(std::begin(tumor_ids), std::end(tumor_ids), read_id.substr(0, read_id.find('/'))) == std::end(tumor_ids)){
                normal = true;
            }else if (std::find(std::begin(tumor_ids), std::end(tumor_ids), read_id.substr(0, read_id.find('/'))) != std::end(tumor_ids)){
                tumor = true;
            }
        }
        gfa_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // write last node if it's tumor
    if (tumor && !normal && num_reads >= read_thresh){
        out_tumor_utg << unitig << '\n';
        out_tumor_fa << '>' << unitig << '\n' << segment << '\n';
    }

    // close input and output files
    gfa_file.close();
    out_tumor_utg.close();
    out_tumor_fa.close();

    out_healthy_utg.close();
    out_healthy_fa.close();

    out_mixed_utg.close();
    out_mixed_fa.close();

    return true;
}

bool preprocess::filter_regions(ArgumentParser& user_args){
    // alignment to reference
    std::cout << "[PREPROCESS] aligning unitigs to reference genome\n\n";
    std::string cmd{"minimap2 -ax map-hifi -s50 -t " + user_args.args["-t"] + " " + user_args.args["--reference"] + " " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs.fa -o " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs_aln.sam"};
    system(cmd.c_str());

    // filter unwanted regions
    if (user_args.args.count("--filter") > 0){
        std::cout << "\n[PREPROCESS] filtering out regions in " << user_args.args["--filter"] << "\n";
        cmd = "samtools view -h -L " + user_args.args["--filter"] + " -U " + user_args.args["-o"] + "/intermediate_output/filtered_unitigs_aln.sam -o /dev/null " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs_aln.sam";
    }else{
        cmd = "mv " + user_args.args["-o"] + "/intermediate_output/tumor_only_unitigs_aln.sam " + user_args.args["-o"] + "/intermediate_output/filtered_unitigs_aln.sam";
    }
    system(cmd.c_str());

    // sort final result
    cmd = "samtools view -bS " + user_args.args["-o"] + "/intermediate_output/filtered_unitigs_aln.sam | samtools sort -o " + user_args.args["-o"] + "/intermediate_output/filtered_unitigs_aln.srt.bam -@" + user_args.args["-t"] + "-m4g ";
    system(cmd.c_str());
    return true;
}
