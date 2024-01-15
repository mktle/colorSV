#include "argument_parser.h"
#include "snv.h"

#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <string>

/* Checks that the user input all required flags */
bool snv::check_args(ArgumentParser& user_args){
    // check required flags
    std::list<std::string> preprocess_required {"-o"};
    if (!user_args.check_required_flags(preprocess_required)){
        std::cout << "[ERROR] snv command missing required flags. See sv-caller --help\n";
        return false;
    }

    if (user_args.args.count("--pileup") == 0 && (user_args.args.count("--normal-reads") == 0 ||  user_args.args.count("--tumor-reads") == 0 || user_args.args.count("--reference") == 0)){
        std::cout << "[ERROR] snv command requires either existing pileup file, or paths to tumor/normal reads and reference genome\n";
        return false;
    }

    // add default flag values if user did not specify values
    if (user_args.args.count("-m") == 0){
        user_args.args.insert({"-m", "1"});
    }

    if (user_args.args.count("--frac") == 0){
        user_args.args.insert({"--frac", "0.1"});
    }

    // check if pileup file is readable
    if (user_args.args.count("--pileup") > 0 && !(user_args.check_file("--pileup", ".vcf"))){
        return false;
    }

    // check if read files are readable
    if (user_args.args.count("--normal-reads") > 0){
        if (!(user_args.check_file("--normal-reads", ".bam"))){
            return false;
        }

        assert(user_args.args.count("--tumor-reads") > 0 && user_args.args.count("--reference") > 0);

        if (!(user_args.check_file("--tumor-reads", ".bam"))){
            return false;
        }

        if (!(user_args.check_file("--reference", ".fa"))){
            return false;
        }
    }
    return true;
}

bool snv::run_pileup(ArgumentParser& user_args){
    // only run pileup if user has not specified an existing file
    if (user_args.args.count("--pileup") > 0){
        return true;
    }
    std::cout << "[SNV] Running pileup (this may take some time)\n";
    std::string cmd{"htsbox pileup -v -c -f " + user_args.args["--reference"] + " -C -Q20 -q5 " + user_args.args["-o"] + "/intermediate_output/filtered_unitigs_aln.srt.bam " + user_args.args["--normal-reads"] + " " + user_args.args["--tumor-reads"] + " > " + user_args.args["-o"] + "/intermediate_output/pileup.vcf"};
    user_args.args.insert({"--pileup", user_args.args["-o"] + "/intermediate_output/pileup.vcf"});
    system(cmd.c_str());
    return true;
}
