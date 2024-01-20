#include "argument_parser.h"
#include "snv.h"

#include <assert.h>
#include <fstream>
#include <iostream>
#include <limits>
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

bool snv::call_snvs(ArgumentParser& user_args){
    assert(TUMOR_IDX == NORMAL_IDX + 1);
    assert(ADF_ALLELE2_IDX == ADF_ALLELE1_IDX + 1);
    assert(ADR_ALLELE2_IDX == ADR_ALLELE1_IDX + 1);

    std::cout << "[SNV] Calling SNVs\n";
    // open input and output files
    std::ifstream pileup_file(user_args.args["--pileup"]);
    std::ofstream out_vcf(user_args.args["-o"] + "/called_snvs.vcf");

    int read_thresh{std::stoi(user_args.args["-m"])};
    float frac_thresh{stof(user_args.args["--frac"])};

    // copy meta-information and header into output .vcf file
    std::string pileup_line;
    while (std::getline(pileup_file, pileup_line) && pileup_line[0] == '#' && pileup_line[1] == '#'){
        out_vcf << pileup_line << '\n';
    }
    out_vcf << pileup_line << '\n';

    std::string chr;
    std::string pos;
    std::string id;
    std::string ref;
    std::string alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::string format;

    // stores info for three categories (unitig/normal sample/tumor sample) with 5 values each
    std::string read_info[3][5];

    int total_support[3][2] = {0, 0, 0, 0, 0, 0};

    while (pileup_file >> chr){
        pileup_file >> pos >> id >> ref >> alt;
        if (chr.length() <= 5 && ref.length() == 1 && alt.length() == 1){
            pileup_file >> qual >> filter >> info >> format;
            pileup_file >> std::ws;

            // check that the unitig genotype supports an alt allele
            std::getline(pileup_file, read_info[UNITIG_IDX][GENOTYPE_IDX], ':');
            if (read_info[UNITIG_IDX][GENOTYPE_IDX] != "./." && read_info[UNITIG_IDX][GENOTYPE_IDX] != "0/0"){
                std::getline(pileup_file, read_info[UNITIG_IDX][ADF_ALLELE1_IDX], ',');
                std::getline(pileup_file, read_info[UNITIG_IDX][ADF_ALLELE2_IDX], ':');

                std::getline(pileup_file, read_info[UNITIG_IDX][ADR_ALLELE1_IDX], ',');

                pileup_file >> read_info[UNITIG_IDX][ADR_ALLELE2_IDX];
                pileup_file >> std::ws;

                for (int i{NORMAL_IDX}; i <= TUMOR_IDX; i++){
                    std::getline(pileup_file, read_info[i][GENOTYPE_IDX], ':');
                    std::getline(pileup_file, read_info[i][ADF_ALLELE1_IDX], ',');
                    std::getline(pileup_file, read_info[i][ADF_ALLELE2_IDX], ':');

                    std::getline(pileup_file, read_info[i][ADR_ALLELE1_IDX], ',');
                    pileup_file >> read_info[i][ADR_ALLELE2_IDX];
                    pileup_file >> std::ws;
                }


                // calculate total read support per sample, per allele
                for (int i{0}; i < 3; i++){
                    for (int j{0}; j < 2; j++){
                        total_support[i][j] = std::stoi(read_info[i][ADF_ALLELE1_IDX + j]) + std::stoi(read_info[i][ADR_ALLELE1_IDX + j]);

                    }
                }

                /*
                for (int i{0}; i < 5; i++){
                    std::cout << read_info[0][i] << ' ' << read_info[1][i] << ' ' << read_info[2][i] << '\n';
                }
                std::cout << "****************\n";

                for (int i{0}; i < 3; i++){
                    for (int j{0}; j < 2; j++){
                        std::cout << total_support[i][j] << ' ';
                    }
                    std:: cout << '\n';
                }

                assert(0);
                */

                bool to_write = false;
                // check each allele to see if:
                // 1) there is unitig support
                // 2) there is tumor support and no normal support
                for (int i{0}; i < 2; i++){
                    // unitig support check
                    if (total_support[UNITIG_IDX][i] > 0
                        // normal support check
                        && total_support[NORMAL_IDX][i] == 0
                        // tumor support check
                        && total_support[TUMOR_IDX][i] >= read_thresh){

                        // check that the proportion of tumor reads is greater than "frac_thresh"
                        if ( (float)(total_support[TUMOR_IDX][i]) / (float)(total_support[TUMOR_IDX][i] + total_support[NORMAL_IDX][i]) >= frac_thresh){
                            to_write = true;
                        }
                    }
                }
                if (to_write){
                    out_vcf << chr << '\t' << pos << '\t' << id << '\t' << ref << '\t' << alt << '\t' << qual << '\t' << filter << '\t' << info << '\t' << format << '\t';
                    for (int i{0}; i < 3; i++){
                        out_vcf << read_info[i][GENOTYPE_IDX] << ':'
                                << read_info[i][ADF_ALLELE1_IDX] << ','
                                << read_info[i][ADF_ALLELE2_IDX] << ':'
                                << read_info[i][ADR_ALLELE1_IDX] << ','
                                << read_info[i][ADR_ALLELE2_IDX];
                        if (i != 2){
                            out_vcf << '\t';
                        }
                    }
                    out_vcf << '\n';
                }
            }else{
                pileup_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }else{
            pileup_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
    return true;
}
