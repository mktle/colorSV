#include "argument_parser.h"

#ifndef SNV_H
#define SNV_H

namespace snv{
	bool check_args(ArgumentParser& user_args);
	bool run_pileup(ArgumentParser& user_args);
	bool call_snvs(ArgumentParser& user_args);
}

const int GENOTYPE_IDX{0};
const int ADF_ALLELE1_IDX{1};
const int ADF_ALLELE2_IDX{2};
const int ADR_ALLELE1_IDX{3};
const int ADR_ALLELE2_IDX{4};

const int UNITIG_IDX{0};
const int NORMAL_IDX{1};
const int TUMOR_IDX{2};

#endif

