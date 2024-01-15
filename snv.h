#include "argument_parser.h"

#ifndef SNV_H
#define SNV_H

namespace snv{
	bool check_args(ArgumentParser& user_args);
	bool run_pileup(ArgumentParser& user_args);
}

#endif

