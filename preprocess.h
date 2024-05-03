#include "argument_parser.h"

#ifndef PREPROCESS_H
#define PREPROCESS_H

namespace preprocess{
	bool file_setup(ArgumentParser& user_args);
	bool align_unitigs(ArgumentParser& user_args);
	bool filter_unitigs(ArgumentParser& user_args, bool is_r_utg);
	bool check_args(ArgumentParser& user_args);
}

#endif

