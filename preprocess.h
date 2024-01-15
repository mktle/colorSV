#include "argument_parser.h"

#ifndef PREPROCESS_H
#define PREPROCESS_H

namespace preprocess{
	bool file_setup(ArgumentParser& user_args);
	bool filter_regions(ArgumentParser& user_args);
	bool filter_unitigs(ArgumentParser& user_args);
	bool check_args(ArgumentParser& user_args);
}

#endif

