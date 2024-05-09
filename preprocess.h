#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "argument_parser.h"

namespace preprocess{
	bool file_setup(ArgumentParser& user_args);
	bool align_unitigs(ArgumentParser& user_args);
	bool filter_unitigs(ArgumentParser& user_args);
	bool check_args(ArgumentParser& user_args);
}

#endif

