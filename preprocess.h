#include "argument_parser.h"

#ifndef PREPROCESS_H
#define PREPROCESS_H

namespace preprocess{
	bool file_setup(ArgumentParser& args);
	bool filter_unitigs(ArgumentParser& args);
	bool check_args(ArgumentParser& args);
}

#endif

