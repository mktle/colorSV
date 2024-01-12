#include <map>
#include <string>

#ifndef ARGUMENT_PARSER_H
#define ARGUMENT_PARSER

class ArgumentParser{
	public:
		ArgumentParser(int &argc, char**argv);

		std::map<std::string, std::string> args;
};

void print_help();

#endif

