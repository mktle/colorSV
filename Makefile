CC = g++

CFLAGS = -std=c++11 -O2 -DNEDEBUG -pedantic-errors -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion

colorSV: main.cpp argument_parser.cpp preprocess.cpp topology_search.cpp
	$(CC) $(CFLAGS) -o colorSV main.cpp argument_parser.cpp preprocess.cpp topology_search.cpp

clean: 
	rm colorSV
