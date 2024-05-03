CC = g++

CFLAGS = -std=c++11 -O2 -DNEDEBUG -pedantic-errors -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion

sv-caller: main.cpp argument_parser.cpp preprocess.cpp topology_search.cpp
	$(CC) $(CFLAGS) -lboost_filesystem -lboost_system -o sv-caller main.cpp argument_parser.cpp preprocess.cpp topology_search.cpp

clean: 
	rm sv-caller
