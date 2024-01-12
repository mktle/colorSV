CC = g++

CFLAGS = -std=c++0x -O2 -DNEDEBUG -pedantic-errors -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion

sv-caller: main.cpp argument_parser.cpp
	$(CC) $(CFLAGS) -o sv-caller main.cpp argument_parser.cpp

clean: 
	rm *.o sv-caller
