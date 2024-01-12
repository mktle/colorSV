CC = g++

CFLAGS = -std=c++11 -O2 -DNEDEBUG -pedantic-errors -Wall -Weffc++ -Wextra -Wconversion -Wsign-conversion

sv-caller: main.o argument_parser.o
	$(CC) $(CFLAGS) -o sv-caller main.o argument_parser.o 

clean: 
	rm *.o sv-caller
