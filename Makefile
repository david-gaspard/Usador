## Created on 2025-06-27 at 16:57:59 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Makefile of the Usador program.
##########################
## BUILD MAIN PROGRAM
##########################
PROGNAME = usador
SRCLIST = $(shell find src/ -name "*.cpp" ! -name "Main.cpp")
BINLIST = $(SRCLIST:src/%.cpp=bin/%.o)

CC=g++
CFLAGS=-W -Wall -Wextra -std=c++17 -Isrc -Jbin -fopenmp -O2
LIBS=-lumfpack -lpng

all: directories $(PROGNAME)

$(PROGNAME): $(BINLIST) src/Main.cpp
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

bin/%.o: src/%.cpp 
	$(CC) $(CFLAGS) -c $< -o $@

directories:
	mkdir -p bin/

##########################
## BUILD TEST EXECUTABLE
##########################
TESTSRCLIST = $(shell find test/ -name "Test*.cpp")
TESTEXELIST = $(TESTSRCLIST:test/Test%.cpp=%.test)

test: directories $(TESTEXELIST)

%.test: $(BINLIST) test/Test%.cpp
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

##################################
## CLEAN ALL BUILDS AND TESTS
##################################
clean:
	rm -rfv $(PROGNAME) bin/*.o $(TESTEXELIST)

###### END OF FILE ######
