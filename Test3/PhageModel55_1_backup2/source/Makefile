CC = gcc
CXX = g++

# Name of your program
PROJECT = demo

#############################################
# DON'T FORGET TO CHANGE DEPENDENCY LINES!! #
#############################################
# C++ both source (.cpp) and header (.hpp)
CCBOTH = cash-display automaton runtime
# C++ only source (.cpp)
CCSOURCE = main
# C++ only header (.hpp)
CCHEADER = assert cellular-automata
# C both source (.c) and header (.h)
CBOTH = 
# C only source (.c)
CSOURCE = arithmetic basic color filter io logical margolus movie neighbors noise png ps random shift x11
# C only header (.h)
CHEADER = cash
# Other files to be archived
OTHERS = Makefile

# Options to compiler (both for CC and CXX)
CCOPT = -O3 -std=c++11 -Wall -DNDEBUG
COPT = -O3 -Wall -DNDEBUG
LIBS =  -lpng -lX11


# All object files that should be generated
OBJALL = $(addsuffix .o, $(CCBOTH) $(CCSOURCE) $(CBOTH) $(CSOURCE))

# Link all files to generate a program
all: $(OBJALL) source.tar.gz
	$(CXX) $(OBJALL) $(CCOPT) -o $(PROJECT) $(LDFLAGS) $(LIBS) $(LDIR)

# Dependency of files. Add/modify if necessarly (all object files depend on Makefile)
$(OBJALL): Makefile

# Cash
arithmetic.o: cash.h
basic.o: cash.h
color.o: cash.h
filter.o: cash.h
io.o: cash.h
logical.o: cash.h
margolus.o: cash.h
movie.o: cash.h
neighbors.o: cash.h
noise.o: cash.h
png.o: cash.h
ps.o: cash.h
random.o: cash.h
shift.o: cash.h
x11.o: cash.h

# My Libraries
cash-display.o: Makefile cash-display.hpp
runtime.o: Makefile runtime.hpp

# Other sources
COMMON = Makefile assert.hpp cellular-automata.hpp 
main.o: $(COMMON)


# Make an archive containing EVERYTHING
EVERYTHING = $(addsuffix .cpp, $(CCBOTH) $(CCSOURCE)) $(addsuffix .hpp, $(CCBOTH) $(CCHEADER)) $(addsuffix .c, $(CBOTH) $(CSOURCE)) $(addsuffix .h, $(CBOTH) $(CHEADER)) $(OTHERS)
source.tar.gz: $(EVERYTHING)
	tar -zcf source.tar.gz $(EVERYTHING)

# Make rules
%.o : %.c
	$(CC) -c $(COPT) $(IDIR) $< -o $@

%.o : %.cpp
	$(CXX) -c $(CCOPT) $(IDIR) $< -o $@

clean:
	rm *.o demo
