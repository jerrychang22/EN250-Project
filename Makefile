SRC :=       \
    mod.cpp  \
#    stat.c \
#	decode.c \

TARGETS	:= mod2 #stat mod-water decode

D1 = "NODEFINED"
D2 = "NODEFINED"

#################################
# Generic part
BUILD_ROOT	:= $(shell pwd)
#include $(BUILD_ROOT)/Makefile.inc
#################################

#gcc -Wall mod.c mod.h structs.h constants.h -fopenmp -lm -O3 -o work3
CC	:= g++
#CC	:= gcc
#OPTFLAGS := 
#OPTFLAGS := -v
#OPTFLAGS := -g
OPTFLAGS := -O3
#OPTFLAGS	:= -O3 -fopenmp 
CFLAGS	:= -Wall -std=c++0x -lstdc++ -lboost_program_options -Iinclude
#LDFLAGS	:= -lm -Wl,-rpath . -Wl,--dynamic-linker=include/ld-linux.so.2
LDFLAGS	:= -lm
#PNGWRITERFLAGS := -Ipngwriter/pngwriter/build/include -L. -L/usr/lib/i386-linux-gnu -lfreetype -lz -I/usr/include/freetype2 -lpng -lpngwriter 
PNGWRITERFLAGS := -L. -Linclude -L/usr/lib/i386-linux-gnu -lfreetype -lz -I/usr/include/freetype2 -lpng -lpngwriter 
Q           = @
OBJ	:= $(SRC:%.cpp=%.o)
RM		:= rm -f

mod2: mod.cpp mod.hpp structs.hpp species.o simulation.o simulation.hpp species.hpp gitversion_tmpfile.h 
	$(CC) mod.cpp simulation.hpp simulation.o species.hpp species.o mod.hpp structs.hpp gitversion_tmpfile.h -D $(D1) -D $(D2) $(CFLAGS) $(LDFLAGS) $(PNGWRITERFLAGS) $(OPTFLAGS) -o $@
	@rm gitversion_tmpfile.h

gitversion_tmpfile.h:
	@echo "#define GIT_VERSION \"`git describe  --dirty --always --long`\"" > gitversion_tmpfile.h

simulation.o: simulation.hpp simulation.cpp species.hpp species.o gitversion_tmpfile.h
	$(CC) species.hpp -D $(D1) -D $(D2) $(CFLAGS) $(PNGWRITERFLAGS) $(OPTFLAGS) -c simulation.cpp

species.o: simulation.hpp species.hpp species.cpp
	$(CC) -D $(D1) -D $(D2) $(CFLAGS) $(PNGWRITERFLAGS) $(OPTFLAGS) -c species.cpp
	

clean:
	$(RM) gitversion_tmpfile.h *.o .*.swp *~ *.bak *.dep *.gch tags $(TARGETS)

#################################
# Include Makefile.dep if it exist
#-include Makefile.dep
#################################
