EXEC = tssp
CC = icpc
SRCDIR = src/
INIDIR = inih/
MKLROOT = /opt/intel/mkl
MKLINCLUDE = $(MKLROOT)/include
LOCALINCLUDE = /usr/local/include
LOCALLIB = /usr/local/lib
LIB = lib/

SRCEXT := cpp
SRC_FILES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
INI_SRC := $(shell find $(INIDIR) -type f -name *.$(SRCEXT))

CFLAGS = -c -g -fvar-tracking  -traceback -Wall -DMKL_ILP64 -qopenmp -fast -O3 -xhost -ip -fbuiltin -ipo -no-ftz -static-intel -std=c++11 -mkl=parallel -I$(MKLINCLUDE) -I$(LOCALINCLUDE)
LFLAGS = -L$(MKLROOT)/lib/intel64 -L$(LOCALLIB) -parallel -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lcfitsio -lm -L$(LIB) -linih

O_FILES = $(SRC_FILES:.cpp=.o) $(INI_SRC:.cpp=.o) $(INIREADER_SRC)

print-%  : ; @echo $* = $($*)

$(EXEC): $(O_FILES)
	$(CC) -o $@ $(O_FILES) $(LFLAGS)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
