EXEC = tssp
CC = icpx
SRCDIR = src/
MKLROOT = /opt/intel/oneapi/mkl/latest
MKLINCLUDE = $(MKLROOT)/include
LOCALINCLUDE = /usr/local/include
LOCALLIB = /usr/local/lib
LIB = lib/

SRCEXT := cpp
SRC_FILES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))

CFLAGS = -c -g -traceback -Wall -DMKL_ILP64 -qopenmp -O3 -xhost -fbuiltin -ipo -no-ftz -static-intel -std=c++11 -qmkl=parallel -I$(MKLINCLUDE) -I$(LOCALINCLUDE)
LFLAGS = -ipo -L$(MKLROOT)/lib/intel64 -L$(LOCALLIB) -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lcfitsio -lm -L$(LIB)

O_FILES = $(SRC_FILES:.cpp=.o) $(INI_SRC:.cpp=.o)

$(EXEC): $(O_FILES)
	$(CC) -o $@ $(O_FILES) $(LFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(O_FILES) $(EXEC)

print-%  : ; @echo $* = $($*)
