CC=g++
SDIR=./src
IDIR=./include
ODIR=./obj
BDIR=./bin
CFLAGS=-I$(IDIR) -std=c++11 -pthread

_DEPS = arguments.h bg.h tabfile.h smooth.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = find_tss.o arguments.o bg.o tabfile.o smooth.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

find_tss: $(OBJ)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
