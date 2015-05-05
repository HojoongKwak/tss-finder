CC=g++
SDIR=./src
IDIR=./include
ODIR=./obj
BDIR=./bin
CFLAGS=-I$(IDIR) -std=c++11 -pthread

_SOBJ = find_tss.o
_DEPS = arguments.h bg.h tabfile.h smooth.h
_OBJS = arguments.o bg.o tabfile.o smooth.o

SOBJ = $(patsubst %,$(ODIR)/%,$(_SOBJ))
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

find_tss: $(SOBJ) $(OBJS)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
