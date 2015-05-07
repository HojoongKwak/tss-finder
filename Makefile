CC=g++
SDIR=./src
IDIR=./include
ODIR=./obj
BDIR=./bin
CFLAGS=-I$(IDIR) -std=c++11 -pthread

_SRCC = find_tss find_tss_func
_SDEP = find_tss
_DEPS = arguments bg tabfile smooth bed

SSRC = $(patsubst %,$(SDIR)/%.cpp,$(_SRCC))
SOBJ = $(patsubst %,$(ODIR)/%.o,$(_SRCC))
SDEP = $(patsubst %,$(IDIR)/%.h,$(_SDEP))
DEPS = $(patsubst %,$(IDIR)/%.h,$(_DEPS))
OBJS = $(patsubst %,$(ODIR)/%.o,$(_DEPS))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS) $(SDEP)
	$(CC) -c -o $@ $< $(CFLAGS)

find_tss: $(SOBJ) $(OBJS)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
