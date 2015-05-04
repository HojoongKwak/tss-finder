CC=g++
IDIR=./include
ODIR=./obj
CFLAGS=-I$(IDIR) -std=c++11 -pthread

_DEPS = arguments.h bg.h tabfile.h smooth.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = arguments.o bg.o tabfile.o smooth.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

find_tss: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
