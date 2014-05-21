ifdef HTSLIB
	HTSARGS=-I $(HTSLIB)/htslib -D_USESAM=1 $(HTSLIB)/libhts.a
	LIBS=-lpthread -lz -lm
else
	HTSARGS=
	LIBS=-lz -lm
endif

ifdef DEBUG
	OPT=-O0 -g -ggdb -DDEBUG=1
else
	OPT=-O2
endif

all: readsim

readsim: readsim.c seq_file.h stream_buffer.h
	$(CC) -Wall -Wextra $(OPT) -o readsim readsim.c $(HTSARGS) $(LIBS)

plot: data/PhiX.1K.1.pdf

data/PhiX.1K.1.out: readsim
	./readsim data/PhiX.1.fq out.1.fq.gz out.2.fq.gz > data/PhiX.1K.1.out

data/PhiX.1K.1.csv: data/PhiX.1K.1.out
	grep '^QProb:' data/PhiX.1K.1.out | sed 's/^QProb: //g' > data/PhiX.1K.1.csv

data/PhiX.1K.1.pdf: plot.R data/PhiX.1K.1.csv
	R -f plot.R --args data/PhiX.1K.1.csv data/PhiX.1K.1.pdf

clean:
	rm -rf readsim *.greg *.dSYM data/PhiX.1K.1.out data/PhiX.1K.1.csv data/PhiX.1K.1.pdf

.PHONY: all clean plot
