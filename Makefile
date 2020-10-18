CC=         gcc
#CC=			clang --analyze
CFLAGS=     -g -O2 -Wall -Wno-unused-function
DFLAGS=     -DUSE_MALLOC_WRAPPERS
AOBJS=      kstring.o SAMCheck.o malloc_wrap.o getReadsInfo.o \
 			sam2sfq.o check2Files.o temp.o multiThread.o \
 			samop.o wrapper.o sam2pe_sfq.o merge_fastq.o \
 			harc2fq.o

PROG=		toy-1.0
INCLUDES=
LIBS=       -lpthread
SUBDIRS=	.

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

$(PROG):$(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main.o -o $@ -L. $(LIBS)

clean:
	rm -f *.o *.out $(PROG)
check:
	rm -f SAMCheck.o
mul:
	rm -f multiThread.o

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.