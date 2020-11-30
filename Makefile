############################################# globel config #############################################
CC       = gcc
AR       = ar
CFLAGS   =  -g -O2 -Wno-unused-function
DFLAGS   =  -DUSE_MALLOC_WRAPPERS
LDFLAGS  =
LIBS     = -lz -lpthread
ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

############################################# objects #############################################
AOBJS   =   kstring.o \
			malloc_wrap.o \
			hfastq.o \
 			check2Files.o \
 			parallel_test.o \
 			samop.o \
 			sam2sfq.o \
 			progress.o \
 			eva2sam.o \
 			sync_pe.o \
 			kopen.o \
 			utils.o \
 			hsfq.o \
 			table.o \
 			kthread.o \
 			generate_cs.o \
 			reads2fa.o \
 			chr_extract.o \
 			proc_stat.o

PROGRAM = toy
all: $(PROGRAM)

############################################# dependency #############################################
.SUFFIXES:.c .o
.c.o:
		$(CC) $(CFLAGS) $(DFLAGS) -c -o $@ $<
$(PROGRAM): $(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) -o $@ $(AOBJS) main.o $(LDFLAGS) $(LIBS)

############################################# clean #############################################
clean:
	rm -f *.o $(PROGRAM) *.out
depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )
