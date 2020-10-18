############################################# 通用编译配置 #############################################
CC       = gcc
AR       = ar
CFLAGS   =  -g -O2 -Wno-unused-function
DFLAGS   =  -DUSE_MALLOC_WRAPPERS
LDFLAGS  =
LIBS     = -lz -lpthread
ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif
# {CFLAGS：c编译选项} {CPPFLAGS：c++编译选项} {DFLAGS：预编译选项} {LDFLAGS：libs directory} {LIBS：libs}

############################################# toy objects #############################################
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

############################################# All programs #############################################
PROGRAMS = toy-1.1
all: $(PROGRAMS)

############################################# toy program #############################################
.SUFFIXES:.c .o

.c.o:
		$(CC) $(CFLAGS) $(DFLAGS) -c -o $@ $<

$(PROGRAMS): $(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) -o $@ $(AOBJS) main.o $(LDFLAGS) $(LIBS)

############################################# clean #############################################
clean:
	rm -f *.o $(PROGRAMS) *.out

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.