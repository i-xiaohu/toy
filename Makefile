CC = g++

CPPFLAGS = -g -pg -O2 -Wall -Wno-unused-function

LIBS =

SRCS =  calcCov.cpp ctgals.cpp DecCheck.cpp fetch.cpp getReads_n.cpp main.cpp metaCheck.cpp \
		SAMCheck.cpp samop.cpp splitPE.cpp

OBJS = $(SRCS: .cpp = .o)

EXEC = toy

%.o : %.cpp
	$(CC) -c $(CPPFLAGS) $<

$(EXEC) : $(OBJS)
	$(CC) $(CPPFLAGS)  $^ -o $@ $(LIBS)

clean:
	rm -f *.o $(EXEC) *.out

check:
	rm -f SAMCheck.o