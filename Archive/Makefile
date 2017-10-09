.IGNORE:

SHELL = /bin/sh

# PATHS
prefix = /usr/local
bindir = $(prefix)/bin
libdir = $(prefix)/lib

PAMR_HOME = /d/bh4/home/fransp/rd_map/pamr/

# Programs
CTOF = mfi
PERL = perl
RPCGEN = rpcgen
RANLIB = ranlib
AR = ar
CC = mpicc
F77 = mpif77

BBH_SYSTEM=LINUX

# Compiler Commands
F77_TRANSFORM = touch

LIBS       =   -lpamr -lm -lbbhutil -lsv
LDFLAGS    =   -L/usr/local/lib  -L/usr/lib -L. -L$(PAMR_HOME)/lib

CFLAGS   = -c 
CPPFLAGS =  -DLINUX -D__ -DPORT_GROUP -DWant_c_files   -I/usr/local/include -I$(PAMR_HOME)/include/

CC_COMP  = $(CC) $(CFLAGS) $(CPPFLAGS)
CC_LOAD  = $(CC) $(LDFLAGS) -L$(PAMR_HOME)/lib/

F77FLAGS = -g

F77_COMP   = $(F77) -c $(F77FLAGS)
F77_LOAD   = $(F77) $(F77FLAGS) $(LDFLAGS)

.f.o:
	$(F77_COMP) -c $*.f

.c.o:
	$(CC_COMP) -c $*.c

AUX = src

EXECUTABLES = test1
LIBFILES =

# targets

all: 
	echo; echo "Making in $(AUX)"
	for f in $(AUX); do \
		(cd $$f; make) \
	done
	make test1

test1: test1.o $(PAMR_HOME)/include/pamr.h $(PAMR_HOME)/lib/libpamr.a $(PAMR_HOME)/include/internal_opts.h
	$(CC_COMP) test1.c
	$(CC_LOAD) test1.o $(LIBS) -o test1

clean:
	echo; echo "Cleaning up"
	for f in $(AUX); do \
		(cd $$f; make clean) \
	done
	(cd bin; rm $(EXECUTABLES) > /dev/null 2>&1) 
	rm $(EXECUTABLES) > /dev/null 2>&1
	rm *.o > /dev/null 2>&1
	rm *~ > /dev/null
	rm *.sdf > /dev/null 2>&1
#	rm config.cache config.log config.status > /dev/null 2>&1
#	rm Makefile > /dev/null 2>&1
