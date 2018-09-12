.IGNORE:

SHELL = /bin/sh

prefix = /usr/local
BBH_SYSTEM = DARWIN
bindir = ../../bin

BBH_SYSTEM = DARWIN

RNPL   = rnpl

AR     = ar
RANLIB = ranlib

F77_TRANSFORM = touch

LIBS       = -lpamrw -lamrdw $(MPILIB) -lm -lmpi   -lbbhutil -lg2c -lm 
LDFLAGS    = -L../../lib  -L/usr/X11R6/lib -L/sw/lib -L. -L$(prefix)/lib  \
               -L/usr/local/lib

CC       = mpicc
CFLAGS     = -g -O2 
CPPFLAGS = -I.  -DDARWIN -D__L -DWant_c_files -I../../include -I$(prefix)/include \
             -I/usr/local/include

CC_COMP  = $(CC) $(CPPFLAGS) -c $(CFLAGS)
CC_LOAD  = $(CC) $(CFLAGS) $(LDFLAGS) 

F77      = mpif77
F77FLAGS = -O6 -fno-second-underscore

F77_COMP   = $(F77) -c $(F77FLAGS) 
F77_LOAD   = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS) 

EXECUTABLES     = wave_w
# Miscellaneous files to clean up
MISCDATAFILES   =

SRC = *.f *.inc

all: $(EXECUTABLES)

.f.o:
	$(F77_COMP) $*.f 

.c.o:
	cat $*.c > $*_.c
	$(CC_COMP) -c $*_.c
	mv $*_.o $*.o

all: $(EXECUTABLES)

install: all

full: install confidence_tests

confidence_tests:
	echo "Not implemented yet"

translate: 
	touch translate

wave_w.o: wave_w.h 

WAVE_OBJS = wave_w.o sf_evo_w.o 
wave_w: $(WAVE_OBJS) 
	$(F77_LOAD) $(WAVE_OBJS) $(LIBS) -o wave_w
	cp wave_w $(bindir)

########################################################################
# Clean-up
########################################################################
clean:
	/bin/rm $(EXECUTABLES)
	/bin/rm *_.c > /dev/null 2>&1 
	/bin/rm *.o 
	/bin/rm *.a
	/bin/rm *.sdf
	/bin/rm *.segdat
	/bin/rm *~
	/bin/rm  Makefile config.cache config.log config.status
	/bin/rm $(MISCDATAFILES)
