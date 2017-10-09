# Generated automatically from Makefile.in by configure.
.IGNORE:

SHELL = /bin/sh

# PATHS

prefix  = /home/bstephen/pamr/w

BBH_SYSTEM = LINUX_IA32

CC     = mpicc
AR     = ar
RANLIB = ranlib

LIBS       = -lpamrw -lamrdw $(MPILIB) -lm   -lbbhutil -lg2c -lm 
LDFLAGS    = -L../../lib  -L/usr/X11R6/lib  -L/usr/X11R6/lib -L. -L$(prefix)/lib    -L/home/bstephen/bin/rnpl.new/lib -L/usr/local/openmpi/1.2.5/gcc-ib/x86_64/lib64
CPPFLAGS   = -I.  -DLINUX -DPORT_GROUP -DWant_c_files -I../../include -I$(prefix)/include   -I/home/bstephen/bin/rnpl.new/include -I/usr/local/openmpi/1.2.5/gcc-ib/x86_64/include  -DHAVE_LIBM=1 -DHAVE_STRDUP=1   -DLINUX -DF77_HAS_SYSTEM -DF77_HAS_CHDIR -I$(prefix)/amrd/
CFLAGS     =  -O0 -g

CTOF     = mfi
PERL     = perl 

F77      = mpif77
F77FLAGS =  -O0 -g 

CC_COMP    = $(CC) $(CPPFLAGS) -c $(CFLAGS)
CC_LOAD    = $(CC) $(CFLAGS) $(LDFLAGS) 

F77_COMP   = $(F77) $(F77FLAGS) 
F77_LOAD   = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS) 

F77_TRANSFORM = touch

.f.o:
	$(F77_COMP) -c $*.f

.c.o:
	cat $*.c > $*_.c
	$(CC_COMP) -c $*_.c
	mv $*_.o $*.o

EXECUTABLES = srhd_w eos_test
OBJECTS = srhd_w.o num_w.o inversion_interface.o u2p_util.o utoprim_2d.o shen.o table_inversion.o eos_test.o
OBJECTS2 = num_w.o inversion_interface.o u2p_util.o utoprim_2d.o shen.o table_inversion.o eos_test.o

# targets

all: $(EXECUTABLES)

srhd_w.o: srhd_w.c num_w.h num_w.f shen.f table_inversion.f
eos_test.o: shen.f table_inversion.f eos_test.f

srhd_w: $(OBJECTS)
	$(F77_LOAD)  $(OBJECTS) $(LIBS) -o srhd_w
	cp srhd_w ../../bin
	cp srhd_w run/.

eos_test: $(OBJECTS2)
	$(F77_LOAD) $(OBJECTS2) $(LIBS) -o eos_test

clean:
	echo; echo "Cleaning up"
	rm *_.c > /dev/null 2>&1 
	rm *.o > /dev/null 2>&1
	rm $(EXECUTABLES) > /dev/null 2>&1
	rm *~ > /dev/null


