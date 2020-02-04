.IGNORE:

SHELL = /bin/sh

# PATHS

prefix  = /Users/pau/local

BBH_SYSTEM = DARWIN

CC     = mpicc
AR     = ar
RANLIB = ranlib

LIBS       = -lpamrw -lamrdw $(MPILIB) -lm -lmpi   -lbbhutil -lg2c -lm
LDFLAGS    = -L../../lib  -L/usr/X11R6/lib -L/sw/lib  -L/usr/X11R6/lib -L/sw/lib -L. -L$(prefix)/lib    -L/Users/pau/local/lib 
CPPFLAGS   = -I.  -DDARWIN -D__L -DWant_c_files -I../../include -I$(prefix)/include   -I/Users/pau/local/include -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DHAVE_LIBMPI=1 -DHAVE_LIBM=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STRDUP=1  -DF77_HAS_SYSTEM -DF77_HAS_CHDIR
CFLAGS     = -g -O2 

CTOF     = mfi
PERL     = perl 

F77      = mpif77
F77FLAGS = -O6 -fno-second-underscore

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

EXECUTABLES = nbs_w
OBJECTS = nbs_w.o num_w.o

# targets

all: $(EXECUTABLES)

nbs_w.o: nbs_w.c num_w.h num_w.f 

nbs_w: $(OBJECTS)
	$(F77_LOAD) $(OBJECTS) $(LIBS) -o nbs_w
	cp nbs_w ../../bin

clean:
	echo; echo "Cleaning up"
	rm *_.c > /dev/null 2>&1 
	rm *.o > /dev/null 2>&1
	rm $(EXECUTABLES) > /dev/null 2>&1
	rm *~ > /dev/null
	rm *.sdf > /dev/null 2>&1
	rm config.cache config.log config.status > /dev/null 2>&1
	rm Makefile > /dev/null 2>&1
