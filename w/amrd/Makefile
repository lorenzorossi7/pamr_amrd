.IGNORE:

SHELL = /bin/sh

# PATHS

prefix  = /Users/pau/local

BBH_SYSTEM = DARWIN

CC     = mpicc
AR     = ar
RANLIB = ranlib

LIBS       = -lpamrw $(MPILIB) -lm -lmpi   -lbbhutil -lg2c -lm 
LDFLAGS    = -L../lib  -L/usr/X11R6/lib -L/sw/lib -L.   -L/Users/pau/local/lib 
CPPFLAGS   = -I. -I../include  -DDARWIN -D__L -DWant_c_files -I$(prefix)/include   -I/Users/pau/local/rnpletal/rnpl/include/ -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DHAVE_LIBMPI=1 -DHAVE_LIBM=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STRDUP=1  -DF77_HAS_SYSTEM -DF77_HAS_CHDIR
CFLAGS     = -g -O2 

CTOF     = mfi
PERL     = perl 

F77      = mpif77
F77FLAGS = -O6 -fno-second-underscore

CC_COMP    = $(CC) $(CPPFLAGS) -c $(CFLAGS)
CC_LOAD    = $(CC) $(CFLAGS) $(LDFLAGS) 

F77_COMP   = $(F77) $(F77FLAGS) 
F77_LOAD   = $(F77) $(F77FLAGS) $(LDFLAGS) 

F77_TRANSFORM = touch

OBJECTS = amrd_w.o globals_w.o io_w.o util_w.o evolve_w.o mg_w.o cls_w.o regrid_script_w.o
HEADERS = globals_w.h io_w.h util_w.h evolve_w.h mg_w.h cls_w.h ../include/amrd_w.h regrid_script_w.h ../include/m_util_r8_w.h

.f.o:
	$(F77_COMP) -c $*.f

.c.o:
	cat $*.c > $*_.c
	$(CC_COMP) -c $*_.c
	mv $*_.o $*.o

# targets

all: libamrdw.a 

amrd_w.o: amrd_w.c $(HEADERS)
globals_w.o: globals_w.c $(HEADERS)
io_w.o: io_w.c $(HEADERS)
mg_w.o: mg_w.c $(HEADERS)
cls_w.o: cls_w.c $(HEADERS)
evolve_w.o: evolve_w.c $(HEADERS)
util_w.o: util_w.c $(HEADERS)
regrid_script_w.o: regrid_script_w.c $(HEADERS)

libamrdw.a: $(OBJECTS)
	$(AR) r libamrdw.a $(OBJECTS) 
	cp libamrdw.a ../lib/

clean:
	echo; echo "Cleaning up"
	rm libamrd.a > /dev/null 2>&1
	rm *_.c > /dev/null 2>&1
	rm *.o > /dev/null 2>&1
	rm *~ > /dev/null
	rm *.sdf > /dev/null 2>&1
	rm config.cache config.log config.status > /dev/null 2>&1
	rm work.pc work.pcl > /dev/null 2>&1
	rm Makefile > /dev/null 2>&1
