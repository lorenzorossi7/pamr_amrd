.IGNORE:

SHELL = /bin/sh

# PATHS

prefix  = /Users/pau/local

BBH_SYSTEM = DARWIN

CC     = mpicc
AR     = ar
RANLIB = ranlib

LIBS       = -lpamrw $(MPILIB) -lm -lmpi   -lbbhutil -lg2c -lm 
LDFLAGS    = -L../lib  -L/usr/X11R6/lib -L/sw/lib -L. -L$(prefix)/lib   -L/Users/pau/local/lib 
CPPFLAGS   = -I. -I../include  -DDARWIN -D__L -DWant_c_files -I$(prefix)/include   -I/Users/pau/local/rnpletal/rnpl/include/ -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DHAVE_LIBMPI=1 -DHAVE_LIBM=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STRDUP=1  -DF77_HAS_SYSTEM -DF77_HAS_CHDIR
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

EXECUTABLES = test1_w test2_w test3_w

# targets

all: $(EXECUTABLES)

.c.o:
	cat $*.c > $*_.c
	$(CC_COMP) -c $*_.c
	mv $*_.o $*.o

test1_w: test1_w.o 
	$(F77_LOAD) test1_w.o $(LIBS) -o test1_w
	cp test1_w ../bin

test2_w: test2_w.o 
	$(F77_LOAD) test2_w.o $(LIBS) -o test2_w
	cp test2_w ../bin

test3_w: test3_w.o 
	$(F77_LOAD) test3_w.o $(LIBS) -o test3_w
	cp test3_w ../bin

#test1_w.o: ../include/pamr.h ../lib/libpamr.a ../include/internal_opts.h test1.c

clean:
	echo; echo "Cleaning up"
	for f in $(AUX); do \
		(cd $$f; make clean) \
	done
	(cd bin; rm $(EXECUTABLES) > /dev/null 2>&1) 
	rm $(EXECUTABLES) > /dev/null 2>&1
	rm *_.c > /dev/null 2>&1
	rm *.o > /dev/null 2>&1
	rm *~ > /dev/null
	rm *.sdf > /dev/null 2>&1
	rm config.cache config.log config.status > /dev/null 2>&1
	rm work.pc work.pcl > /dev/null 2>&1
	rm Makefile > /dev/null 2>&1
