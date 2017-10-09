# Generated automatically from Makefile.in by configure.
.IGNORE:

SHELL = /bin/sh

# PATHS
prefix  = /home/users/matt
bindir = $(prefix)/bin
includedir = $(prefix)/include
libdir = $(prefix)/lib 

BBH_SYSTEM = OSF1

CC     = cxx
AR     = ar
RANLIB = ranlib

DEFS       =  -DHAVE_LIBMPI=1 -DHAVE_LIBM=1   -std1 -DDEC_ALPHA -DF77_HAS_SYSTEM -DF77_HAS_CHDIR

CTOF     = mfi

PERL     = perl

F77      = f90
F77FLAGS = -fast

LIBS       =   -lpamr $(MPILIB) -lm -lmpi   -lbbhutil  -lsv -lfor -lFutil -lots -lUfor -lm
LDFLAGS    =    -L$(PAMR_HOME)/lib -L.   -L/home/users/matt/lib 

CFLAGS     = -c -g  -x cxx
CPPFLAGS = -I. -I$(PAMR_HOME)/include/  -std1 -DDEC_ALPHA -I$(prefix)/include   -I/home/users/matt/include

CC_COMP    = $(CC) $(CFLAGS) $(CPPFLAGS)
CC_LOAD    = $(CC) $(CFLAGS) $(LDFLAGS)

F77_COMP   = $(F77) $(F77FLAGS)
F77_LOAD   = $(F77) $(F77FLAGS) $(LDFLAGS)

F77_TRANSFORM = touch

.f.o:
	   $(F77_COMP) -c $*.f

.c.o:
	   $(CC_COMP) -c $*.c

# targets

all: libpamr.a

OBJECTS = pamr.o gh.o io.o transfer.o misc.o m_util_r8.o
INCS = gh.h misc.h $(PAMR_HOME)/include/pamr.h $(PAMR_HOME)/include/internal_opts.h io.h transfer.h $(PAMR_HOME)/include/m_util_r8.h

pamr.o: $(INCS) pamr.c
gh.o: $(INCS) gh.c
io.o: $(INCS) io.c
misc.o: $(INCS) misc.c
transfer.o: $(INCS) transfer.c

install: libpamr.a 
	if test ! -d $(prefix); then\
		mkdir $(prefix);\
	fi
	if test ! -d $(bindir); then\
		mkdir $(bindir);\
	fi
	if test ! -d $(includedir); then\
		mkdir $(includedir);\
	fi
	if test ! -d $(libdir); then\
		mkdir $(libdir);\
	fi
	cp ../lib/* $(libdir)
	cp ../include/* $(includedir)


libpamr.a: $(OBJECTS)
	$(AR) r libpamr.a $(OBJECTS)
	/bin/cp libpamr.a $(PAMR_HOME)/lib

clean:
	echo; echo "Cleaning up"
	rm *.o > /dev/null 2>&1
	rm *~ > /dev/null
	rm libpamr.a > /dev/null 2>&1
	rm *.sdf > /dev/null 2>&1
	rm config.cache config.log config.status > /dev/null 2>&1
	rm Makefile > /dev/null 2>&1
