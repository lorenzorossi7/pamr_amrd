.IGNORE:

SHELL = /bin/sh

# PATHS
prefix  = /usr/local
bindir = $(prefix)/bin
includedir = $(prefix)/include
libdir = $(prefix)/lib 

BBH_SYSTEM = DARWIN

CC     = mpicc
AR     = ar
RANLIB = ranlib


CTOF     = mfi

PERL     = perl

F77      = mpif77
F77FLAGS = -O6 -fno-second-underscore

LIBS       =   -lpamrw $(MPILIB) -lm -lmpi   -lbbhutil -lg2c -lm
LDFLAGS    =    -L/usr/X11R6/lib -L/sw/lib -L../lib -L.   -L/usr/local/lib 

CFLAGS     = -g -O2
CPPFLAGS = -I. -I../include/  -DDARWIN -D__L -DWant_c_files -I$(prefix)/include   -I/usr/local/include -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DHAVE_LIBMPI=1 -DHAVE_LIBM=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STRDUP=1  -DF77_HAS_SYSTEM -DF77_HAS_CHDIR

CC_COMP    = $(CC) $(CPPFLAGS) -c $(CFLAGS)
CC_LOAD    = $(CC) $(CFLAGS) $(LDFLAGS)

F77_COMP   = $(F77) $(F77FLAGS)
F77_LOAD   = $(F77) $(F77FLAGS) $(LDFLAGS)

F77_TRANSFORM = touch

.f.o:
	$(F77_COMP) -c $*.f

.c.o:
	cat $*.c > $*_.c
	$(CC_COMP) -c $*_.c
	mv $*_.o $*.o

# targets

all: libpamrw.a

OBJECTS = pamr_w.o gh_w.o io_w.o transfer_w.o misc_w.o m_util_r8_w.o m_util_r8_c_w.o
INCS = gh_w.h misc_w.h ../include/pamr_w.h ../include/internal_opts_w.h io_w.h transfer_w.h ../include/m_util_r8_w.h

pamr_w.o: $(INCS) pamr_w.c
gh_w.o: $(INCS) gh_w.c
io_w.o: $(INCS) io_w.c
misc_w.o: $(INCS) misc_w.c
transfer_w.o: $(INCS) transfer_w.c

install: libpamrw.a 
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
	cp ../lib/*.a $(libdir)
	cp ../include/*.h $(includedir)


libpamrw.a: $(OBJECTS)
	$(AR) r libpamrw.a $(OBJECTS)
	/bin/cp libpamrw.a ../lib

clean:
	echo; echo "Cleaning up"
	rm *_.c > /dev/null 2>&1
	rm *.o > /dev/null 2>&1
	rm *~ > /dev/null
	rm libpamr.a > /dev/null 2>&1
	rm *.sdf > /dev/null 2>&1
	rm config.cache config.log config.status > /dev/null 2>&1
	rm work.pc work.pcl > /dev/null 2>&1
	rm Makefile > /dev/null 2>&1
