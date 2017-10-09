#! /bin/sh

P=`basename $0`

MPIRUN=c2run

usage() {
cat<<END
usage: $P <nprocs>

       Runs PAMR wave example on <nprocs> processors
   
       Must be invoked from build directory, which in turn must be 
       globally visible from nodes (i.e. within /d/vnfe[123]/{/home,/home2})

       Uses $MPIRUN, which automatically selects 
       idlest nodes, and starts 'top' on first and last machines in the 
       machine file
END
exit 1
}

die() {
   echo "$P: $1"
   exit 1
}

case $# in
1) NPROCS=$1;;
*) usage;;
esac

test -x wave || die "Executable 'wave' not found."
test -d run_2d || die "Working directory 'run_2d' not found."
cd run_2d
cat ../wave.fparam wave.rtparam > wave.param
$MPIRUN -n $NPROCS ../wave wave.param
