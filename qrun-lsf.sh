#!/bin/bash  -f 

# this script is used to submit jobs to the queuing system on Beijing Computing Center(BCC)
# the available options are explained by 
#          qrun-lsf.sh 
# without any command-line options
echo ' ---------------------------------------------------------'
echo ' *                      qrun-lsf.sh                      *'
echo ' ---------------------------------------------------------'

myname=`basename $0`
myname_full=`which $myname`

if [ $# -eq "0" ]  # Script invoked with no command-line args?
then
  cat <<EOF
Usage: $myname [options] <job to be submitted and its command lines options> 
  Options:
   -e                 # empty run 
   -j <job>          # the job to be run 
   -m <mpi>          # 0 -- non-MPI jobs, 1 -- openmpi job
   -n <npnode>       # the number of processes per node (useful to reduce inter-node communication)  
   -p <nproc>        # number of processes
   -q <queue>        # choose which queue to submit the job 
   -s <scrdir>       # set the scratch directory
   -t <thread>       # the number of threads (default 1, used when using OpenMP)
   -w <w2k>          # option for wien2k parallelization  
                        -1 -- not a wien2k job 
                         0 -- wien2k with k-points parallelization only
                        >0 -- wien2k job with <w2k> processes for MPI 
   -N <jobname >     # job name 
Examples:

  1. run GW code
   $myname -p 8 -n 8 -N Si-gw-e10-k64 -m 1 gwsv6-mpi

  2. run g09: 
  g09 does not use MPI, but you can accelerate the calculation using multi-threading parallelation
   $myname -t 4 -p 4 -m 0 -j 'g09 < test.gjf'     # if you have %nproc=8 in test.gjf 
EOF
  exit 
fi

#MPIRUN=mpirun
env_setup=''
emptyrun=0
while getopts "ej:m:n:N:p:q:s:t:w:" Option
do
  case $Option in
    j     ) job=$OPTARG;;
    m     ) mpi=$OPTARG;;
    n     ) npnode=$OPTARG;;
    N     ) jobname=$OPTARG;;
    p     ) nproc=$OPTARG;;
    q     ) queue=$OPTARG;;
    s     ) scrdir=$OPTARG;;
    t     ) thread=$OPTARG;;
    w     ) w2k=$OPTARG;;
    *     ) echo "Unimplemented option $Option chosen "; exit 1 ;;
  esac
done
shift $(($OPTIND - 1))
job=${job:=$*}

workdir=`pwd`
dirname=`basename $workdir`		
date=`date  +%m-%d`

nproc=${nproc:=1}
npnode=${npnode:=0}
mpi=${mpi:=0}
thread=${thread:=1}
queue=${queue:=scivip}
scrdir=${scrdir:=$SCRATCH}
w2k=${w2k:=-1}
jobname=${jobname:=$dirname-$date}
init_job=''
if [ $w2k -ge 0 ]; then 
  init_job="w2k_initpara -mpi $w2k -q lsf"
fi 

if [ $mpi -eq 1 ]; then 
  echo "Run the MPI job by openmpi" 
  MPIRUN="mpirun -np $nproc"
else
  echo "Run the job without using MPI"
  MPIRUN=""
fi 

npnode_info=''
if [ $npnode -gt 0 ]; then 
  npnode_info="NP_PER_NODE=$npnode"
fi 

echo " The job to be run:" $job
echo " Working directory:" $workdir
echo " nproc required   :" $nproc
echo ''

subsh=q-$jobname.sh
echo ' submission script:' $subsh

cat <<EOF >$subsh 
APP_NAME=$queue
MY_MPI_TYPE=openmpi
MY_MPI_HOME=$MPI
OMP_NUM_THREADS=$thread
NP=$nproc
$npnode_info
RUN="RAW"
EOF

if [ $w2k -ge 0 ]; then
  if [ $w2k -eq 0 ]; then
    echo "w2k_initpara -mpi 1 -mpi0 0 -q lsf" >>$subsh
  else
    echo "w2k_initpara -mpi $w2k -mpi0 1 -q lsf" >>$subsh
  fi

  cat <<EOF >>$subsh 
# check whether scratch is accessable 
touch $scrdir/.tmp
if [ \$? -ne 0 ]  
then 
  export SCRATCH=./scr-$jobname 
else 
  export SCRATCH=$scrdir/scr-$jobname
fi 
mkdir -p \$SCRATCH
$job >& $jobname.log 
rm -rf \$SCRATCH
EOF

else
  echo "$MPIRUN $job >& $jobname.log" >> $subsh
fi 

chmod +x $subsh

if [ $emptyrun -eq 1 ]; then
  echo "an empty run: stop here !"
  exit 0
fi

bsub -J $jobname $workdir/$subsh 
