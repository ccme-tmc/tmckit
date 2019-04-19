#!/bin/bash  -f 

# this script is used to submit jobs to the queuing system on 
# shenteng.sccas.cn
#
# the available options are explained by 
#          qrun-b.sh 
# without any command-line options
echo ' ---------------------------------------------------------'
echo ' *                      qrun-b.sh                        *'
echo ' ---------------------------------------------------------'

myname=`basename $0`
myname_full=`which $myname`

if [ $# -eq "0" ]  # Script invoked with no command-line args?
then
  cat <<EOF
Usage: $myname [options] <job to be submitted and its command lines options> 
  Options:
   -e              # empty run, stop before the script is submitted (for debug)  
   -h <cpuhours>   # required cpu hours 
   -i <input_file> # if the command to be executed needs an input file in the formmat 
                   #  job.x < input_file
   -m <0/1/2/3>  # 0 -- not using mpi; 
                   1 -- intelmpi
                   2 -- mvapich
                   3 -- openmpi 
   -n <nnodes>        # number of nodes
   -o <output_file>   # file name for the standard output  
   -p <npnode>        # number of processes per node, nproc=nnodes*npnode
   -q <queue>         # choose which queue to submit the job, if not set, it will be 
                      # determined automatically in terms of nproc, if set as "local" 
                      # the job will be run on the local computer  
   -t <thread>        # the number of threads (default 1, used when using OpenMP)
   -v <name=value>    # set up additional environmental values (this can appear several times) 
   -N <jobname >      # job name 
Examples:

  1. run GW code
   $myname -p 8 -N Si-gw-e10-k64 gwsv6-mpi
   $myname -p 8 -n 2 -a “export MPI_NPROC_COL=8” -N Si-gw-e10-k64  gwsv6-mpi 

  2. run g09: 
  g09 does not use MPI, but you can accelerate the calculation using multi-threading parallelation
   $myname -i test.gjf -o test.log -t 4 -p 4 -m 0 -N CeO2 g09       # if you have %nproc=8 in test.gjf 
EOF
  exit 
fi

#MPIRUN=mpirun
env_setup=''
emptyrun=0
while getopts "ei:h:m:n:N:o:p:q:t:v:" Option
do
  case $Option in
    e     ) emptyrun=1;;
    i     ) input=$OPTARG;;
    h     ) hours=$OPTARG;;
    m     ) mpi=$OPTARG;;
    n     ) nnodes=$OPTARG;;
    N     ) jobname=$OPTARG;;
    o     ) output=$OPTARG;;
    p     ) npnode=$OPTARG;;
    q     ) queue=$OPTARG;;
    t     ) thread=$OPTARG;;
    v     ) env_setup="export $OPTARG; $env_setup";;
    *     ) echo "Unimplemented option $Option chosen "; exit 1 ;;
  esac
done
shift $(($OPTIND - 1))
exec_basename=`basename $1`
exec_fullname=`which $1`
shift

workdir=`pwd`
dirname=`basename $workdir`		
date=`date  +%m-%d-%H`

npnode=${npnode:=1}
nnodes=${nnodes:=1}
mpi=${mpi:=1}
hours=${hours:=360}
resources=${resources:=""}
thread=${thread:=1}
nproc=$(( $npnode*$nnodes ))
queue=${queue:=x64_tmc}
jobname=${jobname:=$dirname-$date}
output=${output:=$jobname.log}

if [ $mpi -eq 1 ]; then 
  echo "Run the MPI job by intelmpi" 
  MPIRUN=mpijob.lsf
elif [ $mpi -eq 2 ]; then 
  echo "Run the job by mvaipch" 
  MPIRUN=mpijob.lsf 
elif [ $mpi -eq 3 ]; then 
  echo "Run the MPI job by openmpi" 
  MPIRUN=mpijob.openmpi
else
  echo "Run the job without using MPI"
  MPIRUN=""
fi 


if [ $input ];  then
  echo "input=" $input
  runjob="$MPIRUN $exec_fullname < $input >& $output"
else 
  runjob="$MPIRUN $exec_fullname $* >& $output"
fi

echo " The job to be run:" $runjob
echo " Working directory:" $workdir
echo " nproc required   :" $nproc
echo " Output file name :" $output 
echo ''

if [ "$queue" == "local" ]; then
  echo "Run the job on the local computer:"
  nohup $runjob &
  exit 0
fi    

TMPDIR=${TMPDIR:=.}
subsh="$TMPDIR/q-$jobname-$$.sh"
echo ' submission script:' $subsh

cat <<EOF >$subsh 
#!/bin/sh -f
#BSUB -n $nproc
#BSUB -R "span[ptile=$npnode]"
#BSUB -o q-$jobname.out
#BSUB -e q-$jobname.err
#BSUB -W $hours:00
#BSUB -a $mpitype
#BSUB -q $queue

export OMP_NUM_THREADS=$thread
$env_setup 
cd $workdir 
$runjob
EOF

if [ $emptyrun -eq 1 ]; then
  echo "an empty run: stop here !"
  exit 0
fi

bsub < $subsh 
