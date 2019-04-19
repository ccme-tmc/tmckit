#!/bin/bash  -f 

# this script is used to submitted jobs to PBS-type queuing system 
# the available options are explained by 
#          qrun_pbs.sh 
# without any command-line options
echo ' ---------------------------------------------------------'
echo ' *                      qrun_pbs.sh                      *'
echo ' ---------------------------------------------------------'

OPENMPI_MPIRUN=mpirun
MPICH2_MPIEXEC=mpiexec

if [ $# -eq "0" ]  # Script invoked with no command-line args?
then
  cat <<EOF
Usage: `basename $0` [options] <job to be submitted and its command lines options> 
  Options:
   -a <addon>                   # a string of commands (in bash format) before running the job:
                                    this can be used, for example, set additional environmental variable
   -c <customize=/ /w2k/qe/>    # customize for a particular programs   
   -h <cpuhours>                # required cpu hours 
   -i <input_file=NA>           # provide the input file name, this is needed 
                                    if the command to be executed needs an input file in the formmat 
                                      job.x < input_file
   -j <'job'>                   # the job and related optins 
                                #  use this in case that the job to be run has confilciting command-line options as qrun.sh 
   -m <mpi=0>                   # the type of parallel job 
                                    0 -- not using mpi; 
                                    1 -- use openmpi  
   -n <nnodes=1>                # number of nodes to use 
   -o <output=NA>               # file name for the standard output  
   -p <npnode=1>                # number of processes per node, nproc=nnodes*npnode
   -q <queue>                   # choose which queue to submit the job, if not set, it will be 
                                    determined automatically in terms of nproc, if set as "local" 
                                    the job will be run on the local computer  
   -r <additional resource>     # c.f. qsub
   -s <scrdir>                  # set scratch directroy
   -t <thread>                  # the number of threads (default 1, used when using OpenMP) 
                                   in most cases, use the default value 1 
                                   (unless when running g09 with nproc option)
   -N <jobname >                # job name 
Examples:

  1. run GW code
   qrun.sh -p 8 -N Si-gw-e10-k64 gwsv6-mpi
   qrun.sh -p 8 -n 2 -a “export MPI_NPROC_COL=8” -N Si-gw-e10-k64  gwsv6-mpi 

  2. run g09: 
  g09 does not use MPI, but you can accelerate the calculation using multi-threading parallelation
   qrun.sh -i test.gjf -o test.log -t 4 -p 4 -m 0 -N CeO2 g09       # if you have %nproc=4 in test.gjf 
EOF
  exit 
fi

while getopts "i:a:h:j:m:n:o:p:q:r:s:t:N:" Option
do
  case $Option in
    a     ) addon="$OPTARG";; 
    i     ) input=$OPTARG;;
    h     ) hours=$OPTARG;;
    m     ) mpi=$OPTARG;;
    n     ) nnodes=$OPTARG;;
    o     ) output=$OPTARG;;
    p     ) npnode=$OPTARG;;
    q     ) queue=$OPTARG;;
    r     ) resources="$OPTARG,";;
    s     ) scrdir=$OPTARG;;
    t     ) thread=$OPTARG;;
    N     ) jobname=$OPTARG;;
    j     ) job="$OPTARG";;
    *     ) echo "Unimplemented option $Option chosen "; exit 1 ;;
  esac
done
shift $(($OPTIND - 1))
job=${job:=$*}
#exec_basename=`basename $job`
exec_fullname=$job

#shift


workdir=`pwd`
dirname=`basename $workdir`		
date=`date  +%m-%d-%H`

npnode=${npnode:=1}
nnodes=${nnodes:=1}
customize={customize:=''}
mpi=${mpi:=0}
hours=${hours:=960}
resources=${resources:=""}
thread=${thread:=1}
nproc=$(( $npnode*$nnodes ))

SCRATCH=${SCRATCH:=.}
scrdir=${scrdir:=$SCRATCH}
echo "input=" $input
jobname=${jobname:=$dirname-$date}

# set default queue 
if test $nproc -gt 1   
then 
  def_queue=high
else
  def_queue=low
fi 

queue=${queue:=$def_queue}

if [ "$customize" == "w2k" ] 
then 
  mpi=0
fi 

# set runjob 

if [ $mpi -eq 1 ]
then 
  echo "Run the job by MPI (openmpi)" 
  MPIRUN="$OPENMPI_MPIRUN  -np $nproc "
elif [ $mpi -eq 2 ]
then 
  echo "Run the job by MPI-II (mvaipch2)" 
  MPIRUN="$MPICH2_MPIEXEC  -n $nproc "
else
  echo "Run the job without using MPI"
  MPIRUN=""
fi 

if [ $nproc -eq 1 ]
then
  if [ $input ]
  then  
    echo "Input file:" $input 
    runjob="$exec_fullname < $input "
  else
    runjob="$exec_fullname "
  fi

else
  if [ $input ]
  then
    echo "Input file:" $input 
    runjob="$MPIRUN $exec_fullname < $input "
  else 
    runjob="$MPIRUN $exec_fullname "
  fi
fi

output=${output:=q-$jobname.out}
runjob="$runjob >& $output"

# check whether the SCRATCH space accessible
touch $scrdir/.tmp
if [ $? -ne 0 ] 
then 
  scrdir=.
fi 

echo " The job to be run:" $runjob
echo " Working directory:" $workdir
echo " scratch space    :" $scrdir
echo " nproc required   :" $nproc
echo ''

if [ "$queue" == "local" ]
then
  echo "Run the job on the local computer:"
  cd $workdir
  if [ "$customize" == "w2k" ] 
  then 
    if [ $nproc -gt 1 ]; then
      w2k_mpi0=${w2k_mpi0:=0}
      w2k_mpi=${w2k_mpi:=1}      
      host=`hostname `
      echo "Create the .machines file using the local host " $host
      echo "# .machines created by $myname" > .machines

      if [ $w2k_mpi0 -gt 0 ]; then
        echo -n 'lapw0:' >> .machines
        i=1
        while [ $i -le $nproc ]
        do
          echo -n "$host "  >>.machines
          (( i = i + w2k_mpi0 ))
        done
        echo '' >> .machines
      fi

      #example for k-point and mpi parallel lapw1/2
      i=1
      while [ $i -le $nproc ]
      do
        echo -n '1:' >>.machines
        (( i1 = i + w2k_mpi ))
        (( i2 = i1 - 1 ))
        j=$i
        while [ $j -le $i2 ]
        do
          echo -n "$host " >>.machines
          (( j = j + 1 ))
        done

        echo  '' >>.machines
        i=$i1
      done
      echo '' >>.machines
      echo 'granularity:1' >>.machines
      echo 'extrafine:1' >>.machines
      echo 'lapw2_vector_split:1 ' >>.machines

      cat .machines
    fi
  fi 

  export OMP_NUM_THREADS=$thread

  $runjob  &

  exit 0
fi    

TMPDIR=${TMPDIR:=$HOME/tmp/qtmp}; mkdir -p $TMPDIR
subsh="$TMPDIR/q-$jobname-$$.sh"
echo ' submission script:' $subsh

cat <<EOF >$subsh 
#!/bin/sh -f
#PBS -o $TMPDIR/q-$jobname.out
#PBS -e $TMPDIR/q-$jobname.err
#PBS -l ${resources}nodes=$nnodes:ppn=$npnode,walltime=$hours:00:00
#PBS -q $queue
#PBS -m ae
#PBS -M $MAIL
#PBS -N $jobname
#PBS -V

cd $workdir 
echo 'Job is running on node(s): '
cat \$PBS_NODEFILE
export OMP_NUM_THREADS=$thread
$addon   # some additional setup 
EOF

# some customized setup for some special programs
if [ "$customize" == 'w2k' ] 
then 
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

  echo 'Running WIEN2k programs using ' \`which x_lapw\`

  # some parameters needed to set up the machines
  w2k_mpi0=${w2k_mpi0:=0}
  w2k_mpi=${w2k_mpi:=1}      

  # setup wien2k parallel
  # get processes list
  proclist=\`cat \$PBS_NODEFILE\`
  nproc=\$#proclist
  echo "proclist:" \$proclist
  echo "number of processors: \$nproc"

  if [ \$nproc -gt 1 ]
  then
    echo '#' > .machines

    if [ \$w2k_mpi0 -gt 0 ]
    then
      echo -n 'lapw0:' >> .machines
      i=1
      while [ \$i -le \$nproc ]
      do
        echo -n "\$proclist[\$i] " >>.machines
        (( i = i + w2k_mpi0 ))
      done 
      echo ' ' >>.machines
    fi

    #example for k-point and mpi parallel lapw1/2
    i=1
    while [ \$i -le \$nproc ]
    do 
      echo -n '1:' >>.machines
      (( i1 = i + w2k_mpi ))
      (( i2 = i1 - 1  ))
      echo \$proclist[\$i-\$i2] >>.machines
      i=\$i1
    done 
  
    echo '' >>.machines
    echo 'granularity:1' >>.machines
    echo 'extrafine:1' >>.machines
    #  echo 'lapw2_vector_split:1 ' >>.machines

    # --------- end of .machines file

    # --------- writing  .processes and .machine1 file (only needed in case you start with lapw2)
    echo  -n 'init:' > .processes
    echo \$proclist >> .processes
    echo '1 : ' \$proclist[1] " :  1 : $w2k_mpi : 1" >> .processes

    if test -e .machine1 
      then  rm .machine1
    fi 

    i=1
    while [ \$i -le \$nproc ]
    do 
      echo \$proclist[\$i] >>.machine1
      (( i = i + 1 ))
    done 
  fi

EOF
fi

echo $runjob >> $subsh

if [ "$customize" == "w2k" ] 
then  
  echo "rm -rf \$SCRATCH" >> $subsh
fi 

chmod +x $subsh
qsub $subsh 
