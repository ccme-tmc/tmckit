#!/bin/bash  -f 


echo ' ---------------------------------------------------------'
echo ' *                      qw2k.sh                      *'
echo ' ---------------------------------------------------------'
myname=`basename $0`
myname_full=$0 

if [ $# -eq "0" ]  # Script invoked with no command-line args?
then
  cat <<EOF
Usage: `basename $0` [options]  <job_name [job_options]>
  Options:
   -a <addon>        # a string of commands (in c-shell format) that is inserted into qsub script directly:
                     # this can be used, for example, set additional environmental variable
   -e                # an "empty" run: generates relevant script, but not really run 
   -h <hours=960>    # requested walltime in unit of hour
   -i                # run the job interactively
   -m <mpi=1>        # the number of precesses used for MPI parallelization
   -n <nnodes=1>     # number of nodes
   -p <npnode=1>     # number of processes per node
   -q <queue>        # choose which queue to submit the job, if not set, it is automatically set in terms of nproc
                     #  if set "local", run the job locally
   -s <scrdir>       # set scratch directroy, if not set, it will use environmental variable SCRATCH
                     # this is very useful if one needs to run several wien2k jobs in different directory but with same case name 
   -r <additional resources> # c.f. qsub 
   -t <thread=1>     # the number of threads (default 1, used when using OpenMP)
   -v <w2kver=091>   # wien2k version (073/091/091.1/10)
   -N <jobname >     # a informative name for the job 

Examples:
   A typical use is like:
     qw2k.sh -N FePO4-lda-u04  -p 8 runsp_lapw -p -ec 0.000001 -orb
  which means it will use 8-processes k-point parallelization. 
     qw2k.sh -N FePO4-lda-u04  -p 8 -m 2 runsp_lapw -p -ec 0.000001 -orb
  which means totally 8 processes are used, with 2 processes are used in MPI parallelization, the   
  i.e., all processes are divided into 4 groups, two processes in each group take care of a sub-set of k-points
  Please be sure that the number of processed used for k-point parallelization, nproc/nmpi, should never larger 
  than the number of irreducible k-points available in actually calulations.      

Please be sure that the number of 

  Q: how to run a job with string option like min_lapw 
      min_lapw -j “runsp_lapw -p -ec 0.000001 -orb”
  A: put the command line into a file named, say, “min.sh”, then chmod +x min.sh, and then 
      qw2k.sh -N FePO4-min-lda-u04 -p 8 ./min.sh
EOF
  exit 1
fi

WIENDIR=/export/home/tmc/programs/wien2k
mpi0=1
interactive=0
emptyrun=0

while getopts "a:d:eh:im:n:o:p:q:r:s:t:v:N:" Option
do
  case $Option in
    a     ) addon="$OPTARG";; 
    d     ) debig=$OPTARG;; 
    e     ) emptyrun=1;;
    h     ) hours=$OPTARG;;
    i     ) interactive=1;;
    m     ) mpi=$OPTARG;; 
    n     ) nnodes=$OPTARG;;
    o     ) output=$OPTARG;;
    p     ) npnode=$OPTARG;;
    q     ) queue=$OPTARG;;
    r     ) resources=$OPTARG;;
    s     ) scrdir=$OPTARG;;
    t     ) thread=$OPTARG;;
    v     ) w2kver=$OPTARG;;
    N     ) jobname=$OPTARG;;
    *     ) echo "Unimplemented option chosen.";;
  esac
done
shift $(($OPTIND - 1))

w2kver=${w2kver:=091}
WIENROOT=$WIENDIR/v$w2kver
export PATH=$WIENROOT:$PATH
exec_basename=`basename $1`
exec_fullname=`which $1`

workdir=`pwd`
dirname=`basename $workdir`		
date=`date  +%m-%d-%H`
hours=${hours:=960}
scrdir=${scrdir:=$SCRATCH}
npnode=${npnode:=1}
nnodes=${nnodes:=1}
resources=${resources:=''}
mpi=${mpi:=1}
thread=${thread:=1}
jobname=${jobname:=$dirname-$date}
output=${output:=$jobname.log}
nproc=$(( $npnode*$nnodes ))
subsh="$scrdir/q-$jobname-$$.sh"
runjob="$exec_fullname $* "

# set default queue in terms of nproc
if test $nproc -gt 1 
then 
  def_queue=high
else 
  def_queue=low
fi

queue=${queue:=$def_queue}

if [ $interactive -eq 1 ]; then  
  queue="local"
fi 

if [ "$queue" == "local" ]; then 
  echo "nproc=" $nproc
  echo "mpi="   $mpi
  echo "SCRATCH=" $SCRATCH

  export SCRATCH=$scrdir

  if [ $nproc -gt 1 ]; then   
    host=`hostname `
    echo "Create the .machines file using the local host " $host

    echo "# .machines created by $myname" > .machines

    if [ $mpi0 -gt 0 ]; then  
      echo -n 'lapw0:' >> .machines
      i=1
      while [ $i -le $nproc ]
      do 
        echo -n "$host "  >>.machines
        (( i = i + mpi0 ))
      done
      echo '' >> .machines
    fi 

    #example for k-point and mpi parallel lapw1/2
    i=1
    while [ $i -le $nproc ]
    do
      echo -n '1:' >>.machines
      (( i1 = i + mpi ))
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

  if [ $emptyrun -eq 1 ]  
  then 
    echo "an empty run: stop here !"
    exit 0
  fi 

  export OMP_NUM_THREADS=$thread
  if [ $interactive -eq 1 ]; then  
    $runjob
  else 
    nohup $runjob >& $output &
  fi 
  exit $status
fi 

echo " Program to be run:" $exec_fullname
echo " Working directory:" $workdir
echo " scratch spacke   :" $scrdir
echo " nproc required   :" $nproc
echo " submission script:" $subsh 
echo ''

cat <<EOF >$subsh 
#!/bin/tcsh -f
#PBS -o q-$output
#PBS -e q-$jobname.err
#PBS -l ${resources}nodes=$nnodes:ppn=$npnode,walltime=$hours:00:00
#PBS -q batch
#PBS -m ae
#PBS -M $MAIL

setenv PATH $WIENROOT:\$PATH
setenv OMP_NUM_THREADS $thread
setenv SCRATCH $scrdir
cd $workdir 
$addon
echo 'Running WIEN2k programs using ' \`which x_lapw\` 

# setup wien2k parallel
# get processes list
set proclist=\`cat \$PBS_NODEFILE\`
set nproc=\$#proclist
echo "proclist:" \$proclist
echo "number of processors: \$nproc"

if ( \$nproc > 1 ) then 
  echo '#' > .machines

  if ( $mpi0 > 0 ) then 
    echo -n 'lapw0:' >> .machines
    set i=1
    while (\$i <= \$nproc )
      echo -n "\$proclist[\$i] " >>.machines
      @ i = \$i + $mpi0
    end
    echo ' ' >>.machines
  endif

  #example for k-point and mpi parallel lapw1/2
  set i=1
  while (\$i <= \$nproc )
    echo -n '1:' >>.machines
    @ i1 = \$i + $mpi
    @ i2 = \$i1 - 1
    echo \$proclist[\$i-\$i2] >>.machines
    set i=\$i1
  end
  echo '' >>.machines
  echo 'granularity:1' >>.machines
  echo 'extrafine:1' >>.machines
#  echo 'lapw2_vector_split:1 ' >>.machines

  # --------- end of .machines file

  # --------- writing  .processes and .machine1 file (only needed in case you start with lapw2)
  echo  -n 'init:' > .processes
  echo \$proclist >> .processes
  echo '1 : ' \$proclist[1] " :  1 : $mpi : 1" >> .processes

  if ( -e .machine1 ) rm .machine1
  set i=1
  while (\$i <= \$nproc )
    echo \$proclist[\$i] >>.machine1
    @ i = \$i + 1
  end

endif 
# -------------  end of setup --------------

$runjob
EOF

chmod +x $subsh

if [ $emptyrun -eq "1" ]
then
  echo "an empty run: stop here !"
  exit 0
fi

qsub -q $queue -N $jobname $subsh 

echo "a job has been submitted ==> " $runjob
exit 0
