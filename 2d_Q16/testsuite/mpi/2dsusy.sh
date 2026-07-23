#!/bin/bash
#SBATCH --mem=4G  # or any appropriate value
#./2dsusy.sh 0 10 10 10 N rt zeta target?
#./2dsusy.sh 10 40 10 10 ?
if [ $# -lt 8 ]; then
  echo "Usage: $0 <first> <end> nsteps:{<fermion> <gauge>} "
  exit 1
fi

# Input parameters
first=$1
last=$2
fsteps=$3
gsteps=$4
N=$5 # Number of colors
rt=$6
zeta=$7 
target=$8

nodes=1
cpus=8
nx=4
nt=4
tag="N_${N}_${nx}x${nt}_rt_g_${rt}_${zeta}_${target}"
#subtag="nstep_${fsteps}"

# Paths
dir=/home/bana/susy/2d_Q16/susy/$tag  # /$subtag
bin=/home/bana/susy/2d_Q16/susy/susy_$target

# Create the main directory and subdirectories
mkdir -p "$dir/Out" "$dir/Configs"

out_dir=$dir/Out
config_dir=$dir/Configs


out_dir=$dir/Out
config_dir=$dir/Configs

# Handle actions based on the value of 'first'
if [ "$first" -eq 0 ]; then
  # Remove all previous files in Out and Configs directories
  echo "Cleaning up previous files in $out_dir and $config_dir..."
  rm -rf "$out_dir"/* "$config_dir"/*

  # Ensure line "reload_serial $lat.$i" is commented and "fresh" is uncommented
  fresh_uncommented=true
else  
#elif [ "$first" -ne 10 ]; then
  # Ensure line "fresh" is commented and "reload_serial $lat.$i" is uncommented
  fresh_uncommented=false
#else
#  echo "Invalid 'first' value. Only 0 or 10 is supported."
#  exit 1
fi


  # Change N in 4d
  #cd ../../../4d_Q16/susy/
  sed -i -E "s/(#define NCOL) .*/\1 $N/" /home/bana/susy/4d_Q16/include/susy.h
  echo "Changed N in 4d: $PWD"
  
  # Compile in 4d
  # echo "Compiling susy_$target..."
  # if ! make -f Make_scalar susy_$target >& /dev/null ; then
  #   echo "ERROR: susy_$target compilation failed"
  #   make -f Make_scalar susy_$target
  #   exit
  # fi
  

  # Change N in 2d
  #cd ../../susy/
  sed -i -E "s/(#define NCOL) .*/\1 $N/" /home/bana/susy/2d_Q16/include/susy.h
  echo "Changed N in 2d: $PWD"

  #cd to susy_hmc_pg from mpi

  cd ../../susy/
  echo "pwd : $PWD"
  # Compile susy_hmc_pg in 2d (mandatory)
  echo "Compiling susy_$target..."
  if ! make -f Make_mpi susy_$target >& /dev/null ; then
    echo "ERROR: susy_$target compilation failed"
    make -f Make_mpi susy_$target
    exit
  fi

  # Change directory to testsuite
  cd ../testsuite/mpi
  echo "Current directory: $PWD"



echo "Will submit $iter jobs and end up at $count MDTU"

# Adjustable parameters
#nodes=1
#cpus=6
#nx=6
#nt=6
#lambda=1.0
lambda=$(echo "scale=8; ($rt / $nt)^2" | bc)
echo "lambda = $lambda"
#bmass=0.1
bmass=$(echo "scale=8; $zeta * ($rt / $nt)" | bc)
echo "bmass = $bmass"
fmass=0.0
kappa_u1=0.0
G=0.0
Ntraj=10
traj_length=1
skip=10
#skip=`echo $Ntraj | awk -v tau="$traj_length" '{print($1*tau)}'` # skip = Ntraj*traj_length
# Common parameters for all jobs





cd $dir
echo "#!/bin/sh" > temp
echo "cd $dir" >> temp


# Check that we're not going to break anything,
# either through this job or the subsequent jobs it will submit
<<comment # and refresh
lat=$dir/Configs/gauge.$first
if [ ! -f $lat ]; then
  echo "ERROR: LATTICE $lat NOT FOUND, SUBMISSION ABORTED"
  rm -f temp
  exit 1
fi
comment
for(( i=$first ; $i<$last ; i+=$skip )); do
  next=$[$i + $skip]
  out=$dir/Out/out.$i-$next
  lat=$dir/Configs/gauge.$next
  if [ -f $out ]; then
    echo "ERROR: OUTPUT FILE $out EXISTS, but SUBMISSION not ABORTED"
   # rm -f temp
   # exit 1
  fi
  if [ -f $lat ]; then
    echo "ERROR: LATTICE $lat EXISTS, but SUBMISSION not ABORTED"
   # rm -f temp
   # exit 1
  fi
done

# Write this job's evolution tasks to run in a single job
iter=0
for(( i=$first ; $i<$last ; i+=$skip )); do
  iter=$[$iter + 1]
  next=$[$i + $skip]
  out=$dir/Out/out.$i-$next
  lat=$dir/Configs/gauge

  echo "echo \"Job HMC_2d temp run started \"\`date\`\" jobid \$SLURM_JOBID\" >> $out" >> temp
  echo "echo \"=== Running MPI application on $cpus cpus ===\" >> $out" >> temp
  echo "echo \"mpirun -np $cpus $bin\" >> $out" >> temp
  echo "mpirun -n $cpus $bin << EOF >> $out" >> temp
  echo "prompt 0" >> temp
  echo "nx $nx" >> temp
  echo "nt $nt" >> temp
  echo "PBC -1" >> temp
  #echo "iseed ${last}41$i" >> temp
  echo "iseed 41" >> temp

  echo "Nroot 1" >> temp
  echo "Norder 15" >> temp

  echo "warms 0" >> temp
  echo "trajecs $Ntraj" >> temp
  echo "traj_length $traj_length" >> temp
  echo "nstep $fsteps" >> temp # $3 step size = traj_length/ nsteps = 1/10= 0.1 , nsteps = fermion  field steps
  echo "nstep_gauge $gsteps" >> temp # $4 gauge steps / fermion steps 
  echo "traj_between_meas $Ntraj" >> temp

  echo "lambda $lambda" >> temp
  echo "kappa_u1 $kappa_u1" >> temp
  echo "bmass $bmass" >> temp
  echo "fmass $fmass" >> temp
  echo "G $G" >> temp

  echo "max_cg_iterations 15000" >> temp
  echo "error_per_site 1e-7" >> temp
  
  if $fresh_uncommented; then
    echo "fresh" >> temp
  else
    echo "reload_serial $lat.$i" >> temp
  fi

  #echo "fresh" >> temp  # with comment out lat above 
  #echo "reload_serial $lat.$i" >> temp
  echo "save_serial $lat.$next" >> temp
  echo "EOF" >> temp

  echo "echo \"=== MPI application finished at \"\`date\`\" ===\" >> $out" >> temp
  echo "" >> temp
done

bash temp
rm -f temp



