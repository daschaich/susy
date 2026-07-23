#!/bin/bash

# Script used for 2d runs on Param Smriti
# Original Script used at some other cluster https://github.com/daschaich/susy_scripts/blob/master/runLiv  

if [ $# -lt 10 ]; then
  echo "Usage: $0 <first> <last> <batch> <fsteps> <gsteps> <N> <rt> <gamma> <target> <time>"
  exit 1
fi


# Input parameters
first=$1
last=$2
batch=$3
fsteps=$4
gsteps=$5
N=$6
rt=$7
gamma=$8
target=$9
time=${10}

# Check whether we've correctly set $batch to evenly divide ($last-$first)
iter=0
count=$first
for(( i=$first ; $i<$last ; i+=$batch )); do
  iter=$[$iter + 1]
  count=$[$count + $batch]
done
echo "Will submit $iter jobs and end up at $count MDTU"

# Adjustable parameters
nodes=1
cpus=8
Nx=4
Nt=4
lambda=`echo $rt | awk -v nt="$Nt" '{print($1/nt)*($1/nt)}'`
bmass=`echo $lambda | awk -v gamma="$gamma" '{print(sqrt($1)*gamma)}'`
fmass=0.0
kappa=0.0
G=0.0
B=0.0
Ntraj=10
traj_length=1
skip=`echo $Ntraj | awk -v tau="$traj_length" '{print($1*tau)}'`

# Common parameters for all jobs

tag="N_${N}_${Nx}x${Nt}_rt_g_${rt}_${gamma}_${target}"

dir=/home/vamika/bana/susy/2d_Q16/susy/$tag
bin=/home/vamika/bana/susy/2d_Q16/susy/susy_$target

# Create the main directory and subdirectories
mkdir -p "$dir/Out" "$dir/Configs"

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

# Change N in 4dQ16
sed -i -E "s/(#define NCOL) .*/\1 $N/" /home/vamika/bana/susy/4d_Q16/include/susy.h
echo "Changed N in 4dQ16: as $N, pwd = $PWD"

# Change N in 2dQ16
sed -i -E "s/(#define NCOL) .*/\1 $N/" /home/vamika/bana/susy/2d_Q16/include/susy.h
echo "Changed N in 2dQ16: as $N, pwd = $PWD"

# Compile susy_hmc_pg in 2d (mandatory)
echo "Compiling susy_$target..."
if ! make -f Make_mpi susy_$target >& /dev/null ; then
  echo "ERROR: susy_$target compilation failed"
  make -f Make_mpi susy_$target
  exit
fi

cd $dir
 
echo "#!/bin/sh" > temp
echo "#SBATCH --partition=standard" >> temp
echo "#SBATCH  --ntasks=$cpus" >> temp
echo "#SBATCH --nodes=$nodes" >> temp
#echo "#SBATCH --nodelist=cn040,cn041" >> temp
#echo "#SBATCH --ntasks-per-node=36" >> temp
#echo "#SBATCH --mem=1G" >> temp
echo "#SBATCH --time=$time" >> temp
echo "#SBATCH -o job.%j.out" >> temp
echo "#SBATCH -e job.%j.err" >> temp
echo "#SBATCH -J sk_$tag" >> temp

echo "cd $dir" >> temp

# Check that we're not going to break anything,
# either through this job or the subsequent jobs it will submit
<<comment
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
    echo "ERROR: OUTPUT FILE $out EXISTS, file $out removed"
   # rm -f temp
   # exit 1
   rm -f $out
  fi
  if [ -f $lat ]; then
    echo "ERROR: LATTICE $lat EXISTS, file $lat removed"
   # rm -f temp
   # exit 1
   rm -f $lat
  fi
done

# Write this job's evolution tasks to run in a single job
iter=0
this_sub=$[$first + $batch]
for(( i=$first ; $i<$this_sub ; i+=$skip )); do
  iter=$[$iter + 1]
  next=$[$i + $skip]
  out=$dir/Out/out.$i-$next
  lat=$dir/Configs/gauge

  echo "echo \"Job HMC_${Nx}nt${Nt}_$tag started \"\`date\`\" jobid \$SLURM_JOBID\" >> $out" >> temp
  echo "echo \"=== Running MPI application on $cpus cpus ===\" >> $out" >> temp
  echo "echo \"mpirun -quiet -np $cpus $bin\" >> $out" >> temp
  echo "mpirun -quiet -np $cpus $bin << EOF >> $out" >> temp
  echo "prompt 0" >> temp
  echo "nx $Nx" >> temp
  echo "nt $Nt" >> temp
  echo "PBC -1" >> temp
  echo "iseed 41" >> temp

  echo "Nroot 1" >> temp
  echo "Norder 15" >> temp

  echo "warms 0" >> temp
  echo "trajecs $Ntraj" >> temp
  echo "traj_length $traj_length" >> temp
  echo "nstep $fsteps" >> temp
  echo "nstep_gauge $gsteps" >> temp
  echo "traj_between_meas $Ntraj" >> temp

  echo "lambda $lambda" >> temp
  echo "kappa_u1 $kappa" >> temp
  echo "bmass $bmass" >> temp
  echo "fmass $fmass" >> temp
  echo "G $G" >> temp

  echo "max_cg_iterations 5000" >> temp
  echo "error_per_site 1e-05" >> temp

  if $fresh_uncommented; then
    echo "fresh" >> temp
  else
    echo "reload_serial $lat.$i" >> temp
  fi
  echo "save_serial $lat.$next" >> temp
  echo "EOF" >> temp

  echo "echo \"=== MPI application finished at \"\`date\`\" ===\" >> $out" >> temp
  echo "" >> temp
done

# Submit next job, if applicable
# Warned above about possibility of ending up between $last and $last+$batch
if [ $this_sub -lt $last ] ; then
  echo "echo \"./2dsusy.sh $this_sub $last $batch $fsteps $gsteps $time\"" >> temp
  echo "./2dsusy.sh $this_sub $last $batch $fsteps $gsteps $time" >> temp
fi

sbatch temp
rm -f temp
echo "Requested $time to save $iter configs ($first--$this_sub by $skip)"
