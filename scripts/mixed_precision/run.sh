
# setup Julia 
# e.g. export PATH=/PATH/TO/Julia/bin/:$PATH
# or module load julia and other modules as needed
# this script assumes that the JuliaChem package is installed in the project path and the LIBINT wrapper is compiled 

module load julia/1.11.4

# export PATH=/pscratch/sd/j/jhayes1/software/julia-1.11.4/bin/:$PATH

export JC_PATH=/global/cfs/cdirs/m4265/mixed_precision/JuliaChem.jl

export script_path=$JC_PATH/Benchmarks/scripts/mixed_precision

cd $script_path
export run_script=$script_path/run.jl
export JULIACHEM_SYSIMG_PATH=$JC_PATH/perlmutter_JC_sysimg.so

export j_project_path=$JC_PATH/mixed_precision_env_perl


export threads_per_socket=64
export MKL_NUM_THREADS=$threads_per_socket
export OPENBLAS_NUM_THREADS=$threads_per_socket
export JULIA_NUM_THREADS=$threads_per_socket

export input_file_path=$JC_PATH/example_inputs/S22


export basis_id="cc-pvdz-ri"
# export basis_id="cc-pvdz-jk"
# export basis_id="6-311++G_2d_2p_"
export run_start_index=1
export number_of_runs=1
export file_start_index=1
export file_end_index=5
export settings_id=DF_RHF_denseCPU

export output_full_path=$script_path/S22/${settings_id}_MIXED
mkdir -p $output_full_path

export DO_MIXED=true

# srun -n 1  --cpu-bind=socket --cpu-bind=v 
julia --sysimage=$JULIACHEM_SYSIMG_PATH \
 --project=$j_project_path --threads=$threads_per_socket --check-bounds=no --math-mode=fast --optimize=3 --inline=yes \
 $run_script $input_file_path $output_full_path $settings_id $basis_id \
  $run_start_index $number_of_runs $file_start_index $file_end_index \
  &>> $output_full_path/output.log

# export DO_MIXED=false
# export output_full_path=$script_path/S22/${settings_id}_DOUBLE

# julia --sysimage=$JULIACHEM_SYSIMG_PATH \
#  --project=$j_project_path --threads=$threads_per_socket --check-bounds=no --math-mode=fast --optimize=3 --inline=yes \
#  $run_script $input_file_path $output_full_path $settings_id $basis_id \
#   $run_start_index $number_of_runs $file_start_index $file_end_index \
#   &>> $output_full_path/output.log