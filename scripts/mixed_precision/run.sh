
# # setup Julia 
# # e.g. export PATH=/PATH/TO/Julia/bin/:$PATH
# # or module load julia and other modules as needed
# # this script assumes that the JuliaChem package is installed in the project path and the LIBINT wrapper is compiled 

# module load julia/1.11.4

# setup Julia 
module load julia/1.11.4

if [ -z "$JC_PATH" ]; then
    echo "Error: JC_PATH environment variable is not set."
    echo "Please set JC_PATH to the path of your JuliaChem installation."
    exit 1
fi

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
export run_start_index=1
export number_of_runs=2
export file_start_index=22
export file_end_index=22

# # options are divide_total_aux e.g. 1,2,4,8,9,10 , range_of_range_counts e.g. 1:10, list_of_range_counts e.g. 1,2,5,10 
# export mixed_df_Q_range_mode="divide_total_aux"
# export mixed_df_Q_ranges="10"
# options for mixed precision

export mixed_df_Q_range_mode="list_of_range_counts"
export mixed_df_Q_ranges="10"

run_julia() {
    local settings_id=$1
    local contraction_float_type=$2
    local use_sym=$3
    local mixed_df_Q_range_mode=$4
    local mixed_df_Q_ranges=$5

    export settings_id  
    export contraction_float_type
    export use_sym
    export mixed_df_Q_range_mode
    export mixed_df_Q_ranges

    export output_full_path=$script_path/${settings_id}_${basis_id}_${contraction_float_type}_use_sym_${use_sym}
    mkdir -p $output_full_path
    # srun -N 1 -n 1  --cpu-bind=socket --cpu-bind=v julia \
    julia \
     --project=$j_project_path --threads=$threads_per_socket \
     $run_script $input_file_path $output_full_path $settings_id $basis_id \
      $run_start_index $number_of_runs $file_start_index $file_end_index \
      &> $output_full_path/output.log
}

# DF_RHF_screenedCPU_mixed
run_julia DF_RHF_screenedCPU_mixed single true list_of_range_counts 4
run_julia DF_RHF_screenedCPU_mixed double true list_of_range_counts 2
# run_julia DF_RHF_screenedCPU_mixed single false list_of_range_counts 4
# run_julia DF_RHF_screenedCPU_mixed double false list_of_range_counts 4

# DF_RHF_denseCPU_mixed
# run_julia DF_RHF_denseCPU_mixed single true list_of_range_counts 4
# run_julia DF_RHF_denseCPU_mixed double true list_of_range_counts 4
# run_julia DF_RHF_denseCPU_mixed single false list_of_range_counts 4
# run_julia DF_RHF_denseCPU_mixed double false list_of_range_counts 4

# DF_RHF_denseCPU (no contraction_float_type, just use last set value)
# run_julia DF_RHF_denseCPU double false list_of_range_counts 4
# run_julia DF_RHF_denseCPU double true list_of_range_counts 4

# DF_RHF_screenedCPU (no contraction_float_type, just use last set value)
run_julia DF_RHF_screenedCPU double true list_of_range_counts 1
# run_julia DF_RHF_screenedCPU double false list_of_range_counts 1

# test refactor full precision vs current dev branch full precision 