include("/global/cfs/cdirs/m4265/cfundell/JuliaChem.jl/example_scripts/jc_timings_read.jl") # JuliaChem/example_scripts/jc_timings_read.jl
include("/global/cfs/cdirs/m4265/cfundell/JuliaChem.jl/example_scripts/Scripts/input_output_file_helper.jl")
using HDF5
 
function get_run_scf_energies(path_to_outputs)
    total_SCF_Energies = Dict{String, Float64}() # Energies for the various input file
   
    #files at path
    run_folders, full_folder_paths, full_folder_paths = get_output_folders(path_to_outputs)
    for run_folder_path in full_folder_paths
        files = readdir(run_folder_path)
        full_file_paths = [joinpath(run_folder_path, file) for file in files]
 
        run_name = split(run_folder_path, "/")[end]
        input_SCF_Energies = []
 
        for file_path in full_file_paths #copies for each Run and MPI RANK?
            run_level_data, timings, scf_options, scf_options_user, non_timing_data = jc_timings_read(file_path)
            #create dictionary from Nx2 array
            run_level_data_dict = Dict(run_level_data[:, 1] .=> run_level_data[:, 2])
            push!(input_SCF_Energies, parse(Float64, run_level_data_dict["scf_energy"])) # Assuming the first row contains the SCF energy
        end
        total_SCF_Energies[run_name] = input_SCF_Energies[1] # Store the first SCF energy value for each run configuration
    end
    return total_SCF_Energies
end
 
 
 
function main(paths_to_outputs)
   
    scf_energies_for_run_type = Array{Dict{String, Float64}}(undef, length(paths_to_outputs))
 
    if length(paths_to_outputs) > 2
        error("This script is designed to handle only two paths to outputs at the moment.")
    end
 
    i = 1
    for path in paths_to_outputs
        println("path: $path")
        folder_scf_energies = get_run_scf_energies(path)
        scf_energies_for_run_type[i] = folder_scf_energies
        i += 1
    end
 
    # Calculate the minimum, maximum, and average E differences
    max_diff = -100000.0
    max_diff_input = ""
    min_diff = 100000.0
    min_diff_input = ""
    difference_total = 0
    number_of_inputs = 0
    input_file_keys = keys(scf_energies_for_run_type[1])
    for input_key in input_file_keys
        scf_energy_1 = scf_energies_for_run_type[1][input_key]
        scf_energy_2 = scf_energies_for_run_type[2][input_key]
        diff = abs(scf_energy_1 - scf_energy_2)
        difference_total += diff
        println("SCF ΔE for $input_key: | $(scf_energy_1) - $(scf_energy_2) | = $diff")
        # check if new maximum
        if diff > max_diff
            max_diff = diff
            max_diff_input = input_key
        end
        # check if new minimum 
        if diff < min_diff
            min_diff = diff
            min_diff_input = input_key
        end        

        number_of_inputs += 1
    end
    println("Maximum difference in SCF energies: $max_diff")
    println("Minimum difference in SCF energies: $min_diff") 

    average_diff = difference_total / number_of_inputs
    println("Average difference in SCF energies: $average_diff")

    

    # display(scf_energies_for_run_type)
end
 
# path1 =  "/pscratch/sd/j/jhayes1/source/JuliaChem.jl/JuliaChem-Papers/DF-RHF-Paper/Benchmarks/0.4.3/1hsg/perlmutter/DF_RHF_screenedCPU"
# path2 =  "/pscratch/sd/j/jhayes1/source/JuliaChem.jl/JuliaChem-Papers/DF-RHF-Paper/Benchmarks/0.4.3/1hsg/perlmutter/DF_RHF_screenedCPU_no_sym"
# main([path1, path2])

path1 = "/global/cfs/cdirs/m4265/cfundell/JuliaChem_runs_only/JuliaChem.jl/Benchmarks/scripts/mixed_precision/DF_RHF_screenedCPU_6-311++G_2d_2p_"
path2 = "/global/cfs/cdirs/m4265/mixed_precision/JuliaChem.jl/Benchmarks/scripts/mixed_precision/DF_RHF_screenedCPU_mixed_6-311++G_2d_2p__S22_divide10"

main([path1, path2]) 
