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
 
    max_diff = -100000.0
    max_diff_input = ""
    input_file_keys = keys(scf_energies_for_run_type[1])
    for input_key in input_file_keys
        scf_energy_1 = scf_energies_for_run_type[1][input_key]
        scf_energy_2 = scf_energies_for_run_type[2][input_key]
        diff = abs(scf_energy_1 - scf_energy_2)
        println("SCF ΔE for $input_key: | $(scf_energy_1) - $(scf_energy_2) | = $diff")
        if diff > max_diff
            max_diff = diff
            max_diff_input = input_key
        end
    end
    println("Maximum difference in SCF energies: $max_diff")
 
 
   
 
    # display(scf_energies_for_run_type)
end
 
# path1 =  "/pscratch/sd/j/jhayes1/source/JuliaChem.jl/JuliaChem-Papers/DF-RHF-Paper/Benchmarks/0.4.3/1hsg/perlmutter/DF_RHF_screenedCPU"
# path2 =  "/pscratch/sd/j/jhayes1/source/JuliaChem.jl/JuliaChem-Papers/DF-RHF-Paper/Benchmarks/0.4.3/1hsg/perlmutter/DF_RHF_screenedCPU_no_sym"
# main([path1, path2])
 
path1="/global/cfs/cdirs/m4265/cfundell/JuliaChem.jl/Benchmarks/scripts/mixed_precision/S22/DF_RHF_denseCPU_DOUBLE_JK"
path2="/global/cfs/cdirs/m4265/cfundell/JuliaChem.jl/Benchmarks/scripts/mixed_precision/S22/RHF_staticCPU_DOUBLE_JK"
 
s22_path = "/global/cfs/cdirs/m4265/cfundell/JuliaChem.jl/Benchmarks/scripts/mixed_precision/S22"

#RHF vs DF 
println("6-311")
path1 = joinpath(s22_path, "RHF_staticCPU_DOUBLE_6-311")
path2 = joinpath(s22_path, "DF_RHF_denseCPU_DOUBLE_6-311")
println("RHF vs DF, ")
main([path1, path2])


#DF vs DF mixed
path1 = joinpath(s22_path, "DF_RHF_denseCPU_DOUBLE_6-311")
path2 = joinpath(s22_path, "DF_RHF_denseCPU_mixed_6-311")
println("DF vs DF mixed, 6-311")
main([path1, path2])

println("cc-pvdz JK")
println("RHF vs DF JK")

#DF vs DF mixed JK
path1 = joinpath(s22_path, "RHF_staticCPU_DOUBLE_JK")
path2 = joinpath(s22_path, "DF_RHF_denseCPU_DOUBLE_JK")

main([path1, path2])
println("DF vs DF mixed JK")
path1 = joinpath(s22_path, "DF_RHF_denseCPU_DOUBLE_JK")
path2 = joinpath(s22_path, "DF_RHF_denseCPU_mixed_JK")
main([path1, path2])

#RHF vs DF RI
path1 = joinpath(s22_path, "RHF_staticCPU_DOUBLE_RI")
path2 = joinpath(s22_path, "DF_RHF_denseCPU_DOUBLE_RI")
println("RHF vs DF RI")
main([path1, path2])
#DF vs DF mixed RI
path1 = joinpath(s22_path, "DF_RHF_denseCPU_DOUBLE_RI")
path2 = joinpath(s22_path, "DF_RHF_denseCPU_mixed_RI")
println("DF vs DF mixed RI")
main([path1, path2])