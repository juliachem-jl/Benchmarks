# using MKL
using BLISBLAS
println("starting run.jl");flush(stdout)
using JuliaChem
println("done loading JuliaChem");flush(stdout)
using Dates
using Pkg
using LinearAlgebra
using MPI 
using Base.Threads
using HDF5
# using JuliaChem.Properties

JC_PATH = ""
if haskey(ENV, "JC_PATH")
    JC_PATH = ENV["JC_PATH"]
else
    error("JC_PATH environment variable is not set. Please set it to the path of the JuliaChem repository.")
end
include("$(JC_PATH)/example_scripts/Scripts/input_output_file_helper.jl")
include("$(JC_PATH)/example_scripts/Settings/run_settings.jl")
include("$(JC_PATH)/example_scripts/Settings/get_basis_info.jl")
include("$(JC_PATH)/example_scripts/Scripts/rank_print.jl")
include("$(JC_PATH)/example_scripts/jc_timings_write.jl") # JuliaChem/example_scripts/jc_timings_write.jl
include("$(JC_PATH)/example_scripts/Scripts/replace_keywords.jl")



function run_file(input_file_path, input_file_name, output_dir, scf_keywords, basis, aux_basis, run_index_start,number_of_runs, output_print_level)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank_println("running input file $input_file_path", rank)
  
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file_path; output=output_print_level)
    model_replace_dict = Dict{String,Any}("basis" => basis, "auxiliary_basis" => aux_basis)
    
    replace_model(model, model_replace_dict)
    replace_scf_keywords(keywords, scf_keywords)

    for run_id in run_index_start:run_index_start+number_of_runs-1
        name = split(input_file_name, ".")[1]
        name = "$(name)_run_$run_id"
        rank_println("Running $name", rank)

        mol, basis = JuliaChem.JCBasis.run(molecule, model; output=output_print_level)
        JuliaChem.JCMolecule.run(mol)

        run_time = @elapsed rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"]; output=output_print_level)
        timings = rhf_energy["Timings"]
       
        timings.run_name = name
        timings.run_time = run_time


        hdf5_file = HDF5.h5open(joinpath(output_dir, "$(name)-$(rank).h5"), "w")
        save_jc_timings_to_hdf5(timings, hdf5_file)
        close(hdf5_file)

        rhf_energy = nothing 


    end
end


function main(inputs)
    
    path_to_inputs = inputs[1]
    path_to_outputs = inputs[2]
    settings_id = inputs[3]
    basis_id = inputs[4]
    run_index_start = parse(Int, inputs[5])
    number_of_runs = parse(Int, inputs[6])



    input_files_paths, input_file_names = get_file_pathsV1(path_to_inputs)

    if length(inputs) >= 8
        file_start_index = parse(Int, inputs[7])
        file_end_index = parse(Int, inputs[8])
        input_files_paths = input_files_paths[file_start_index:file_end_index]
        input_file_names = input_file_names[file_start_index:file_end_index]
    elseif length(inputs) == 7
        file_start_index = parse(Int, inputs[7])
        input_files_paths = input_files_paths[file_start_index:end]
        input_file_names = input_file_names[file_start_index:end]
    end
    


    scf_keywords = get_settings(settings_id)

    if haskey(ENV, "df_exchange_n_blocks")
        scf_keywords["df_exchange_n_blocks"] = parse(Int, ENV["df_exchange_n_blocks"])
    end

    #options are divide_total_aux e.g. 1,2,4,8,9,10 , range_of_range_counts e.g. 1:10, list_of_range_counts e.g. 1,2,5,10 
    Q_range_mode_env = "divide_total_aux"
    if haskey(ENV, "mixed_df_Q_range_mode")
        Q_range_mode_env =  ENV["mixed_df_Q_range_mode"]
    end

    Q_ranges_string = "4"

    Q_range_iterable = 4:4 

    Q_divisors = Int[]
    if haskey(ENV, "mixed_df_Q_ranges")
        Q_ranges_string = ENV["mixed_df_Q_ranges"]
    end

    if Q_range_mode_env == "range_of_range_counts"
        Q_range_mode = "num_ranges"

        range_start, range_end = split(Q_ranges_string, ":")
        range_start = parse(Int, range_start)
        range_end = parse(Int, range_end)
        Q_range_iterable = range_start:range_end
        println("doing range of Q ranges $Q_range_iterable")
    elseif Q_range_mode_env == "list_of_range_counts"
        Q_range_mode = "num_ranges"
        Q_range_strings = split(Q_ranges_string, ",")
        Q_range_iterable = parse.(Int, Q_range_strings)
        println("doing list of Q ranges: $Q_range_iterable")
    else #assume Q_range_mode_env == "divide_total_aux"
        Q_range_mode = "divide_total_Q"
        Q_divisor_strings = split(Q_ranges_string, ",")
        Q_range_iterable = parse.(Int, Q_divisor_strings)
        println("doing divide total aux with divisors: $Q_range_iterable")
    end

   


    basis, aux_basis = get_basis_names(basis_id)

    JuliaChem.initialize()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    if rank == 0
        Pkg.status()
    end


    rank_println("starting run starting at index $(run_index_start)", rank)
    rank_println("Date of run: $(Dates.now())", rank)
    rank_println("settings_id: $settings_id", rank)
    rank_println("basis_id: $basis_id", rank)
    rank_println("basis: $basis", rank)
    rank_println("aux_basis: $aux_basis", rank)
    rank_println("BLAS config $(BLAS.get_config())", rank)

    if haskey(ENV, "df_adaptive_basis_limit")
        scf_keywords["df_adaptive_basis_limit"] = parse(Int, ENV["df_adaptive_basis_limit"])
        rank_println("df_adaptive_basis_limit: $(scf_keywords["df_adaptive_basis_limit"])", rank)
    end

    # scf_keywords["df_adaptive_basis_limit"] = 100

    print("scf_keywords: ")
    rank_display(scf_keywords, rank)


    rank_println("input_files_paths", rank)
    rank_println(input_files_paths, rank)
    rank_println("output file path", rank)
    rank_println(path_to_outputs, rank)


    
    rank_println("initialized JuliaChem", rank)
    flush(stdout)
    n_threads = Threads.nthreads()
    rank_println("Threads $n_threads", rank)

    if haskey(ENV, "OPENBLAS_NUM_THREADS") &&  ENV["OPENBLAS_NUM_THREADS"] != nothing
        rank_println("OPENBLAS_NUM_THREADS: $(ENV["OPENBLAS_NUM_THREADS"])", rank)
        open_blas_threads = parse(Int, ENV["OPENBLAS_NUM_THREADS"])
        BLAS.set_num_threads(open_blas_threads)
    end


    rank_println("BLAS Threads $(BLAS.get_num_threads())", rank)
    # BLAS.set_num_threads(n_threads)
    # println("done pinning threads")
    # create output directory if it does not exist
    if !isdir(path_to_outputs)
        mkdir(path_to_outputs)
    end

    scf_keywords["num_devices"] = 1

    output_print_level = 2
    for q_range_value in Q_range_iterable

        if Q_range_mode == "divide_total_Q"
            scf_keywords["Q_ranges_divide_Q_by"] = q_range_value
        else
           scf_keywords["num_Q_ranges"] = q_range_value
        end

        for i in eachindex(input_files_paths)
            #create folder for path_to_outputs/input_file_names[i]
            input_file_path = input_files_paths[i]
            input_file_name = input_file_names[i]
            output_dir = create_output_folderV2(path_to_outputs, input_file_path, rank)
            run_file(input_file_path, input_file_name, output_dir, scf_keywords, basis, aux_basis, run_index_start, number_of_runs, output_print_level)
        end
    end

    # JuliaChem.finalize()

end

main(ARGS)
