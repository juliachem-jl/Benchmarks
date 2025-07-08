file1 = "S22/DF_RHF_denseCPU_cc-pvdz-ri_S22_1_fock/output.log"
file2 = "S22/DF_RHF_denseCPU_mixed_cc-pvdz-ri_divide_total_aux_S22_1_fock_divide8/output.log"

function get_fock_size_and_line(filename::String)
    lines = readlines(filename)  # Read all lines into memory
    for line in enumerate(lines)
        # find line of Fock size 
        if occursin("Matrix{Float64}:", line[2])
            fock_line = line[1] + 1
            line_wo_spaces_array = split(line[2], " ")
            fock_size = parse(Int64, split(line_wo_spaces_array[1], "×")[1]) 
            return fock_size, fock_line
        end
    end
    return nothing, nothing # Return if keyword not found
end

function get_fock_matrix(filename::String)
    fock_size, fock_line = get_fock_size_and_line(filename)
    F = zeros(fock_size, fock_size)
    lines = readlines(filename)

    matrix_line_index = 1
    for line in fock_line:(fock_line + fock_size - 1)
        line_data = split(lines[line])
        F[matrix_line_index, :] .= parse.(Float64, line_data)
        matrix_line_index += 1
    end

    return F
end

function compute_fock_difference(file1::String, file2::String)
    println("File 1 path: $file1")
    println("File 2 path: $file2")
    fock1 = get_fock_matrix(file1)
    fock2 = get_fock_matrix(file2)

    fock_difference = abs.(fock1 - fock2)

    min_diff = 100000
    max_diff = -100000
    avg_diff = 0
    for val in fock_difference
        if val < min_diff
            min_diff = val
        end
        if val > max_diff
            max_diff = val
        end
        avg_diff += val
    end
    avg_diff = avg_diff / length(fock_difference)

    println("Minimum difference between Fock elements: $min_diff")
    println("Maximum difference between Fock elements: $max_diff")
    println("Average difference between Fock elements: $avg_diff")
    println("Differences between Fock matrix elements:")
    display(fock_difference) 
end

display(compute_fock_difference(file1, file2))