using JLD2, NCDatasets

file1 = jldopen("output/sim/thread1_floes.jld2", "r")
file2 = jldopen("output/sim/thread8_floes.jld2", "r")

function compare_floe_data(filename1, filename2)
    file1 = jldopen(filename1, "r")
    file2 = jldopen(filename2, "r")
    @assert length(keys(file1)) == length(keys(file2))
    for k in keys(file1)
        if k!= "metadata"
            println(string("Comparing key: ", k))
            for t in keys(file1[k])
                f1_data = file1[k][t]
                f2_data = file2[k][t]
                same = true
                for i in eachindex(f1_data)
                    same = same && all(f1_data[i] .== f2_data[i])
                    if !same
                        println(string(
                            "Differences starting at time ",
                            t,
                            " and index ",
                            i,
                        ))
                        break
                    end
                end
            end
        end
    end
    close(file1)
    close(file2)
    return
end


function compare_grid_data(filename1, filename2)
    file1 = Dataset(filename1, "r")
    file2 = Dataset(filename2, "r")
    @assert length(keys(file1)) == length(keys(file2))
    for k in keys(file1)
        if k!= "time" && k!= "x" && k!= "y"
            println(string("Comparing key: ", k))
            f1_data = file1[k]
            f2_data = file2[k]
            same = true
            same = same && all(f1_data[:, :, :] .== f2_data[:, :, :])
            if !same
                println(string("Differences in ", k))
                break
            end
        end
    end
    close(file1)
    close(file2)
end

function compare_checkpointer_data(filename1, filename2)
    file1 = jldopen(filename1, "r")
    file2 = jldopen(filename2, "r")

    ocean1 = file1["ocean"]
    ocean2 = file2["ocean"]
    fields =fieldnames(typeof(Ocean))
    println("Comparing Ocean")
    for t in keys(ocean1)
        for f in fields
            same = true
            same = same && all(getfield(ocean1[t], f) .== getfield(ocean2[t], f))
            if !same
                println(string("Differences in field ", f, "and timestep ", t))
                break
            end
        end
    end

    atmos1 = file1["atmos"]
    atmos2 = file2["atmos"]
    fields =fieldnames(typeof(Atmos))
    println("Comparing Atmosphere")
    for t in keys(ocean1)
        for f in fields
            same = true
            same = same && all(getfield(atmos1[t], f) .== getfield(atmos2[t], f))
            if !same
                println(string("Differences in field ", f, "and timestep ", t))
                break
            end
        end
    end

    close(file1)
    close(file2)
end


