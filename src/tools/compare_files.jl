using JLD2, NCDatasets

file1 = jldopen("output/sim/thread1_floes.jld2", "r")
file2 = jldopen("output/sim/thread8_floes.jld2", "r")

for k in keys(file1)
    if k!= "metadata"
        println(k)
        for t in keys(file1[k])
            f1_data = file1[k][t]
            f2_data = file2[k][t]
            same = true
            for i in eachindex(f1_data)
                same = same && all(f1_data[i] .== f2_data[i])
                if !same
                    println(string(t, " ", i))
                    break
                end
            end
        end
    end
end

close(file1)
close(file2)

file1 = Dataset("output/sim/thread1_grid.nc","r")
file2 = Dataset("output/sim/thread8_grid.nc", "r")

for k in keys(file1)
    if k!= "time" && k!= "x" && k!= "y"
        f1_data = file1[k]
        f2_data = file2[k]
        same = true
        same = same && all(f1_data[:, :, :] .== f2_data[:, :, :])
        if !same
            println(k)
            break
        end
    end
end

close(file1)
close(file2)

file1 = jldopen("output/sim/thread1_checkpointer.jld2", "r")
file2 = jldopen("output/sim/thread8_checkpointer.jld2", "r")

ocean1 = file1["ocean"]
ocean2 = file2["ocean"]
for t in keys(ocean1)
    same = true
    same = same && all(ocean1[t].fx .== ocean2[t].fx)
    same = same && all(ocean1[t].fy .== ocean2[t].fy)
    if !same
        println(t)
        break
    end
end

close(file1)
close(file2)


