@info "Loading libraries and possibly one-time compiling";
using DataFrames, CSV;
include("functions.jl");

function main(args::Vector{String})

    if length(args) < 1
        @error("Please provide the name of the file containing the edge list of the network, e.g., networks/star.net");
        exit();
    end
    filein = args[1];

    # reads the file with the edge list
    N = readNetwork(filein);
    
    # convert to a matrix
    A = network2matrix(N);
    

    fileout = "results/"*basename(filein)*".csv";
    node_names = [N.idxmap[i] for i in 1:length(N.idxmap)];

    # estimates the fitness centrality
    if isdirected(N)
        Fout, Fin = asymmetricNHEFC(A);
        df = DataFrame(node=node_names, fitness_out = Fout, fitness_in = Fin);
    else
        F = symmetricNHEFC(A);
        df = DataFrame(node=node_names, fitness = F);
    end

    @info "Saving results into file \"$fileout\"";
    CSV.write(fileout, df);

    return;
end

main(ARGS)
