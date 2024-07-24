@info "Loading libraries";
using DataFrames, CSV;
include("functions.jl");

function main(filein::String)
    # reads the file with the edge list
    N = readNetwork(filein);
    # convert to a matrix
    A = network2matrix(N);
    # estimates the fitness centrality
    steps, F = symmetricNHEFC(A);

    println("The algorithm converged in $steps steps");

    fileout = "results/"*basename(filein)*".csv";

    node_names = [N.idxmap[i] for i in 1:length(N.idxmap)];

    df = DataFrame(node=node_names, fitness = F);

    @info "Saving results into file \"$fileout\"";
    CSV.write(fileout, df);

    return F;
end

main("networks/star.net")