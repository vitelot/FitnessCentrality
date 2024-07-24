struct Network{T}
    directed::Bool;
    weighted::Bool;
    data::DataFrame;
    idxmap::Dict{Int, T};
end

isdirected(N::Network) = N.directed;
isweighted(N::Network) = N.weighted;



"""
Loads the network from an edge list and delivers a DataFrame.
Asks if the network is directed.
Determines if the network is weighted. 
The comment symbol in the file is % 
"""
function readNetwork(filein::String)
    @info "Reading the network";

    df = CSV.read(filein, comment="%", header=false, DataFrame);
    ncolumns = ncol(df);
    println("Found $ncolumns columns in network data file \"$filein\"");
    column_types = eltype.(eachcol(df));

    print("with data types ");
    for t in column_types
        print("$t ");
    end
    println();

    if ncolumns < 2
        println("Only one column found. We need at least two columns with node ids to create links.");
        println("Exiting.");
        exit();
    end

    println("Assuming the first two columns are node ids");

    print("Is the network directed [y/N]? ");

    DIRECTED = false;
    WEIGHTED = false;

    answer = readline();

    if answer == "y"
        println("Considering the network as directed");
        DIRECTED = true;
    else
        println("Considering the network as undirected");
    end

    if ncolumns == 2 
        println("Only two columns found: condidering the network as unweighted");
        rename!(df, [:from, :to]);
    elseif ncolumns == 3 
        println("Three columns found: condidering the network as weighted");
        rename!(df, [:from, :to, :weight]);
        WEIGHTED = true;
    else
        println("More than three columns found: condidering the network as weighted and ignoring all but the first three columns");
        select!(df, 1:3);
        rename!(df, [:from, :to, :weight]);
        WEIGHTED = true;
    end

    return Network(DIRECTED, WEIGHTED, df, Dict{Int, column_types[1]}() );
end

"""
Converts the network into a matrix.
"""
function network2matrix(N::Network)
    @info "Building the adjacency matrix";

    df = N.data;
    nodeList = unique(sort([df.from; df.to]));

    Idx = Dict{eltype(nodeList), Int}();
    inverseIdx = Dict{Int, eltype(nodeList)}();

    # assign a unique integer id to nodes
    idx = 1;
    for i in nodeList
        if !haskey(Idx, i)
            Idx[i] = idx;
            idx += 1;
        end
    end

    # inverseIdx = Dict(value => key for (key, value) in Idx);
    for (key,value) in Idx
        N.idxmap[value] = key;
    end

    A = zeros(length(nodeList), length(nodeList));

    if isweighted(N)
        if isdirected(N)
            for r in eachrow(df)
                from = Idx[r.from]; to = Idx[r.to]; w = r.weight;
                A[from, to] = w;
            end
        else
            for r in eachrow(df)
                from = Idx[r.from]; to = Idx[r.to]; w = r.weight;
                A[from, to] = A[to,from] = w;
            end
        end
    else
        if isdirected(N)
            for r in eachrow(df)
                from = Idx[r.from]; to = Idx[r.to];
                A[from, to] = 1;
            end
        else
            for r in eachrow(df)
                from = Idx[r.from]; to = Idx[r.to];
                A[from, to] = A[to,from] = 1;
            end
        end

    end

    return A;
end


"""
Non-homogeneous fitness and complexity algorithm.\n
Takes the matrix A, the inhomogeneous term δ and the required final relative error ϵ
and delivers the number of iterations needed to converge, the fitness centrality.
"""
function symmetricNHEFC(A::Matrix{T}; δ::Float64=1e-2, ϵ::Float64=1e-2) where T <:Real
    @info "Calculating fitness centrality";
    
    c,p = size(A);
    err::Float64 = 1000.0;
    nr_of_iterations::Int = 0;


    F0 = ones(T,c);
    F1 = ones(T,c);
    while err>ϵ
        nr_of_iterations += 1;
        F1 = δ .+ A  * inv.(F0);
        err = maximum(abs.( (F1 ./ F0) .- 1)); #maximum(abs.(vcat(F1-F0, S1-S0)));
        # @info err
        F0 = copy(F1);
    end

    nr_of_iterations, F1
end

