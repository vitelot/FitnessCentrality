struct Network{T}
    directed::Bool;
    weighted::Bool;
    data::DataFrame;
    nodenames::Vector{T};
end

isdirected(N::Network) = N.directed;
isweighted(N::Network) = N.weighted;

"non-homogeneous term in the map"
const DELTA   = 1e-2; 

"relative error. Map converges if the maximum relative error between two iterations is less than EPSILON."
const EPSILON = 1e-2; # 

"""
Loads the network from an edge list and delivers a DataFrame.
Asks if the network is directed.
Determines if the network is weighted. 
The comment symbol in the file is % 
"""
function readNetwork(filein::String, type::String="ask")
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
        @error("Only one column found. We need at least two columns with node ids to create links.\n\tExiting.");
        exit();
    end

    println("Assuming the first two columns are node ids");

    DIRECTED = false;
    WEIGHTED = false;

    if type == "directed"
        DIRECTED = true;
    elseif type == "undirected"
        DIRECTED = false;
    elseif type == "ask"
        print(stderr, "Is the network directed [y/N]? ");

        answer = readline();

        if answer == "y" || answer == "Y"
            println(stderr,"Considering the network as directed");
            DIRECTED = true;
        else
            println(stderr,"Considering the network as undirected");
        end
    else
        @error("Unknown $type network type");
        exit(1);
    end

    if ncolumns == 2 
        println("Only two columns found: condidering the network as unweighted");
        rename!(df, [:from, :to]);
    elseif ncolumns == 3 
        println("Three columns found: condidering the network as weighted");
        rename!(df, [:from, :to, :weight]);
        WEIGHTED = true;
    else
        @warn("More than three columns found: considering the network as weighted and ignoring all but the first three columns");
        select!(df, 1:3);
        rename!(df, [:from, :to, :weight]);
        WEIGHTED = true;
    end

    return Network(DIRECTED, WEIGHTED, df, Vector{column_types[1]}() );
end

"""
Converts the network into a matrix.
"""
function network2matrix(N::Network)
    @info "Building the adjacency matrix";

    df = N.data;
    nodeList = unique(sort([df.from; df.to]));

    Idx = Dict{eltype(nodeList), Int}();
    nodenames = Vector{eltype(nodeList)}();

    # assign a unique integer id to nodes
    idx = 1;
    for i in nodeList
        if !haskey(Idx, i)
            Idx[i] = idx;
            push!(nodenames, i);
            idx += 1;
        end
    end
      
    append!(N.nodenames, nodenames);

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

function U(A::Matrix{T}, V::Vector{T}, δ::Float64) where T <: Real
    invV = inv.(V);
    interaction = 0.5 * transpose(invV) * A * invV;
    logterm = sum(log.(V));
    deltaterm = δ * sum(invV);

    return interaction + logterm + deltaterm;
end

"""
Non-homogeneous fitness and complexity algorithm.\n
Takes the symmetric matrix A, the inhomogeneous term δ and the required final relative error ϵ
and delivers the fitness centrality.
"""
function symmetricNHEFC(A::Matrix{T}; 
            δ::Float64=DELTA, ϵ::Float64=EPSILON, 
            printpotential::Bool=false, randominit::Bool=false,
            outputonfile::Bool=false,
        )::Vector{T} where T <:Real
    
    @info "Using the plain symmetric algorithm";

    if outputonfile
        OUT = open("output.csv","a");
    else
        OUT = stdout;
    end

    c,p = size(A);
    err::Float64 = 1000.0;
    nr_of_iterations::Int = 0;

    randomID = rand('A':'Z', 6) |> String;
    if randominit
        F0 = rand(c) ./ δ;
    else
        F0 = ones(c);
    end

    while err>ϵ
        nr_of_iterations += 1;
      
        F1 = δ .+ A  * inv.(F0);
        err = maximum(abs.( (F1 ./ F0) .- 1)); #maximum(abs.(vcat(F1-F0, S1-S0)));
        # @info err
        F0 = copy(F1);
        
        printpotential && println(OUT, "#$randomID# $nr_of_iterations $(U(A,F1,δ)) $err" );

        if nr_of_iterations%1000 == 0
            println("Iteration: $nr_of_iterations -- RelativeError: $err");
        end
    end
    
    outputonfile && close(OUT);

    println("The algorithm converged in $nr_of_iterations steps");
    
    F0
end

"""
Non-homogeneous fitness and complexity algorithm.\n
Takes the symmetric matrix A, the inhomogeneous term δ and the required final relative error ϵ
and delivers the fitness centrality.\n
To avoid oscillations it moves two steps a time.
"""
function twostepSymmetricNHEFC(A::Matrix{T}; 
        δ::Float64=DELTA, ϵ::Float64=EPSILON, 
        printpotential::Bool=false, randominit::Bool=false
    )::Vector{T} where T <:Real

    @info "Using the two step algorithm";
    c,p = size(A);
    err::Float64 = 1000.0;
    nr_of_iterations::Int = 0;

    randomID = rand('A':'Z', 6) |> String;
    if randominit
        F0 = rand(c) ./ δ;
    else
        F0 = ones(c);
    end

    F0 = δ .+ A  * inv.(F0);
    while err>ϵ
        nr_of_iterations += 1;
        F00 = δ .+ A  * inv.(F0);
        F1 = δ .+ A  * inv.(F00);
        err = maximum(abs.( (F1 ./ F0) .- 1)); #maximum(abs.(vcat(F1-F0, S1-S0)));
        # @info err
        F0 = copy(F1);
        
        printpotential && println("#$randomID# $nr_of_iterations $(U(A,F1,δ)) $err" );

        if nr_of_iterations%1000 == 0
            println("Iteration: $nr_of_iterations -- RelativeError: $err");
        end
    end

    println("The algorithm converged in $nr_of_iterations steps");
    
    F0
end

"""
Non-homogeneous fitness and complexity algorithm.\n
Takes the symmetric matrix A, the inhomogeneous term δ and the required final relative error ϵ
and delivers the fitness centrality.\n
It uses the gradient of the potential energy using the variable Z = -log(V)
"""
function symmetricGradientNHEFC(A::Matrix{T}; 
            δ::Float64=DELTA, ϵ::Float64=EPSILON, 
            printpotential::Bool=false, randominit::Bool=false
    )::Vector{T} where T <:Real
    
    
    @info "Using the gradient map in Z";

    c,p = size(A);
    err::Float64 = 1000.0;
    nr_of_iterations::Int = 0;

    randomID = rand('A':'Z', 6) |> String;
    if randominit
        Z0 = -log.(rand(c) ./ δ);
    else
        Z0 = zeros(c);
    end

    while err>ϵ
        nr_of_iterations += 1;
        EZ0 = exp.(Z0);
        Z1 = 1 .+ Z0 - δ .* EZ0 - EZ0 .* (A  * EZ0);
        err = maximum(abs.( (Z1 ./ Z0) .- 1));
        # @info err
        Z0 = copy(Z1);
        F1 = exp.(-Z1);
        printpotential && println("#$randomID# $nr_of_iterations $(U(A,F1,δ)) $err" );

        if nr_of_iterations%1000 == 0
            println("Iteration: $nr_of_iterations -- RelativeError: $err");
        end
    end

    println("The algorithm converged in $nr_of_iterations steps");
    
    exp.(-Z0)
end

"""
Non-homogeneous fitness and complexity algorithm.\n
Takes the symmetric matrix A, the inhomogeneous term δ and the required final relative error ϵ
and delivers the fitness centrality.\n
It uses the gradient of the potential energy using the variable V
"""
function symmetricDirectGradientNHEFC(A::Matrix{T}; 
            δ::Float64=DELTA, ϵ::Float64=EPSILON, 
            printpotential::Bool=false, randominit::Bool=false
    )::Vector{T} where T <:Real

    @info "Using the gradient map in V";

    c,p = size(A);
    err::Float64 = 1000.0;
    nr_of_iterations::Int = 0;

    randomID = rand('A':'Z', 6) |> String;
    if randominit
        V0 = rand(c) ./ δ;
    else
        V0 = ones(c);
    end

    while err>ϵ
        nr_of_iterations += 1;
        V0_inv = inv.(V0);
        V1 = V0 .* exp.(-1 .+ (V0_inv .* (δ .+ A*V0_inv)));
        
        # err = abs( (U(A,V0,δ) / U(A,V1,δ)) - 1);
        err = maximum(abs.( (V1 ./ V0) .- 1));
        # @info err
        V0 = copy(V1);
        printpotential && println("#$randomID# $nr_of_iterations $(U(A,V0,δ)) $err" );

        if nr_of_iterations%1000 == 0
            println("Iteration: $nr_of_iterations -- RelativeError: $err");
        end
    end

    println("The algorithm converged in $nr_of_iterations steps");
    
    V0
end

"""
Non-homogeneous fitness and complexity algorithm.\n
Takes the symmetric matrix A, the inhomogeneous term δ and the required final relative error ϵ
and delivers the fitness centrality.
"""
function asymmetricNHEFC(A::Matrix{T}; δ::Float64=1e-2, ϵ::Float64=1e-2) where T <:Real
    @info "Calculating fitness centrality in the directed case";

    c,p = size(A);
    err::Float64 = 1000.0;
    nr_of_iterations::Int = 0;

    F0 = ones(T,c);
    F1 = ones(T,c);
    S0 = ones(T,c);
    S1 = ones(T,c);

    AT = A';

    while err>ϵ
        nr_of_iterations += 1;
        F1 = δ .+ A  * inv.(S0);
        S1 = δ .+ AT * inv.(F0);
        err1 = maximum(abs.( (F1 ./ F0) .- 1)); #maximum(abs.(vcat(F1-F0, S1-S0)));
        err2 = maximum(abs.( (S1 ./ S0) .- 1)); #maximum(abs.(vcat(F1-F0, S1-S0)));
        err = max(err1, err2);
        # @info err

        F0 = copy(F1);
        S0 = copy(S1);
        if nr_of_iterations%1000 == 0
            println("Iteration: $nr_of_iterations -- RelativeError: $err");
        end
    end

    println("The algorithm converged in $nr_of_iterations steps");
    
    F1,S1
end
