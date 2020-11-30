export ZDT, ZDT1, ZDT2, ZDT3, ZDT4, ZDT6;

# Problems taken from 
# Zitzler et al. - 2000 - Comparison of Multiobjective Evolutionary Algorithms
abstract type ZDT <: MOP end;

zdt_ref_string = "Zitzler et al. - 2000 - Comparison of Multiobjective Evolutionary Algorithms";

for i ∈ [1,2,3,4,6] 
    type_name = Symbol("ZDT", i)
    @eval begin 
        $("""
            struct ZDT$i <: ZDT 
                n_vars :: Int = 2;
            end

        Instance of a ZDT problem as given in 
        ( $zdt_ref_string )
        """)
        @with_kw struct $type_name <: ZDT
            n_vars :: Int = 2;
        end
        num_vars( mop :: $type_name ) = mop.n_vars;
        num_objectives( mop :: $type_name ) = 2;


        function get_pareto_set( mop :: $type_name )
            return SamplingFunction( mop, :ParetoSet, [:regular, :random] );
        end

        function get_pareto_front(  mop :: $type_name )
            return SamplingFunction( mop, :ParetoFront, [:regular, :random] );
        end 
    end
    if i ∈ [1,2,3,6]
        @eval constraints( mop :: $type_name ) = Box( zeros(mop.n_vars), ones( mop.n_vars) );
    elseif i == 4
        @eval begin 
            constraints( mop :: $type_name ) = Box( 
                [ 0; fill(-5.0, mop.n_vars - 1) ],
                [ 1.0; fill(5.0, mop.n_vars - 1) ]
            );
        end
    end
end

function get_objectives(mop::ZDT1)
    h(f1,g) = 1 - sqrt( f1 / g );
    g(x) = 1 + 9 / (num_vars(mop) - 1) * sum(x[2:end]);
    f1(x :: Vector{R} where R<:Real) = x[1]
    f2(x :: Vector{R} where R<:Real) = g(x) * h(f1(x), g(x))
    return [f1; f2]
end

function get_objectives(mop::ZDT2)
    h(f1,g) = 1 - ( f1 / g )^2;
    g(x) = 1 + 9 / (num_vars(mop) - 1) * sum(x[2:end]);
    f1(x :: Vector{R} where R<:Real) = x[1]
    f2(x :: Vector{R} where R<:Real) = g(x) * h(f1(x), g(x))
    return [f1; f2]
end

function get_objectives(mop::ZDT3)
    h(f1,g) = 1 - sqrt( f1 / g ) - (f1 / g) * sin(10 * pi * f1);
    g(x) = 1 + 9 / (num_vars(mop) - 1) * sum(x[2:end]);
    f1(x :: Vector{R} where R<:Real) = x[1]
    f2(x :: Vector{R} where R<:Real) = g(x) * h(f1(x), g(x))
    return [f1; f2]
end

function get_objectives(mop::ZDT4)
    h(f1,g) = 1 - sqrt( f1 / g )
    g(x) = 1 + 10 * ( num_vars(mop) -1 ) + sum( x[2:end].^2 -10*cos.( 4*pi *x[2:end]) )
    f1(x :: Vector{R} where R<:Real) = x[1]
    f2(x :: Vector{R} where R<:Real) = g(x) * h(f1(x), g(x))
    return [f1; f2]
end

function get_objectives(mop::ZDT6)
    h(f1,g) = 1 - ( f1 / g )^2;
    g(x) = 1 + 9 * ( sum(x[2:end])/(num_vars(mop)-1) )^0.25;
    f1(x :: Vector{R} where R<:Real) = x[1] - exp( -4*x[1] ) * (sin(6*pi*x[1]))^6;
    f2(x :: Vector{R} where R<:Real) = g(x) * h(f1(x), g(x));
    return [f1; f2]
end


function get_points( sample_func :: SamplingFunction, mop :: Union{ZDT1,ZDT2,ZDT4},
        ::Val{:ParetoSet}, num_points :: Int; method = :regular )
    
    if method == :regular
        x_range = num_points == 1 ? [0.0,] : range(0.0, 1.0; length = num_points );
    elseif method == :random
        x_range = [ rand() for i = 1 : num_points ]
    end
    return [ [x; zeros( num_vars(mop) -1 )] for x in x_range ]
end


function get_points( sample_func :: SamplingFunction, mop :: Union{ZDT6},
        ::Val{:ParetoSet}, num_points :: Int; method = :regular )
    # see 
    # https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt6/index.php
    x0 = .2807753191;
    if method == :regular
        x_range = num_points == 1 ? [x0,] : range(x0, 1.0; length = num_points );
    elseif method == :random
        x_range = [ x0 + (1-x0) * rand() for i = 1 : num_points ]
    end
    return [ [x; zeros( num_vars(mop) -1 )] for x in x_range ]
end

# ZDT 3 is a special snowflake …
function sample_from_range( arr, N, ::Val{:regular} )
    [ arr[1] + (arr[2] - arr[1]) * x for x in range(0, 1;length = N) ]
end

function sample_from_range( arr, N, ::Val{:random} )
    [ arr[1] + (arr[2] - arr[1]) * rand() for i = 1 : N ]
end

function get_points( sample_func :: SamplingFunction, mop :: ZDT3,
        ::Val{:ParetoSet}, num_points :: Int; method = :regular )
    # see 
    # https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt3/index.php

    chunk_size = num_points ÷ 5;
    
    x_ranges = [
        [0, 0.0830015349], 
            [[1e-10, 0.0] .+ v
            for
            v ∈ [
                [0.1822287280, 0.2577623634],
                [0.4093136748, 0.4538821041],
                [0.6183967944, 0.6525117038],
                [0.8233317983, 0.8518328654],
            ]]...
    ];

    x_range = Float64[];
    for i = 1 : 5
        if i == 1
            N = chunk_size + (num_points % chunk_size)
        else
            N = chunk_size
        end
        push!( x_range, sample_from_range( x_ranges[i], N, Val(method) )... )
    end
    return [ [x; zeros( num_vars(mop) -1 )] for x in x_range ]
end
