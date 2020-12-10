export ZDT, ZDT1, ZDT2, ZDT3, ZDT4, ZDT6;

zdt_ref_string = "Zitzler et al. - 2000 - Comparison of Multiobjective Evolutionary Algorithms";
@doc """
    abstract type ZDT <: MOP end

Abstract super type for the class of ZDT problems as defined in
$(zdt_ref_string)

These are biobjective problems with a variable number of decision variables
`n_vars >= 2`.
They share a common structure
```math 
\\begin{aligned}
    f_1(\\mathbf{x}) &= f( x_1 ),\\\\ 
    f_2(\\mathbf{x}) &= g( x_2, …, x_n ) · h( f_1( x_1 ), g( x_2, …, x_n ) )
\\end{aligned}
```
"""
abstract type ZDT <: MOP end;

# some ipmlementations that all problems share 
num_vars( mop :: Z where Z<:ZDT ) = mop.n_vars; # make sure below that every implementation has field n_vars
num_objectives( mop :: Z where Z<:ZDT ) = 2;

get_pareto_set( mop :: Z where Z<:ZDT ) = SamplingFunction( mop, :ParetoSet, [:regular, :random] );
get_pareto_front( mop :: Z where Z<:ZDT ) = SamplingFunction( mop, :ParetoFront, [:regular, :random] );

function get_objectives( mop :: Z where Z<:ZDT ) 
    f1, g, h = zdt_helper_functions( mop )
    f2(x :: Vector{R} where R<:Real) = g(x) * h(f1(x), g(x))
    return [f1; f2]
end

########## ZDT 1 ##############
@doc """
    ZDT1( n_vars :: Int )

Bi-Objective problem ZDT1 with convex Pareto Front as given in
$(zdt_ref_string).

The objectives are given as 
```math 
\\begin{aligned}
    f_1( x_1 ) &= x_1,\\\\ 
    f_2( x_1, …, x_n) &= g( x_2, …, x_n ) · h( f_1( x_1 ), g( x_2, …, x_n ) ), \\\\
    g( x_2, …, x_n ) &= 1 + (9/(n-1)) · Σ_i x_i, \\\\
    h( f1, g ) &= 1 - \\sqrt{f_1/g}.
\\end{aligned}
```
"""
@with_kw struct ZDT1 <: ZDT 
    n_vars :: Int = 2;
end

constraints( mop :: ZDT1 ) = Box( zeros(mop.n_vars), ones( mop.n_vars) );

function zdt_helper_functions( mop :: ZDT1 )
    h(f1,g) = 1 - sqrt( f1 / g );
    g(x) = 1 + 9 / (num_vars(mop) - 1) * sum(x[2:end]);
    f1(x :: Vector{R} where R<:Real) = x[1]
    return [f1, g, h]
end

function get_ideal_point(mop::ZDT1)
    # min of f1 = x is mop.lb[1] = 0
    # min of g is 1, min of h is 1 - sqrt(ub[1]) = 0
    return zeros(2);
end

########## ZDT 2 ##############
@doc """
    ZDT2( n_vars :: Int )

Bi-Objective problem ZDT2 with **non**convex Pareto Front as given in
$(zdt_ref_string).

The objectives are given as 
```math 
\\begin{aligned}
    f_1( x_1 ) &= x_1,\\\\ 
    f_2( x_1, …, x_n) &= g( x_2, …, x_n ) · h( f_1( x_1 ), g( x_2, …, x_n ) ), \\\\
    g( x_2, …, x_n ) &= 1 + (9/(n-1)) · Σ_i x_i, \\\\
    h( f1, g ) &= 1 - (f_1/g)².
\\end{aligned}
```
"""
@with_kw struct ZDT2 <: ZDT 
    n_vars :: Int = 2;
end

constraints( mop :: ZDT2 ) = Box( zeros(mop.n_vars), ones( mop.n_vars) );

function zdt_helper_functions(mop::ZDT2)
    h(f1,g) = 1 - ( f1 / g )^2;
    g(x) = 1 + 9 / (num_vars(mop) - 1) * sum(x[2:end]);
    f1(x :: Vector{R} where R<:Real) = x[1]
    return [f1, g, h]
end

get_ideal_point(mop::ZDT2) = zeros(2);

########## ZDT 3 ##############
@doc """
    ZDT3( n_vars :: Int )

Bi-Objective problem ZDT3 with discrete Pareto Front as given in
$(zdt_ref_string).

The objectives are given as 
```math 
\\begin{aligned}
    f_1( x_1 ) &= x_1,\\\\ 
    f_2( x_1, …, x_n) &= g( x_2, …, x_n ) · h( f_1( x_1 ), g( x_2, …, x_n ) ), \\\\
    g( x_2, …, x_n ) &= 1 + (9/(n-1)) · Σ_i x_i, \\\\
    h( f1, g ) &= 1 - \\sqrt{f_1/g}- (f_1/g)·\\sin( 10 π f_1 )
\\end{aligned}
```
"""
@with_kw struct ZDT3 <: ZDT 
    n_vars :: Int = 2;
end

constraints( mop :: ZDT3 ) = Box( zeros(mop.n_vars), ones( mop.n_vars) );

function zdt_helper_functions(mop::ZDT3)
    h(f1,g) = 1 - sqrt( f1 / g ) - (f1 / g) * sin(10 * pi * f1);
    g(x) = 1 + 9 / (num_vars(mop) - 1) * sum(x[2:end]);
    f1(x :: Vector{R} where R<:Real) = x[1]
    
    return [f1, g, h]
end
get_ideal_point(mop::ZDT3) = zeros(2);

########## ZDT 4 ##############
@doc """
    ZDT4( n_vars :: Int )

Bi-Objective problem ZDT4 with many local Pareto Fronts as given in
$(zdt_ref_string).

The multimodality is due to ``g`` being a version of Rastrigin's function.

The objectives are given as 
```math 
\\begin{aligned}
    f_1( x_1 ) &= x_1,\\\\ 
    f_2( x_1, …, x_n) &= g( x_2, …, x_n ) · h( f_1( x_1 ), g( x_2, …, x_n ) ), \\\\
    g( x_2, …, x_n ) &= 1 + (9/(n-1)) · Σ_i ( x_i^2 - \\cos( 4 π x_i ) ), \\\\
    h( f1, g ) &= 1 - \\sqrt{f_1/g}
\\end{aligned}
```
"""
@with_kw struct ZDT4 <: ZDT 
    n_vars :: Int = 2;
end

constraints( mop :: ZDT4 ) = Box( 
    [ 0; fill(-5.0, mop.n_vars - 1) ],
    [ 1.0; fill(5.0, mop.n_vars - 1) ]
);

function zdt_helper_functions(mop::ZDT4)
    h(f1,g) = 1 - sqrt( f1 / g )
    g(x) = 1 + 10 * ( num_vars(mop) -1 ) + sum( x[2:end].^2 -10*cos.( 4*pi .* x[2:end]) )
    f1(x :: Vector{R} where R<:Real) = x[1]

    return [f1, g, h]
end
get_ideal_point(mop::ZDT4) = zeros(2);

get_critical_set( mop :: Z where Z<:ZDT ) = SamplingFunction( mop, :ParetoCriticalSet, [:regular, :random] );


########## ZDT 5 ##############
# not implemented yet due to it having binary variables 

########## ZDT 6 ##############
@doc """
    ZDT6( n_vars :: Int )

Bi-Objective problem ZDT6 with biased Pareto Front as given in
$(zdt_ref_string).

The objectives are given as 
```math 
\\begin{aligned}
    f_1( x_1 ) &= 1 - \\exp(-4 x_1 ) \\sin^6( 6 π x_1) ,\\\\ 
    f_2( x_1, …, x_n) &= g( x_2, …, x_n ) · h( f_1( x_1 ), g( x_2, …, x_n ) ), \\\\
    g( x_2, …, x_n ) &= 1 + 9 · Σ_i ( x_i / (n-1) )^{1/4}, \\\\
    h( f1, g ) &= 1 - (f_1/g)²
\\end{aligned}
```
"""
@with_kw struct ZDT6 <: ZDT 
    n_vars :: Int = 2;
end

constraints( mop :: ZDT6 ) = Box( zeros(mop.n_vars), ones( mop.n_vars) );

function zdt_helper_functions(mop::ZDT6)
    h(f1,g) = 1 - ( f1 / g )^2;
    g(x) = 1 + 9 * ( sum(x[2:end])/(num_vars(mop)-1) )^0.25;
    f1(x :: Vector{R} where R<:Real) = x[1] - exp( -4*x[1] ) * (sin(6*pi*x[1]))^6;
    
    return [f1, g, h]
end
get_ideal_point( mop::ZDT6 ) = zeros(2);

####### Point Sampling Functions #########################

#helper
function zdt_x_range(mop :: Z where Z<:ZDT, method :: Symbol, num_points::Int, 
    x1 :: Real, x2 :: Union{Nothing,Vector{R} where R<:Real} = nothing )
    
    if isnothing(x2)
        x2 = zeros( num_vars(mop) -1 );
    end

    @assert length(x2) == num_vars(mop) - 1;

    if method == :regular
        x_range = num_points == 1 ? [x1,] : range(x1, 1.0; length = num_points );
    elseif method == :random
        x_range = [ x1 + (1-x1) * rand() for i = 1 : num_points ] 
    end
    return [ [x; x2] for x in x_range ]
end

@doc """
    get_points( sample_func, zdt, Val(:ParetoSet), N; method = :regular )

Return a list of points ``x ∈ ℝⁿ`` where ``x₁ = … = xₙ = 0``.
Methods are `:regular` or `:random`.
"""
function get_points( sample_func :: SamplingFunction, mop :: Union{ZDT1,ZDT2,ZDT4},
        ::Val{:ParetoSet}, num_points :: Int; method = :regular )
    zdt_x_range(mop, method, num_points, 0.0)   
end

function get_points( sample_func :: SamplingFunction, mop :: Union{ZDT6},
        ::Val{:ParetoSet}, num_points :: Int; method = :regular )
    # see 
    # https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt6/index.php
    x1 = .2807753191;
    zdt_x_range( mop, method, num_points, x1)   
end

# ZDT 3 is a special snowflake …
function sample_from_range( arr, N, ::Val{:regular} )
    [ arr[1] + (arr[2] - arr[1]) * x for x in range(0, 1;length = N) ]
end

function sample_from_range( arr, N, ::Val{:random} )
    [ arr[1] + (arr[2] - arr[1]) * rand() for i = 1 : N ]
end

function get_points( :: SamplingFunction, mop :: ZDT3, ::Val{:ParetoSet}, 
        num_points :: Int; method = :regular )
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

function get_points( :: SamplingFunction, mop :: ZDT4, ::Val{:ParetoCriticalSet}, 
        num_points :: Int; x0 = nothing, method = :regular )
    #@eval using NLopt;

    if isnothing(x0) 
        x̂ = zeros(mop.n_vars - 1) # global optimum;
    else
        if length(x0) == mop.n_vars 
            x0 = x0[2:end];
        end
        
        if mop.n_vars == 2 && isa( x0, Real )
            x0 = [x0,]
        end

        _, g, _ = zdt_helper_functions( mop );
        # gradient via autodiff
        # g does not depend on x1 = 0.0
        ∇g = x -> ForwardDiff.gradient( g, [0.0; x] )[2:end];

        opt_func = function( x::Vector, grad::Vector)
            if length(grad) > 0
                grad[:] = ∇g(x)
            end
            return g(x)
        end

        box = constraints( mop );

        opt = NLopt.Opt( :LD_MMA, length(x0) );
        opt.lower_bounds = box.lb[2:end];
        opt.upper_bounds = box.ub[2:end];

        opt.xtol_rel = 1e-5;
        opt.maxeval = 200 * length(x0);

        println("Starting search for local minimum of g for ZDT4 with x0:")
        
        opt.min_objective = opt_func;
        _, x̂, ret =  NLopt.optimize( opt, x0 );
        #@show x̂
        println("Finished with $(opt.numevals) eval(s).")
    end
    
    return zdt_x_range( mop, method, num_points, 0.0, x̂ );
end

get_critical_front( mop :: Z where Z<:ZDT ) = SamplingFunction( mop, :ParetoCriticalFront, [:regular, :random] );
