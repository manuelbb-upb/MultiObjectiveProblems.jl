export MHT, MHT1, MHT2, MHT3,MHT4,MHT5,MHT6;
export CheapFunction, ExpensiveFunction;

# Wrappers for Heterogenous Problems
struct CheapFunction <: Function 
    f :: Function 
end

struct ExpensiveFunction <: Function 
    f :: Function 
end

(F::CheapFunction)( args...; kwargs... ) = F.f( args...; kwargs... );
(F::ExpensiveFunction)( args...; kwargs... ) = F.f( args...; kwargs...);

####### General Setup ###########
mht_ref_string = "Thomann, 2019 - Doctoral Thesis: A Trust Region Approach for Multi-Objective Heterogenous Optimization";
@doc """
    abstract type MHT <: MOP end

Abstract super type for the set of test problems MHT1-MHT7 
(called T1-T7 in the source) from 
$(mht_ref_string)
"""
abstract type MHT <: MOP end;

# some ipmlementations that some problems share 

function get_pareto_set( mop ::Union{MHT1, MHT3, MHT4,MHT5,MHT6} )
    return SamplingFunction( mop, :ParetoSet, [:regular, :random] );
end
function get_pareto_front( mop :: Union{MHT1, MHT3, MHT4,MHT5,MHT6} )
    return SamplingFunction( mop, :ParetoFront, [:regular, :random] );
end
########## MHT1 ##############
@doc """
    MHT1()

Bi-Objective, bi-variate problem MHT1 (T1) from
$(mht_ref_string).

Unconstrained, convex objectives, convex Pareto Front.
``f₁`` expensive, ``f₂`` cheap.

The objectives are given as 
```math 
\\begin{aligned}
    f_1( \\mathbf{x} ) &= \\frac{1}{2}x_1^2 + x_2^2 - 10 x_1 - 100,\\\\ 
    f_2( \\mathbf{x} ) &= x_1^2 + \\frac{1}{2}x_2^2 - 10 x_2 - 100.
\\end{aligned}
```
"""
struct MHT1 <: MHT end;

constraints( mop :: MHT1 ) = nothing;
num_vars( mop :: MHT1 ) = 2;
num_objectives( mop :: MHT1 ) = 2;

function get_objectives( mop :: MHT1 )
    f1 = ExpensiveFunction(
        function( x :: Vector{R} where R<:Real)
            0.5 * x[1]^2 + x[2]^2 - 10 * x[1] - 100
        end
    );
    f2 = CheapFunction( 
        function( x :: Vector{R} where R<:Real)
            x[1]^2 + 0.5 * x[2]^2 - 10 * x[2] - 100
        end
    );
    return [f1, f2]
end

get_ideal_point( mop :: MHT1 ) = fill(-150.0, 2);

mht1_points = [
    [-34.24; 47.06],
    [45.72; -1.46],
    [30.03; -35.81],
    [-7.82; 41.57],
    [29.22; 45.95],
    [15.75; -46.43],
    [34.91; 43.4],
    [17.87; 25.77],
    [24.31; -10.78],
    [15.55; -32.88],
]
function get_starting_point( mop :: MHT1, N :: Int64 )
    n = 1 + (N-1) % 10
    return mht1_points[n]
end

########## MHT2 ##############
@doc """
    MHT2()

Bi-Objective, bi-variate problem MHT2 (T2) from
$(mht_ref_string).

Unconstrained, non-convex objectives, non-convex Pareto Front.
``f₁`` cheap, ``f₂`` expensive.

The objectives are given as 
```math 
\\begin{aligned}
    f_1( \\mathbf{x} ) &= \\sin( x₂ ),\\\\ 
    f_2( \\mathbf{x} ) &= 1 - 
        \\exp \\left( 
            - \\left( x₁ - \\frac{1}{\\sqrt{2}} \\right)^2 
            - \\left( x₂ - \\frac{1}{\\sqrt{2}} \\right)^2 
        \\right)
\\end{aligned}
```
"""
struct MHT2 <: MHT end;

constraints( mop :: MHT2 ) = nothing;
num_vars( mop :: MHT2 ) = 2;
num_objectives( mop :: MHT2 ) = 2;

function get_objectives( mop :: MHT2 )
    f1 = CheapFunction(
        function( x :: Vector{R} where R<:Real)
            sin( x[2] )
        end
    );
    f2 = ExpensiveFunction( 
        function( x :: Vector{R} where R<:Real)
            1 - exp( -( x[1] - 1/sqrt(2) )^2 - ( x[2] - 1/sqrt(2) )^2 )
        end
    );
    return [f1, f2]
end

get_ideal_point( mop :: MHT2 ) = [-1.0; 0.0];

mht2_points = [
    [20.6; -46.82],
    [-22.31; -45.38],
    [-40.29; 32.35],
    [-.4; 1.0],
    [45.02; -46.56],
    [-6.13; -11.84],
    [-0.01; -0.51],
    [-31.31; -1.02],
    [0.51; 1.00],
    [.9; 1.51],
]
get_starting_point( ::MHT2, N :: Int64) = mht2_points[ 1 + (N-1) % 10 ];

########## MHT3 ##############
@doc """
    MHT3()

Bi-Objective, bi-variate problem MHT3 (T3) from
$(mht_ref_string).

Box constraints [-2,2]², convex objectives, convex Pareto Front.
``f₁`` cheap, ``f₂`` expensive.

The objectives are given as 
```math 
\\begin{aligned}
    f_1( \\mathbf{x} ) &= x₁ + 2,\\\\ 
    f_2( \\mathbf{x} ) &= x₁ - 2 + x₂.
\\end{aligned}
```
"""
struct MHT3 <: MHT end;

constraints( mop :: MHT3 ) = Box(
    fill(-2.0, 2),
    fill(2.0, 2)
);
num_vars( mop :: MHT3 ) = 2;
num_objectives( mop :: MHT3 ) = 2;

function get_objectives( mop :: MHT3 )
    f1 = CheapFunction(
        function( x :: Vector{R} where R<:Real)
            x[1] + 2
        end
    );
    f2 = ExpensiveFunction( 
        function( x :: Vector{R} where R<:Real)
            x[1] - 2 + x[2]
        end
    );
    return [f1, f2]
end

get_ideal_point( mop :: MHT3 ) = [ 0.0; -6.0];

mht3_points = [
    2.5767 -0.0793 -0.3193 0.0511 1.9058 0.8659 1.8695 -0.8956 2.2557 0.7349;
    1.6543 -0.3848 -1.1619 0.0646 1.7690 -0.7283 0.197 2.634 0.3009 0.5223
]
get_starting_point( ::MHT3, N :: Int64) = mht3_points[:,1 + (N-1) % 10 ];

########## MHT4 ##############
@doc """
    MHT4(n)

Bi-Objective, `n`-variate problem MHT4 (T4) from
$(mht_ref_string).

Box constraints [-10,10]ⁿ, convex objectives, convex Pareto Front.
``f₁`` expensive, ``f₂`` cheap.

The objectives are given as 
```math 
\\begin{aligned}
    f_1( \\mathbf{x} ) &= Σ_{i=1}^{n-1} xᵢ² + 2, \\\\ 
    f_2( \\mathbf{x} ) &= Σ_{i=1}^{n} xᵢ - 2
\\end{aligned}
```
"""
@with_kw struct MHT4 <: MHT
    n_vars :: Int64 = 2;
end;

constraints( mop :: MHT4 ) = Box(
    fill(-10.0, mop.n_vars),
    fill(2.0, mop.n_vars)
);
num_vars( mop :: MHT4 ) = mop.n_vars;
num_objectives( mop :: MHT4 ) = 2;

function get_objectives( mop :: MHT4 )
    f1 = ExpensiveFunction( 
        function( x :: Vector{R} where R<:Real )
            sum( x[1:end-1].^2 .+ 2.0 )
        end
    );
    f2 = CheapFunction(
        function( x :: Vector{R} where R<:Real)
            sum( x .- 2.0 )
        end
    );
    return [f1, f2]
end

get_ideal_point( mop :: MHT4 ) = [
    (mop.n_vars - 1)*2.0;
    -12.0 * mop.n_vars
];

mht4_points = [
    8.2401 -5.8961 0.0861 3.7051 -0.4188 2.4755 -9.0298 9.8070 2.5461 -8.1104;
    2.7963 9.4970 -4.333 -8.8077 0.7261 6.2651 3.0539 -0.4131 5.1145 -2.6839
]
function get_starting_point( mop::MHT4, N :: Int64) 
    if mop.n_vars == 2 
        return mht4_points[:,1 + (N-1) % 10 ];
    else
        return get_random_point( mop )
    end
end

########## MHT5 ##############
@doc """
    MHT5()

Bi-Objective, bi-variate problem MHT5 (T5) from
$(mht_ref_string).

Box constraints [1e-12,30] × [0,30],
convex objectives, convex Pareto Front.

``f₁`` expensive, ``f₂`` cheap.

The objectives are given as 
```math 
\\begin{aligned}
    f_1( \\mathbf{x} ) &= x₁ \\ln( x₁ ) + x₂²,\\\\ 
    f_2( \\mathbf{x} ) &= x₁² + x₂⁴.
\\end{aligned}
```
"""
struct MHT5 <: MHT end;

constraints( mop :: MHT5 ) = Box(
    [1e-12; 0.0],
    [30.0; 30.0]
);
num_vars( mop :: MHT5 ) = 2;
num_objectives( mop :: MHT5 ) = 2;

function get_objectives( mop :: MHT5 )
    f1 = ExpensiveFunction(
        function( x :: Vector{R} where R<:Real)
            x[1] * log( x[1] ) + x[2]^2
        end
    );
    f2 = CheapFunction( 
        function( x :: Vector{R} where R<:Real)
            x[1]^2 + x[2]^4
        end
    );
    return [f1, f2]
end

get_ideal_point( mop :: MHT5 ) = [ -1/ℯ; 0.0];

mht5_points = [
    25.9816 4.2130 7.7573 2.1431 10.4985 22.1153 18.7311 29.6945 25.6892 14.5115;
    2.6608 19.5013 24.4573 8.8787 6.5516 9.8276 9.1069 23.0979 25.2217 3.8922
]
get_starting_point( ::MHT5, N :: Int64) = mht5_points[:,1 + (N-1) % 10 ];

########## MHT6 ##############
@doc """
    MHT6()

Bi-Objective, bi-variate problem MHT6 (T6) from
$(mht_ref_string).

Box constraints [1e-12,100]²,
convex objectives, convex Pareto Front.

``f₁`` expensive, ``f₂`` cheap.

The objectives are given as 
```math 
\\begin{aligned}
    f_1( \\mathbf{x} ) &= -\\ln( x₁ ) - \\ln( x₂ ),\\\\ 
    f_2( \\mathbf{x} ) &= x₁² + x₂.
\\end{aligned}
```
"""
struct MHT6 <: MHT end;

constraints( mop :: MHT6 ) = Box(
   fill(1e-12,2),
   fill(100.0,2)
);

num_vars( mop :: MHT6 ) = 2;
num_objectives( mop :: MHT6 ) = 2;

function get_objectives( mop :: MHT6 )
    f1 = ExpensiveFunction(
        function( x :: Vector{R} where R<:Real)
            - log( x[1] ) - log( x[2] )
        end
    );
    f2 = CheapFunction( 
        function( x :: Vector{R} where R<:Real)
            x[1]^2 + x[2]
        end
    );
    return [f1, f2]
end

get_ideal_point( mop :: MHT6 ) = [-2*log(100);0.0];

mht6_points = [
   9.5896 95.6433 61.6860 23.3495 44.4905 5.8716 72.6812 25.441 42.2296 44.4367;
   9.4335 20.8035 50.6963 6.1688 2.7671 52.5858 58.64 7.8104 35.9188 17.2381
]
get_starting_point( ::MHT6, N :: Int64) = mht6_points[:,1 + (N-1) % 10 ];

##########################################################################
# Pareto Sets when individual minima / other points are connected by a line

function individual_minima(::MHT1)
    x_min1 = [ -10.0 ;  0 ];
    x_min2 = [ 0; -10.0 ];
    x_diff = x_min2 .- x_min1
    return x_min1, x_min2, x_diff
end

function individual_minima(::MHT3)
    x_min1 = [ -2.0 ;  2.0 ];
    x_min2 = [ -2.0 ; -2.0 ];
    x_diff = x_min2 .- x_min1
    return x_min1, x_min2, x_diff
end

function individual_minima(mop::MHT4)
    x_min1 = [ zeros(mop.n_vars-1); -10.0 ];
    x_min2 = fill( -10.0, mop.n_vars );
    x_diff = x_min2 .- x_min1
    return x_min1, x_min2, x_diff
end

function individual_minima(mop::MHT5)
    x_min1 = [ 1/ℯ, 0.0];
    x_min2 = [1e-12, 0.0];
    x_diff = x_min2 .- x_min1
    return x_min1, x_min2, x_diff
end

function individual_minima(mop::MHT6)
    x_min1 = [100.0; 100.0];
    x_min2 = [0.0; 100.0]
    x_diff = x_min2 .- x_min1
    return x_min1, x_min2, x_diff
end

function get_points( sample_func :: SamplingFunction, mop :: Union{MHT1,MHT3,MHT4,MHT5,MHT6}, 
    ::Val{:ParetoSet}, num_points :: Int; method = :regular )

    if method == :regular
        t_range = num_points == 1 ? [0.5,] : range(0, 1.0; length = num_points );
    elseif method == :random
        t_range = rand(num_points);
    end

    x_min1, x_min2, x_diff = individual_minima( mop );   

    return [ x_min1 .+ τ .* x_diff for τ in t_range ]
end

function get_points( sample_func :: SamplingFunction, mop :: Union{MHT1,MHT3,MHT4,MHT5,MHT6},
    ::Val{:ParetoFront}, num_points :: Int; method = :regular )

    pset = get_points( sample_func, mop, Val(:ParetoSet), num_points; method );
    F = get_vector_objective( mop );
    return F.( pset );
end
