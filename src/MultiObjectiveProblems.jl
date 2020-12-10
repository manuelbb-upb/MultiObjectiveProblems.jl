module MultiObjectiveProblems

import JuMP
using OSQP
import ForwardDiff
using Parameters: @with_kw

import NLopt;    # for ZDT4 local critical set

export MOP, SamplingFunction, FixedPointSet, Constraints, Box
export get_objectives, get_vector_objective, get_gradients, get_omega_function, constraints,
    get_pareto_set, get_pareto_front, get_points, get_scatter_arrays, get_scatter_points,
    get_random_point, get_ideal_point, get_critical_set, get_critical_front;
export sampling_methods;

abstract type Constraints end;
struct Box <:Constraints 
    lb :: Vector{Float64}
    ub :: Vector{Float64}
end

@doc """
        abstract type MOP end;

    A supertype for all problem types.
    Each problem must implement the following methods:
    * `num_vars`
    * `num_objectives`
    * `get_objectives`
"""
abstract type MOP end;
Broadcast.broadcastable(m::M where M<:MOP) = Ref(m)

# required
@doc "Return on array of multivariate objective functions."
function get_objectives( :: MOP ) end;
@doc "Return the number of decision variables."
function num_vars( :: MOP) end;
@doc "Return the number of objective functions."
function num_objectives( :: MOP ) end;
    
# non-mandatory but generically implemented


@doc "Return a constraint object for `mop` or nothing."
function constraints( mop :: M where M <: MOP ) end;
function get_vector_objective( mop :: M where M<:MOP ) 
    func_array = get_objectives( mop );
    return (x::Vector{R} where R<:Real) -> [ f(x) for f in func_array ];
end

function get_gradients( mop :: M where M <: MOP )
    return [ (x :: Vector{R} where R<:Real ) -> ForwardDiff.gradient( f, x )
        for f in get_objectives(mop) ]
end

function get_omega_function( mop :: M where M <: MOP)
    omega_function = function( x :: Vector{R} where {R<:Real} )
        try
            ∇F = transpose( hcat( ( ∇f(x) for ∇f in get_gradients( mop ) )... ) );

            prob = JuMP.Model( OSQP.Optimizer );
            
            JuMP.set_silent(prob);
            JuMP.set_optimizer_attribute(prob,"eps_rel",1e-5);
            JuMP.set_optimizer_attribute(prob,"eps_abs", 1e-5);
            JuMP.set_optimizer_attribute(prob,"polish",true);

            JuMP.@variable(prob, α );     # negative of marginal problem value
            JuMP.@variable(prob, d[1:num_vars(mop)] );   # direction vector
        
            JuMP.@objective(prob, Min, α);

            JuMP.@constraint(prob, descent_contraints, ∇F*d .<= α);
            JuMP.@constraint(prob, norm_constraints, -1.0 .<= d .<= 1.0);
            
            global_constraints = constraints(mop);
            if isa( global_constraints, Box )
                JuMP.@constraint(prob, box_constraints,  global_constraints.lb .<= x .+ d .<= global_constraints.ub );
            end

            JuMP.optimize!(prob)
            ω = -JuMP.value(α);
            return ω 
        catch e
            return -Inf
        end
    end
    return omega_function
end

function get_random_point( mop:: M where M <: MOP )
    if isa( constraints(mop), Box )
        lb = constraints(mop).lb; 
        w = constraints(mop).ub .- lb;
        return lb .+ w .* rand( length(lb) );
    end
end

# non-mandatory and not implemented by default
@doc "Return the ideal (or an utopia) point of the problem."
function get_ideal_point( ::MOP ) end;

@doc """
    get_pareto_set( mop )

Get an object representing the Pareto Set of `mop`.
Possible return types are `FixedPointSet` and `SamplingFunction`.

See also: [`get_points`](@ref), [`get_scatter_points`](@ref)
"""
function get_pareto_set( :: MOP ) end;

@doc """
    get_pareto_front( mop )

Get an object representing the Pareto Front of `mop`.
Possible return types are `FixedPointSet` and `SamplingFunction`.

See also: [`get_points`](@ref), [`get_scatter_points`](@ref)
"""
function get_pareto_front( :: MOP ) end;

function get_critical_set( :: MOP ) end;
function get_critical_front( :: MOP ) end;

@doc "Abstract Super Type for Pareto Set and Pareto Frontier."
abstract type CompareSet end;
function get_points( ::CompareSet, ::Int; kwargs...) end;

@doc "A wrapper around a `list_of_points`."
struct FixedPointSet <: CompareSet
    list_of_points :: Vector{Vector{R}} where R<:Real;
end;     # for calculated point clouds

@doc """
    A `SamplingFunction` enables you to sample points from the Pareto Set
    or the Pareto Front.
    The basic usage for `samp_func :: SamplingFunction` is 
    ```
    num_points = 10;
    list_of_points = get_points( samp_func, num_points );
    ```
    
    Check `sampling_methods(samp_func)` to see whether multiple sampling techniques
    are available.
    You can pass a sampling method via `method` keyword argument:
    ```
    list_of_points = get_points( samp_func, num_points; method = :regular );
    ```

    See also: [`sampling_methods`](@ref), [`get_points`](@ref), [`get_scatter_points`](@ref),
"""
@with_kw struct SamplingFunction <: CompareSet
    mop :: M where M <: MOP;
    type :: Symbol;
    methods :: Vector{Symbol};

    @assert type ∈ [:ParetoSet, :ParetoFront, :ParetoCriticalSet, :ParetoCriticalFront]
end

@doc "Return a list of Symbols that can be pased to `get_points` for sampling from 
a `SamplingFunction`."
sampling_methods( samp_func :: SamplingFunction ) = samp_func.methods;

@doc """
    get_points( samp_func, N; kwargs... )

Sample `N` points using `samp_func::SamplingFunction`.
Additional `kwargs` are passed down to the problem-specfic implementations.
"""
function get_points( samp_func :: SamplingFunction, N :: Int; kwargs... )
    return get_points( samp_func, samp_func.mop, Val(samp_func.type), N; kwargs...);
end

@doc "
    get_points( point_set, N; kwargs... )

Get all points using `point_set::FixedPointSet`."
function get_points( point_set :: FixedPointSet, N :: Int; kwargs... )
    println("Retrieving *all* points from a `FixedPointSet`, all other arguments are ignored.")
    return point_set.list_of_points;
end

@doc "Like `get_points` but returns data suited for plotting.\n
Specify variable axes by kwarg `dims` which default to `[1,2]`.

See also: [`get_points`](@ref)"
function get_scatter_points( samp_func :: SamplingFunction, N :: Int; dims = [1,2], kwargs... )
    get_scatter_arrays( get_points( samp_func, N; kwargs... ); dims = dims )
end

# not called directly
function get_points( samp_func :: SamplingFunction, mop :: M where M<:MOP, 
        ::Val{:ParetoFront}, N; kwargs... )
    pset = get_points( samp_func, samp_func.mop, Val(:ParetoSet), N; kwargs...);
    F = get_vector_objective( mop );
    return F.(pset) 
end 

# not called directly
function get_points( samp_func :: SamplingFunction, mop :: M where M<:MOP, 
    ::Val{:ParetoCriticalFront}, N; kwargs... )
pset = get_points( samp_func, samp_func.mop, Val(:ParetoCriticalSet), N; kwargs...);
F = get_vector_objective( mop );
return F.(pset) 
end 

# also not called directly
@doc "Return x and y data suited for scatter plots when `arr` is a list of vectors."
function get_scatter_arrays( arr :: Vector{Vector{R}} where R <: Real;
        dims :: Vector{Int} = [1, 2] )
    X = [ point[dims[1]] for point in arr ];
    Y = [ point[dims[2]] for point in arr ];
    if length( dims ) == 2
        return X, Y 
    else
        Z = [ point[dims[3]] for point in arr ];
        return X,Y,Z
    end
end

include("two_parabolas.jl");
include("ZDT_problems.jl");
include("DTLZ_problems.jl");
include("lovison_problems.jl")
    
end