export TwoParabolasConstrained, TwoParabolasUnconstrained;

struct TwoParabolasUnconstrained <: MOP end;
@with_kw struct TwoParabolasConstrained <: MOP
    lb :: Union{R,Vector{R}} where R<:Real = -4.0;
    ub :: Union{R,Vector{R}} where R<:Real = 4.0;
 end;

num_vars( ::Union{TwoParabolasUnconstrained, TwoParabolasConstrained} ) = 2;
num_objectives( ::Union{TwoParabolasUnconstrained, TwoParabolasConstrained} ) = 2;

function get_objectives( ::Union{TwoParabolasUnconstrained, TwoParabolasConstrained} )
    f1( x :: Vector{R} where R<:Real ) = sum( ( x .- 1.0 ).^2 );
    f2( x :: Vector{R} where R<:Real ) = sum( ( x .+ 1.0 ).^2 );
    return [f1, f2]
end

function constraints( mop::TwoParabolasConstrained )
    Box( 
        isa( mop.lb, Real) ? fill( mop.lb, num_vars(mop) ) : mop.lb,
        isa( mop.ub, Real) ? fill( mop.ub, num_vars(mop) ) : mop.ub,
    )
end

function get_pareto_set( mop ::Union{TwoParabolasUnconstrained, TwoParabolasConstrained} )
    return SamplingFunction( mop, :ParetoSet, [:regular, :random] );
end

function get_pareto_front( mop ::Union{TwoParabolasUnconstrained, TwoParabolasConstrained} )
    return SamplingFunction( mop, :ParetoFront, [:regular, :random] );
end 

function get_points( sample_func :: SamplingFunction, :: Union{TwoParabolasUnconstrained, TwoParabolasConstrained},
        ::Val{:ParetoSet}, num_points :: Int; method = :regular )
    
    if method == :regular
        x_range = num_points == 1 ? [0.0,] : range(-1.0, 1.0; length = num_points );
    elseif method == :random
        x_range = [ -1.0 + 2 * rand() for i = 1 : num_points ]
    end
    return [ [x;x] for x in x_range ]
end

function get_points( sample_func :: SamplingFunction, mop :: Union{TwoParabolasUnconstrained, TwoParabolasConstrained},
        ::Val{:ParetoFront}, num_points :: Int; method = :regular )

    pset = get_points( sample_func, mop, Val(:ParetoSet), num_points; method );
    F = get_vector_objective( mop );
    return F.( pset );
end