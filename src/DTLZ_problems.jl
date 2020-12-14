# Problems taken from 
# Deb et. al - "Scalable Test Problems for Evolutionary Multiobjective Optimization"

# DTLZ8 & DTLZ9 not supported yet due to constraints...
export DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7, DTLZ;
export get_ideal_point;

dtlz_ref_string = """Deb et. al - "Scalable Test Problems for Evolutionary Multiobjective Optimization"."""

#%% -------------- Problem definitions ----------------------------------#
abstract type DTLZ <: MOP end;

dtlz_descriptor( i, n_vars, M) = """
    struct DTLZ$i <: DTLZ
        n_vars :: Int = $n_vars;
        M :: Int = $M ;
        k :: Int = n_vars + 1 - M;
    end

Test problem defined in: $dtlz_ref_string\n
Has `n_vars` decision variables and `M` objectives.
The decision vector `x` is split into two parts: `x[1:M-1]` and `x[M:end]`.
The latter has length `k`.
"""

pareto_front_descriptions = Dict( 
    1 => "DTLZ1 has an affine-linear global Pareto-optimal front with Σfᵢ=0.5 and 
    **many** local minima because a function `g` akin to Rastrigin's is used.",
    2 => "DTLZ2 has spherical global Pareto-optimal front with Σfᵢ²=1.",
    3 => "DTLZ3 has spherical global Pareto-optimal front with Σfᵢ²=1 and 
    **many** local minima because a function `g` akin to Rastrigin's is used.",
    4 => "DTLZ4 has spherical global Pareto-optimal front with Σfᵢ²=1 and 
    solutions tend to lie near the ``f_M-f₁`` plane.",
    5 => "DTLZ5 has spherical Pareto-optimal front that is a curve.",
    6 => "DTLZ6 has spherical Pareto-optimal front that is a curve and hard to reach.",
    7 => "DTLZ7 has ``2^{M-1}`` disconnected Pareto-optimal regions in search space.
    It is similar in construction to the ZDT problems.",
);

# dtlz index => (n_vars, M)
dtlz_standard_settings = Dict( 
    1 => (10,6),
    2 => (15,6),
    3 => (15,6),
    4 => (15,6),
    5 => (15,6),
    6 => (15,6),
    7 => (22,3)
)

for i = 1 : 7    
       
    n_vars, M = dtlz_standard_settings[i];
    help_text = """
    $(dtlz_descriptor(1, n_vars, M))

    $(pareto_front_descriptions[i]).
    """
    
    typename = Symbol("DTLZ", i); 
   
    @eval begin 
        @doc $help_text
        @with_kw struct $typename <: DTLZ
            n_vars :: Int = $(n_vars);
            M :: Int = $(M) ;
            k :: Int = n_vars + 1 - M;
            @assert 1 <= k;
            @assert 1 <= M <= n_vars;
            @assert n_vars == M - 1 + k;
        end
    end
end        

#%% --------------------------------------------------------------------#
# some general properties

num_vars( mop :: D where D<:DTLZ) = mop.n_vars;
num_objectives( mop :: D where D<:DTLZ) = mop.M;

for i=1:6
    typename = Symbol("DTLZ", i)
    @eval get_ideal_point( mop :: $typename ) = zeros(mop.M);
end

function constraints( mop :: D where D<:DTLZ) 
    Box( 
        zeros( num_vars( mop ) ),
        ones( num_vars( mop ) )
    )
end

function get_pareto_set( mop :: D where D<:DTLZ)
    return SamplingFunction( mop, :ParetoSet, [:random] );
end

function get_pareto_front(  mop :: D where D<:DTLZ)
    return SamplingFunction( mop, :ParetoFront, [:random] );
end


function get_points( sample_func :: SamplingFunction, mop :: Union{DTLZ1,DTLZ6,DTLZ7},
    ::Val{:ParetoSet}, num_points :: Int; method = :random )
return [ [ rand(mop.M-1); zeros(mop.k) ] for i = 1 : num_points ]
end

function get_points( sample_func :: SamplingFunction, mop :: Union{DTLZ2,DTLZ3,DTLZ4,DTLZ5},
    ::Val{:ParetoSet}, num_points :: Int; method = :random )
return [ [ rand(mop.M-1); fill(.5, mop.k) ] for i = 1 : num_points ]
end

#%% ----------------------- Objectives --------------------------#
########## DTLZ 1 ##############

# concerning `g(x)` there is a inconsistency in the source: 
# they use other factors but then the Pareto Front
# is no longer the hyperplane with Σ fⱼ = .5
function g(x, mop::Union{DTLZ1, DTLZ3})
    χ = x[mop.M : end];
    return ( 
        .75 * mop.k +
        sum( (χ .- .5 ).^2 ) -
        sum( cos.( (20 * pi) .* ( χ .- .5) ) ) 
    );
end

function get_objectives( mop :: DTLZ1 )
    objectives = Vector{Function}( undef, mop.M )
    for j = 1 : mop.M
        fⱼ = function( x :: Vector{R} where R <: Real)
            return ( 
                .5 * prod( x[1 : mop.M-j] ) *
                (1 + g( x , mop )) *
                (1 - x[mop.M-j+1])^(j>1) 
            )
        end
        objectives[j] = fⱼ
    end
    return objectives
end

########## DTLZ 2 ##############

θ(x, :: Union{DTLZ2, DTLZ3} ) = (pi/2) .* x;
function g(x, mop:: Union{DTLZ2,DTLZ4,DTLZ5} )
    χ = x[mop.M : end];
    return sum( (χ .- 0.5 ).^2 );
end

function get_objectives( mop :: Union{DTLZ2,DTLZ3,DTLZ4,DTLZ5,DTLZ6} )
    objectives = Vector{Function}( undef, mop.M )
    for j = 1 : mop.M
        fⱼ = function( x :: Vector{R} where R <: Real)
            return ( 
                (1 + g( x, mop ) ) *
                prod( cos.( θ(x,mop)[1 : mop.M - j] ) ) *
                sin( θ(x,mop)[mop.M-j+1] )^(j>1)
            )
        end
        objectives[j] = fⱼ
    end
    return objectives
end
########## DTLZ 3 ##############

# `θ` as for DTLZ2
# `g` as for DTLZ1  

# `get_objectives` defined above

########## DTLZ 4 ##############

θ(x, :: DTLZ4) = (pi/2) .* x.^100;
# `g` same as for DTLZ2 & DTLZ5
    
# `get_objectives` defined above

########## DTLZ 5 ##############

# `g` same as for DTLZ2 & DTLZ4

function θ(x, mop :: Union{DTLZ5, DTLZ6} )
    G = g( x, mop );
    return [ 
        (pi/2) * x[1];
        ( pi / ( 4 * (1 + G ) ) ) .* ( 1 .+ (2*G) .* x[2: mop.M - 1] ) ;
        x[mop.M : end] 
    ];
end

# `get_objectives` defined above

########## DTLZ 6 ##############

function g(x, mop::DTLZ6) 
    χ = x[mop.M : end];
    return sum( χ.^0.1 )
end
# `θ` same as for DTLZ5
  
# `get_objectives` defined above

########## DTLZ 7 ##############

function g(x, mop::DTLZ7)
    χ = x[mop.M : end];
    return 1 + (9/mop.k) * sum( χ )
end 

function get_objectives( mop::DTLZ7 )
    objectives = Vector{Function}( undef, mop.M )

    h(F, g) = mop.M - sum( F[1 : mop.M-1] ./ (1 + g) .* ( 1 .+ sin.( (3*pi) .* F[ 1:mop.M-1 ] ) ) );

    for j=1:mop.M-1
        fⱼ = function( x :: Vector{R} where R<:Real )
            x[j]
        end
        objectives[j] = fⱼ;
    end
    
    objectives[mop.M] = function( x :: Vector{R} where R<:Real )
        F = [ objectives[j](x) for j=1:mop.M-1 ];
        G = g( x, mop );
        return (1+G) * h(F,G)
    end
    return objectives
end
 
# TODO check this
function get_ideal_point( mop::DTLZ7)
    [zeros( mop.M - 1); 1.0];
end
##########################################################################
