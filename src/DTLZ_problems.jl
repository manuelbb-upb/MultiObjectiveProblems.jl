# Problems taken from 
# Deb et. al - "Scalable Test Problems for Evolutionary Multiobjective Optimization"

# DTLZ8 & DTLZ9 not supported yet due to constraints...
export DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7, DTLZ;
abstract type DTLZ <: MOP end;

typehelp(i) = """ 
        struct DTLZ$(i) <: DTLZ <: MOP
            n_vars :: Int;
            M :: Int;
            k :: Int = n_vars + 1 - M;
        end

    Test problem defined in 
    Deb et. al - "Scalable Test Problems for Evolutionary Multiobjective Optimization".
    Has `n_vars` decision variables and `M` objectives.
    The decision vector `x` is split into two parts: `x[1:M-1]` and `x[M:end]`.
    The latter has length `k`.
""";

@doc typehelp(1)
@with_kw struct DTLZ1 <: DTLZ
    n_vars :: Int = 10;
    
    M :: Int = 6;   # num objfs
    k :: Int = n_vars + 1 - M;   # size of last variable group
    @assert 1 <= k;
    @assert 1 <= M <= n_vars;
    @assert n_vars == M - 1 + k;
end;

@doc typehelp(7)
@with_kw struct DTLZ7 <: DTLZ
    n_vars :: Int = 22;
    
    M :: Int = 3;
    k :: Int = n_vars + 1 - M;   # size of last variable group
    @assert 1 <= k;
    @assert 1 <= M <= n_vars;
    @assert n_vars == M - 1 + k;
end;

for i = 1 : 7 
    typename = Symbol("DTLZ", i)
    thishelp = typehelp(i)    
    if i ∈ [2,3,4,5,6]
        @eval begin
            @doc $thishelp
            @with_kw struct $typename <: DTLZ
                n_vars :: Int = 15;
                
                M :: Int = 6;   # num objfs
                k :: Int = n_vars + 1 - M;   # size of last variable group
                @assert 1 <= k;
                @assert 1 <= M <= n_vars;
                @assert n_vars == M - 1 + k;
            end
        end
        if i ∉ [1,6,7]
            @eval begin
                function get_points( sample_func :: SamplingFunction, mop :: $typename,
                        ::Val{:ParetoSet}, num_points :: Int; method = :random )
                    return [ [ rand(mop.M-1); fill(.5, mop.k) ] for i = 1 : num_points ]
                end
            end
        end
    end
    @eval begin                
        function get_pareto_set( mop :: $typename )
            return SamplingFunction( mop, :ParetoSet, [:random] );
        end

        function get_pareto_front(  mop :: $typename )
            return SamplingFunction( mop, :ParetoFront, [:random] );
        end
 
        num_vars( mop :: $typename ) = mop.n_vars;
        num_objectives( mop :: $typename ) = mop.M;

        function constraints( mop :: $typename ) 
            Box( 
                zeros( num_vars( mop ) ),
                ones( num_vars( mop ) )
            )
        end       
    end
end


function g(x, mop::Union{DTLZ1, DTLZ3})
    χ = x[mop.M : end];
    return ( 
        .75 * mop.k +
        sum( (χ .- .5 ).^2 ) -
        sum( cos.( (20 * pi) .* ( χ .- .5) ) ) 
    );
end

function g(x, mop::Union{DTLZ2,DTLZ4,DTLZ5} ) 
    χ = x[mop.M : end];
    return sum( (χ .- 0.5 ).^2 );
end
function g(x, mop::DTLZ6) 
    χ = x[mop.M : end];
    return sum( χ.^0.1 )
end

function g(x, mop::DTLZ7)
    χ = x[mop.M : end];
    return 1 + (9/mop.k) * sum( χ )
end 

θ(x, :: Union{DTLZ2,DTLZ3} ) = (pi/2) .* x;
θ(x, :: DTLZ4) = (pi/2) .* x.^100;
function θ(x, mop :: Union{DTLZ5,DTLZ6} )
    G = g( x, mop );
    return [ 
        (pi/2) * x[1];
        ( pi / ( 4 * (1 + G ) ) ) .* ( 1 .+ (2*G) .* x[2: mop.M - 1] ) ;
        x[mop.M : end] 
    ];
end

function get_objectives( mop :: DTLZ1 )
    objectives = Vector{Function}( undef, mop.M )
    
    # concerning `g(x)` there is a inconsistency in the source: 
    # they use other factors but then the Pareto Front
    # is no longer the hyperplane with Σ fⱼ = .5

    for j = 1 : mop.M
        fⱼ = function( x :: Vector{R} where R <: Real)
            return ( .5 * prod( x[1 : mop.M-j] ) *
                (1 + g( x , mop )) *
                (1 - x[mop.M-j+1])^(j>1) )
        end
        objectives[j] = fⱼ
    end
    return objectives
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

function get_points( sample_func :: SamplingFunction, mop :: Union{DTLZ1,DTLZ6,DTLZ7},
        ::Val{:ParetoSet}, num_points :: Int; method = :random )
    return [ [ rand(mop.M-1); zeros(mop.k) ] for i = 1 : num_points ]
end
