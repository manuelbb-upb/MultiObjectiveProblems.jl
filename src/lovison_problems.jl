# Problems taken from 
# A. Lovison - "Singular Continuation: Generating Piecewise 
# Linear Approximations to Pareto Sets via Global Analysis"
# (Thomann refers to these as 'lovison1-4')
# lovison1 is basically Two Parabolas.
# Note: signs are switched because we minimize
# Note: there are more problems in the paper
# TODO: optimal and critical points!

abstract type Lovison <: MOP end;

for i = 1 : 4
    typename = Symbol("Lovison", i);
    typename_constrained = Symbol( "Lovison",i, "Constrained" );
    @eval begin
        @doc $("""
                struct $typename <: Lovison <: MOP end;

            Problem taken from 
            A. Lovison - "Singular Continuation: Generating Piecewise 
            Linear Approximations to Pareto Sets via Global Analysis"
        """)
        struct $typename <: Lovison end;

        @doc $("""
                struct $typename_constrained <: Lovison <: MOP end;

            Problem taken from 
            A. Lovison - "Singular Continuation: Generating Piecewise 
            Linear Approximations to Pareto Sets via Global Analysis"
        """)
        struct $typename_constrained <: Lovison end;
        
        num_vars(::Union{$typename, $typename_constrained}) = 2;
        num_objectives(::Union{$typename, $typename_constrained}) = 2;

        export $typename, $typename_constrained;
    end
end

function get_objectives(::Lovison1)
    return [
        (x :: Vector{R} where R<:Real) -> sum( [1.05; 0.95].* x.^2 );
        (x :: Vector{R} where R<:Real) -> .99*(x[1]-3)^2 + 1.03*(x[2]-2.5)^2
    ];
end

function get_objectives(::Lovison2)
    return [
        (x :: Vector{R} where R<:Real) -> x[2];
        (x :: Vector{R} where R<:Real) -> (x[2]- x[1]^3)/(x[1]+1)
    ];
end 

function get_objectives(::Lovison3)
    return [
        (x :: Vector{R} where R<:Real) -> sum( x.^2 ); 
        (x :: Vector{R} where R<:Real) -> (x[1]-6)^2 - (x[2]+.3)^2
    ];
end 

function get_objectives(::Lovison4)
    return [
        (x :: Vector{R} where R<:Real) -> sum( x.^2 ) + 4 * exp( -(x[1]+2)^2 - x[2]^2 ) -  exp( -(x[1]-2)^2 - x[2]^2 )
        (x :: Vector{R} where R<:Real) -> (x[1]-6)^2 + (x[2]+.5)^2
    ];
end 

# box constraints as given by Thomann (p.XII)
constraints( ::Lovison1Constrained ) = Box(
    zeros(2),
    fill(3.0, 2) 
);

constraints( ::Lovison2Constrained ) = Box(
    fill(-.5, 2),
    [0.0, 0.5]
);

constaints( ::Lovison3Constrained ) = Box(
    [0.0, -4.0],
    [6.0, 4.0]
);

constaints( ::Lovison4Constrained ) = Box(
    [0.0, -1.0],
    [6.0, 1.0]
);