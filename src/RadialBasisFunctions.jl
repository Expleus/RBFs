module RadialBasisFunctions

using LinearAlgebra

export RBF, TPS, Gaussiana, ∇², dX, dr, drr, dX²

abstract type RBF end

function ∇²(f::RBF, ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0 ? f.LaplacianLimit : drr(f, ξ) + (1 / r) * dr(f, ξ)
end

function dX(f::RBF, ξ::Vector{<:Number}, coord::Integer)
    @assert coord > 0 "The coordinate must be a positive integer."
    r = norm(ξ)
    return r == 0 ? 0 : (ξ[coord]/r)*dr(f,ξ)
end

function dX²(f::RBF, ξ::Vector{<:Number}, coord::Integer)
    @assert coord > 0 "The coordinate must be a positive integer."
    r = norm(ξ)
    return r == 0 ? 0 : (ξ[coord]/r)^2 * drr(f,ξ) + (r^2-ξ[coord]^2)/(r^3) * dr(f,ξ)
end

struct TPS <: RBF
    B::Real
    LaplacianLimit::Real
    TPS(B::Float64) = new(B, B == 1 ? -Inf : 0)
end

function (t::TPS)(ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0.0 ? 0.0 : r^(2*t.B) * log(r)
end

function dr(t::TPS, ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0.0 ? 0.0 : r^(2 * t.B - 1) * (2 * t.B * log(r) + 1) 
end

function drr(t::TPS, ξ::Vector{<:Number})
    r = norm(ξ)
    return t.B > 1 ? r^(2 * t.B - 2) * (2 * t.B * (2* t.B - 1) * log(r) + 4 * t.B - 1) : -Inf
end

struct Gaussiana <: RBF
    ε::Real
    LaplacianLimit::Real
    Gaussiana(ε::Float64) = new(ε, -4 * ε^2)
end

function (g::Gaussiana)(ξ::Vector{<:Number})
    r = norm(ξ)
    return exp(-(g.ε * r)^2)
end

function dr(g::Gaussiana, ξ::Vector{<:Number})
    r = norm(ξ)
    return -2 * g.ε^2 * r * g(ξ)
end

function drr(g::Gaussiana, ξ::Vector{<:Number})
    r = norm(ξ)
    return (1 - 2 * (g.ε * r)^2) * dr(g, ξ) / r
end

end # module
