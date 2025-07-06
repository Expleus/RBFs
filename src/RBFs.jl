module RBFs

using LinearAlgebra

export RBF, TPS, Gaussiana, ∇², dr, drr

abstract type RBF end

function ∇²(f::RBF, ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0 ? f.LaplacianLimit : drr(f, ξ) + (1 / r) * dr(f, ξ)
end

struct TPS <: RBF
    LaplacianLimit::Float64
    TPS() = new(0)
end

function (t::TPS)(ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0.0 ? 0.0 : r^2 * log(r)
end

function dr(t::TPS, ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0.0 ? 0.0 : 2 * r * log(r) + r
end

function drr(t::TPS, ξ::Vector{<:Number})
    r = norm(ξ)
    return 2 * log(r) + 3
end

struct Gaussiana <: RBF
    ε::Float64
    LaplacianLimit::Float64
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

struct TPSo2 <: RBF
    LaplacianLimit::Float64
    TPS() = new(0)
end

function (t::TPSo2)(ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0.0 ? 0.0 : r^4 * log(r)
end

function dr(t::TPSo2, ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0.0 ? 0.0 : r^3*(4*log(r)+1)
end

function drr(t::TPS, ξ::Vector{<:Number})
    r = norm(ξ)
    return r == 0.0 ? 0.0 : r^2*(12*log(r)+7)
end

end # module
