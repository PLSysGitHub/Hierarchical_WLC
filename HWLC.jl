"""
The module contains functions used for data generation and plotting
"""

module HWLC

using Distributions, Random
include("WLC.jl") #Force extension for single chains

export force_extensions, stiffnesses, power_law_samples

"""
Given the lengths and persistence lengths of WLC elements in a series,
return extensions for given forces
"""
function force_extensions(forces, lengths, persistence_lengths, kT, Marko_Siggia)
    extensions=zeros(length(forces),length(lengths))
    for i in 1:length(forces)
        extensions[i,:]=wlc_ext.(lengths,persistence_lengths,forces[i],kT,Marko_Siggia)
    end
    return vec(sum(extensions,dims=2)),extensions
end

"""
Numerical differentiation of F(x) to yield F'(F)
"""
function stiffnesses(forces, extensions)
    return map(x->num_differentiate(x,forces), eachcol(extensions))
end

"""
Helper for sampling truncated power-law distribution
"""
function invert_cdf(p,minl,maxl,α)
    return (p*(maxl^(1+α)-minl^(1+α))+minl^(1+α))^(1/(1+α))
end

"""
Return persistence lengths sampled using truncated power-law
"""
function power_law_samples(α,N,minLp,maxLp)
    p_lengths=rand(Uniform(0.,1.),N)
    p_lengths=invert_cdf.(p_lengths,minLp,maxLp,α)
    return p_lengths
end
end
