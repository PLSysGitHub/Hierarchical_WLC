include("HWLC.jl")
using .HWLC
using Plots, Distributions, MAT

#Set plot parameters
gr(xlabel="F(X) (pN)", label="",ylabel="F'(X) (pN/nm)",
    grid=false,guidefont=("Arial",10),
    legend=:topleft,
    legendfont=("Arial",10),scale=:log10)

#Import experiment data for comparison
data=matread("Experimental_data/ForceStiffnessStd.mat")
mean_experiment=vec(data["k_avg"])
std_exp=vec(data["k_err"])
fs_exp=vec(data["f_ax"])

#Set parameters
kT=4.22
Marko_Siggia=true

f_min,f_max,N_f=extrema(fs_exp)...,500
forces=exp.(range(log(f_min),log(f_max),length=N_f))

N=7        #Number of component chains
L=1200.     #Total length of assembly

#Power-law distribution for persistence lengths
slope=0.86        #Slope of stiffness-force graph
# Set parameters for f_c distribution
minfc, maxfc=9.,1000.
λ=160.

"""
Show individual and total curves for one assembly
"""
fcs=rand(truncated(Exponential(λ),minfc,maxfc),N)
persistence_lengths=fcs.^(-1)*kT
#ls=ones(N)*L/N  #All equal length
ls=rand(Normal(L/N,L/N*0.1),N) #Alternatively sample from normal distr

total_ext,extensions=force_extensions(forces,ls,persistence_lengths,kT,Marko_Siggia)
total_stiffness=stiffnesses(forces,total_ext)[1]
individual_stiffnesses=stiffnesses(forces,extensions)

#Plot individual and total response
plot(forces[2:end],total_stiffness, label="Assembly");
plot!([20,200],[10^-1.5,(200/20)^slope*10^-1.5], color="black",label="Slope $slope");
for ind in individual_stiffnesses
    plot!(forces[2:end],ind, color="lightgrey")
end
display(current())
png("Plots/exp_distr_individual_responses")


"""
Sample many assemblies, and compare mean
stiffening behaviour to data
"""
n_samples=15
sim_data=zeros(length(forces)-1,n_samples)
for i in 1:n_samples
    #sample assembly
    sample_ls=rand(Normal(L/N,L/N*0.1),N) #comment if using constant lengths
    sample_ps=rand(truncated(Exponential(λ),minfc,maxfc),N).^(-1)*kT
    tot,exts=force_extensions(forces,sample_ls,sample_ps,kT,Marko_Siggia)
    sim_data[:,i]=stiffnesses(forces,tot)[1]
end

mean_sim=mean(sim_data,dims=2)
std_sim=std(sim_data,dims=2)

plot(fs_exp,mean_experiment, label="Experimental mean",legend=:topleft, color="blue");
plot!(fs_exp,mean_experiment.+std_exp, color="lightblue");
plot!(fs_exp[findall(>(0.),mean_experiment.-std_exp)],(mean_experiment.-std_exp)[findall(>(0.),mean_experiment.-std_exp)], color="lightblue");
plot!(forces[2:end],mean_sim, label="Sim mean", color="orange", ribbon=std_sim)
png("Plots/exp_distr_compare_data_sim")
