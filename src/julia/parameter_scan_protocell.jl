using ClusterManagers

addprocs_sge(30)

require("3_abm_protocell.jl")

using Distributions
using Iterators
using DataFrames


#N=vcat(collect(5:10),5*collect(3:20))
N=vcat(collect(5:10),5*collect(3:20),150,100*collect(2:10),1100)

#N=vcat(150)
#seqtypes=[3:25]

#mu=vcat(0,0.0001*collect(1:50),0.001*collect(6:10))
mu=vcat(0)

#delta=vcat(0.002,0.01*collect(1:10),0.1*collect(2:55))




looptuples = Iterators.product(N,mu)

#----------------------------------
# main simulation
#----------------------------------

pmap(x->SimPop(x[1],x[2]),looptuples)

