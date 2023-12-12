using ClusterManagers
using Distributed

nprocs = parse(Int, ENV["SLURM_NPROCS"])
addprocs(SlurmManager(nprocs))
