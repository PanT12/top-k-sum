global const PROJPATH = pwd();
include(pwd() * "/src/TopKSum.jl")
using .TopKSum
using Test, LinearAlgebra
using Random

n = 10^2;
Random.seed!(n)
x0 = rand(n);
sig = sortperm(x0, rev=true, alg=QuickSort);
x0sort = x0[sig];
k = Int64(ceil(n * 0.5));
r = sum(x0sort[1:k]) * 0.5; # 1

xstar_eips = similar(x0);
out_eips = project_topksum_EIPS!(xstar_eips, x0, r, k)

xstar_esgs = similar(x0);
out_esgs = project_topksum_esgs!(xstar_esgs, x0sort, r, k, true)
xstar_esgs[sig] = xstar_esgs;

xstar_pclp = similar(x0);
out_pclp = project_topksum_plcp!(xstar_pclp, x0sort, r, k, true)
xstar_pclp[sig] = xstar_pclp;

xstar_grid = similar(x0);
out_grid = project_topksum_grid!(xstar_grid, x0sort, r, k, true)
xstar_grid[sig] = xstar_grid;

xstar_gurobi = similar(x0);
out_grbs = project_topksum_guro!(xstar_gurobi, sig, x0, r, k, true)


xstar_eips_sort = sort(xstar_eips, rev=true, alg=QuickSort);
xstar_esgs_sort = sort(xstar_esgs, rev=true, alg=QuickSort);
xstar_pclp_sort = sort(xstar_pclp, rev=true, alg=QuickSort);
xstar_grid_sort = sort(xstar_grid, rev=true, alg=QuickSort);
xstar_gurobi_sort = sort(xstar_gurobi, rev=true, alg=QuickSort);

# test FeasibilityTol
@test sum(xstar_eips_sort[1:k]) ≈ r
@test sum(xstar_esgs_sort[1:k]) ≈ r
@test sum(xstar_pclp_sort[1:k]) ≈ r 
@test sum(xstar_grid_sort[1:k]) ≈ r
@test sum(xstar_gurobi_sort[1:k]) ≈ r


# xstar_gurobi = similar(x0);
# out_grbs = project_topksum_grbs(x0sort, r, k)
# norm(xstar_eips_sort .- out_grbs[:x])


# test difference
tol = 1e-4;
@test norm(xstar_eips .- xstar_gurobi, Inf) <= tol
@test norm(xstar_eips .- xstar_esgs, Inf) <= tol
@test norm(xstar_eips .- xstar_pclp, Inf) <= tol
@test norm(xstar_eips .- xstar_grid, Inf) <= tol
