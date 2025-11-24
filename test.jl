using LinearAlgebra, Random, DelimitedFiles
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")


n = 10^6;
Random.seed!(n)
x0 = rand(n);
x0sort = sort(x0, rev=true);
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
for ki in 1:length(klevel)
  k = max(2, Int64(ceil(n*klevel[ki])))
  tks = sum(x0sort[1:k]);
  for ri in 1:length(rlevel)
    # if (rlevel[ri], klevel[ki]) ∉ expers
    #     continue
    # end
    println("ri is $(ri), ki is $(ki)")
    r = tks * float(rlevel[ri]) # 1
    xstarf = similar(x0);
    # @views flag1, case, Ite, k0k1, t = project_topksum_IPESF_experiment5!(xstarf, x0, r, k, false, false)
    project_topksum_EIPS_experiment!(xstarf, x0, r, k, false, true, false)
    project_topksum_EIPS_experiment!(xstarf, x0, r, k, true, true, false)
    # xbarsort = similar(x0)
    # project_topksum_plcp_experiment!(xbarsort, sig, x0, r, k, tks>r, false, false, false)
    # flag[ri,ki]=flag1
    # if maximum(xstarf .- xbarsort) > 1e-3
    #   @warn("wrong, $(maximum(xstarf .- xbarsort))")
    # end

  end
end


k = max(2, Int64(ceil(n*klevel[2])));
tks = sum(x0sort[1:k]);
r = tks * float(rlevel[2]); # 1
xstarf = similar(x0);
xstarf2 = similar(x0);
project_topksum_EIPS_experiment!(xstarf, x0, r, k, false, true, true)
project_topksum_EIPS_experiment!(xstarf2, x0, r, k, true, true, true)


xstar = similar(x0);
project_topksum_esgs_experiment!(
  xstar, x0sort, x0, r, k, true, false, false, false
)
u = 0.20727306697516568;
m = sum(x0sort .>= u) # number of elements larger than u
b = view(x0, x0 .< u) # elements larger than u
lf = (r - sum(max.(x0 .- u, 0))) / k
project_topksum_g_searching1(b, u, k-m)


x = 0.1:0.001:1.0;
f = zeros(length(x));
g = zeros(length(x));
for i in 1:length(x)
  f[i] = (r - sum(max.(x0 .- x[i], 0))) / k
  m = sum(x0sort .>= x[i]);
  if m > k
    continue
  end
  println(i)
  g[i], _ = project_topksum_g_searching(x0, f[i], x[i], k-m)
end
using Plots
plot(x, f, label="f", xlabel="u", ylabel="lf", title="f vs u")
plot!(x, g, label="g", xlabel="u", ylabel="lf", title="g vs u")




include(pwd() * "/src/TopKSum.jl")
using .TopKSum
using Test, LinearAlgebra
using Random






n = 10^3;
Random.seed!(n)
x0 = rand(n);
sig = sortperm(x0, rev=true, alg=QuickSort);
x0sort = x0[sig];
k = 500;
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
@test norm(xstar_eips .- xstar_gurobi) <= tol
@test norm(xstar_eips .- xstar_esgs, Inf) <= tol
@test norm(xstar_eips .- xstar_pclp, Inf) <= tol
@test norm(xstar_eips .- xstar_grid, Inf) <= tol




norm(xstar_eips .- x0)
norm(xstar_esgs .- x0)
norm(xstar_pclp .- x0)
norm(xstar_grid .- x0)
norm(xstar_gurobi .- x0)




diff = xstar_eips_sort .- xstar_gurobi_sort;
argmax(diff)
argmin(diff)

