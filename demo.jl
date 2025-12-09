using LinearAlgebra, Random, DelimitedFiles
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")


n = 10^2;
Random.seed!(5)
x0 = rand(n);
x0[1] = x0[2];
x0sort = sort(x0, rev=true);
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
x0prepop, verb, hist = false, false, false
for ki in 1:length(klevel)
  print(ki)
  k = max(2, Int64(ceil(n*klevel[ki])))
  tks = sum(x0sort[1:k]);
  for ri in 1:length(rlevel)
    println("ri is $(ri), ki is $(ki)")
    r = tks * float(rlevel[ri]) # 1
    xstar1 = similar(x0);
    xstar2 = similar(x0);
    @views sp, k0k1, nit, t, case, lam = project_topksum_esgs_experiment!(
      xstar1, x0sort, x0, r, k, sum(x0sort[1:k]) > r, x0prepop, verb, hist
    )
    @views flag, Ite, t = project_topksum_EIPS_experiment!(
      xstar2, x0, r, k, false, true, false
    )
    if maximum(xstar1 .- xstar2) > 1e-3
      @warn("wrong, $(maximum(xstar1 .- xstar2)), ri = $(ri), ki = $(ki)")
    end

  end
end


xstar3 = similar(x0);
k = max(2, Int64(ceil(n*klevel[13])))
tks = sum(x0sort[1:k]);
r = tks * float(rlevel[3]) # 1
@views flag, Ite, t = project_topksum_EIPS_experiment!(
  xstar3, x0, r, k, false, true, true
)