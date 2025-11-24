function initialize_output(n, nr, nk)
  alg_names = ["esgs", "eips"]
  out = Dict()
  out["esgs"] = Dict()
  out["eips"] = Dict()
  out["ALL"] = Dict()
  
  # t_init
  out["eips"][:t_init] = fill(+Inf, nr, nk)
  
  # t_run
  out["eips"][:t_run] = fill(+Inf, nr, nk)
  
  # t_primal: primal recovery
  out["eips"][:t_primal] = fill(+Inf, nr, nk)
  
  # obj
  out["eips"][:obj] = fill(+Inf, nr, nk) # obj value 1
  
  # nit
  out["eips"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 
  out["eips"][:nit_p] = zeros(Int64, nr, nk) # number of pivot iterations 
  out["eips"][:nit_e] = zeros(Int64, nr, nk) # number of exact iterations 
  out["eips"][:nit_total] = zeros(Int64, nr, nk) # number of total iterations 

  
  # case
  out["eips"][:case] = zeros(Int64, nr, nk) # case id: 1 = strict, 2 = strict but lcp ez; 0 = inactive; -1 = k==1; -2 = k==n
  
  # bestfeas
  out["eips"][:bestfeas] = zeros(Bool, nr, nk) # does 1 satisfy Mk(xbar) <= r and have lowest objval

  # other
  out["ALL"][:obj_1234] = Array{Vector{Int64}}(undef, nr, nk) # best objective value
  out["ALL"][:active] = zeros(Bool, nr, nk) # is Mk(x0) > r?
  out["ALL"][:x0sort] = zeros(Float64, n) # x0sort
  out["ALL"][:x0] = zeros(Float64, n) # x0
  out["ALL"][:t_psort] = fill(+Inf, nk) # partial sort time

  
  # solution
  out["eips"][:xstar] = zeros(Float64, n)
  return out
  end


function time_test(n::Integer, nrep::Integer,
  rlevel::Vector, klevel::Vector,
  datapath::String,
  distribution::String="Uniform",
  expers=nothing,
  rep_offset=0
  )
  # setup
  nr = length(rlevel)
  nk = length(klevel)
  alg_names = ["without_init", "with_init"]
  out = initialize_output(n, nr, nk)

  hist = false # don't record history

  # test output
  println(datapath * "/$(distribution)/$n/")
  mkpath(datapath * "/$(distribution)/$n/")
  writedlm(datapath*"/$(distribution)/$(n)/rlevel.csv", rlevel)
  writedlm(datapath*"/$(distribution)/$(n)/klevel.csv", klevel)
  

  for repi in 1:nrep
    println("rep progress: $(round((repi-rep_offset) / nrep, digits=2))")
    total = length(klevel) * length(rlevel)
    flag = 0
    
    # generate x0
    if distribution == "Normal"
      x0 = randn(n);
    elseif distribution == "Uniform"
      x0 = rand(n);
    end

    out["ALL"][:x0] .= x0;
    # writedlm(datapath*"/$(distribution)/$(n)/x0_$repi.csv", out["ALL"][:x0])
    xstar = similar(x0);
    out["ALL"][:t_sort] = @elapsed begin
      x0sort = sort(x0, rev=true, alg=QuickSort)
    end
    out["ALL"][:x0sort] .= x0sort;
    writedlm(datapath*"/$(distribution)/$(n)/t_sort_$repi.csv", out["ALL"][:t_sort])
  
    # loop over k and r 
    for ki in eachindex(klevel)
      k = max(2, Int64(ceil(n*klevel[ki])))
      tks = sum(x0sort[1:k])
      out["ALL"][:t_psort][ki] = @elapsed begin
        x0psort = partialsort(x0, 1:k);
      end
      writedlm(datapath*"/$(distribution)/$(n)/t_psort_$repi.csv", out["ALL"][:t_psort])
      
      # GC.enable(false)

      for ri in eachindex(rlevel)
        # println(ri)
        r = Float64(rlevel[ri]) * tks
        if !isnothing(expers)
          if (rlevel[ri], klevel[ki]) ∉ expers
          # println("  --> skipping: $((rlevel[ri], klevel[ki]))")
          continue
          end
        end
        # println("$((rlevel[ri], klevel[ki]))")
        out["ALL"][:active][ri,ki] = sum(x0sort[1:k]) > r
        # eips
        xstar = similar(x0);
        debug = false;
        @views flag, Ite, t = project_topksum_EIPS_experiment!(
          xstar, x0, r, k, false, debug, hist
        )
        D = out["eips"]
        D[:case][ri,ki] = flag
        D[:xstar] .= xstar
        D[:name] = "eips"
        D[:t_init][ri,ki] = t[1]
        D[:t_run][ri,ki] = t[2]
        D[:t_primal][ri,ki] = t[3]
        # D[:tks][ri,ki] = max(sum(sort(xstar, rev=true)[1:k]) - r, 0.0)
        D[:obj][ri,ki] = 0.5 * norm(xstar .- x0, 2)^2
        D[:nit_init][ri,ki] = Ite[1]
        D[:nit_p][ri,ki] = Ite[2]
        D[:nit_e][ri,ki] = Ite[3]
        D[:nit_total][ri,ki] = sum(Ite)
        # D[:case][ri,ki] = case
        writeout_piv(D, n, repi, datapath, distribution)

      end
      # GC.enable(true)

    end
  end
end

function writeout_piv(D::Dict, n::Integer, repi::Integer, datapath::String, distribution::String)
  id=1
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_flag_$repi.csv", D[:case])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_total_$repi.csv", D[:t_init] + D[:t_run] + D[:t_primal])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_run_$repi.csv", D[:t_run])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_primal_$repi.csv", D[:t_primal])
  # writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_tks_$repi.csv", D[:tks])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_init_$repi.csv", D[:t_init])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_obj_$repi.csv", D[:obj])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_init_$repi.csv", D[:nit_init])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_p_$repi.csv", D[:nit_p])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_e_$repi.csv", D[:nit_e])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_total_$repi.csv", D[:nit_total])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_case_$repi.csv", D[:case])
end



using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")
include(PROJPATH * "/experiment/data_process.jl")
global const DATAPATH = PROJPATH * "/complexity"
# global const DATAPATH = PROJPATH * "/complexity_tau_r_near1"
# global const DATAPATH = PROJPATH * "/complexity_tau_r=1"


nlevel = 10 .^(collect(3:7));
nlevel = [10^3; 10^7];
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
ri = [2; 6; 10; 13; 15; 17];
ki = [1; 6; 8; 12];

klevel_sub = klevel[ki];
rlevel_sub = rlevel[ri];

expers = [(x, y) for x in rlevel_sub for y in klevel_sub]

repeat = 100;
rep_offset = 0;
distribution = "Uniform";
for ni in eachindex(nlevel)
  local n = nlevel[ni]
  Random.seed!(n)
  println("===================")
  println("   n=$n  ")
  println("===================")
  time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, expers)
end

root_dir = DATAPATH * "/Uniform";
rm(root_dir * "/1000", recursive=true)

out = load_results(DATAPATH, 1, 100, 1)
process_load_result(out, DATAPATH, 1)

####################################################################################

# # near 1 part 
# ri = [13]
# klevel = [1; 2; 3; 10; 100; 500; 1000; 2000; 3000; 4000; 5000] .// 10^4
# ki = [1; 2; 3; 4; 5]

global const DATAPATH = PROJPATH * "/complexity_tau_r=1"
# = 1 part 
ki = [1; 2; 3; 4; 5; 6; 7; 8; 9];
ri = [13];
nlevel = Int.([1; 10; 50; 100; 500; 1000].*10^5)

klevel_sub = klevel[ki];
rlevel_sub = rlevel[ri];

expers = [(x, y) for x in rlevel_sub for y in klevel_sub]

repeat = 100;
rep_offset = 0;
distribution = "Uniform";
for ni in eachindex(nlevel)
  local n = nlevel[ni]
  Random.seed!(n)
  println("===================")
  println("   n=$n  ")
  println("===================")
  time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, expers)
end

root_dir = DATAPATH * "/Uniform";
rm(root_dir * "/100000", recursive=true)

out = load_results(DATAPATH, 1, 100, 1)

process_load_result(out, DATAPATH, 1)