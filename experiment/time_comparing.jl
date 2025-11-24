function initialize_output(n, nr, nk)
  alg_names = ["EIPS","ESGS", "PLCP", "GRID", "GURO"]
  out = Dict()
  out["EIPS"] = Dict()
  out["ESGS"] = Dict()
  out["PLCP"] = Dict()
  out["GRID"] = Dict()
  out["GURO"] = Dict()
  out["ALL"] = Dict()
  
  # t_init
  out["EIPS"][:t_init] = fill(+Inf, nr, nk)
  out["ESGS"][:t_init] = fill(+Inf, nr, nk)
  out["PLCP"][:t_init] = fill(+Inf, nr, nk)
  out["GRID"][:t_init] = fill(+Inf, nr, nk)
  out["GURO"][:t_init] = fill(+Inf, nr, nk)
  
  # t_run
  out["EIPS"][:t_run] = fill(+Inf, nr, nk)
  out["ESGS"][:t_run] = fill(+Inf, nr, nk)
  out["PLCP"][:t_run] = fill(+Inf, nr, nk)
  out["GRID"][:t_run] = fill(+Inf, nr, nk)
  out["GURO"][:t_run] = fill(+Inf, nr, nk)
  
  # t_primal: primal recovery
  out["EIPS"][:t_primal] = fill(+Inf, nr, nk)
  out["ESGS"][:t_primal] = fill(+Inf, nr, nk)
  out["PLCP"][:t_primal] = fill(+Inf, nr, nk)
  out["GRID"][:t_primal] = fill(+Inf, nr, nk)
  out["GURO"][:t_primal] = fill(0.0, nr, nk)
  # out["IPES"][:flag] = fill(0, nr, nk)

  # t_sort
  out["ALL"][:t_sort] = +Inf # full sort time
  out["ALL"][:t_psort] = fill(+Inf, nk) # partial sort time
  
  # tks infeasibility
  out["EIPS"][:tks] = fill(+Inf, nr, nk) # method 1 violation of Mk(xbar) <= r
  out["ESGS"][:tks] = fill(+Inf, nr, nk) # method 2 violation of Mk(xbar) <= r
  out["PLCP"][:tks] = fill(+Inf, nr, nk) # method 3 violation of Mk(xbar) <= r
  out["GRID"][:tks] = fill(+Inf, nr, nk) # method 4 violation of Mk(xbar) <= r
  out["GURO"][:tks] = fill(+Inf, nr, nk) # method 5 violation of Mk(xbar) <= r

  # obj
  out["EIPS"][:obj] = fill(+Inf, nr, nk) # obj value 1
  out["ESGS"][:obj] = fill(+Inf, nr, nk) # obj value 2
  out["PLCP"][:obj] = fill(+Inf, nr, nk) # obj value 3
  out["GRID"][:obj] = fill(+Inf, nr, nk) # obj value 4
  out["GURO"][:obj] = fill(+Inf, nr, nk) # obj value 5
  

  # nit init
  out["EIPS"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 1
  out["ESGS"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 2
  out["PLCP"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 3
  out["GRID"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 4
  out["GURO"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 5

  # nit pivot
  out["EIPS"][:nit_p] = zeros(Int64, nr, nk) # number of pivot iterations 1
  out["ESGS"][:nit_p] = zeros(Int64, nr, nk) # number of pivot iterations 2
  out["PLCP"][:nit_p] = zeros(Int64, nr, nk) # number of pivot iterations 3
  out["GRID"][:nit_p] = zeros(Int64, nr, nk) # number of pivot iterations 4
  out["GURO"][:nit_p] = zeros(Int64, nr, nk) # number of pivot iterations 5  

  # nit exact
  out["EIPS"][:nit_e] = zeros(Int64, nr, nk) # number of exact iterations 1
  out["ESGS"][:nit_e] = zeros(Int64, nr, nk) # number of exact iterations 2
  out["PLCP"][:nit_e] = zeros(Int64, nr, nk) # number of exact iterations 3
  out["GRID"][:nit_e] = zeros(Int64, nr, nk) # number of exact iterations 4
  out["GURO"][:nit_e] = zeros(Int64, nr, nk) # number of exact iterations 5    

  # nit total
  out["EIPS"][:nit_total] = zeros(Int64, nr, nk) # number of total iterations 1
  out["ESGS"][:nit_total] = zeros(Int64, nr, nk) # number of total iterations 2
  out["PLCP"][:nit_total] = zeros(Int64, nr, nk) # number of total iterations 3
  out["GRID"][:nit_total] = zeros(Int64, nr, nk) # number of total iterations 4
  out["GURO"][:nit_total] = zeros(Int64, nr, nk) # number of total iterations 5    
  
  # case
  out["EIPS"][:case] = zeros(Int64, nr, nk) # case id: 1 = strict, 2 = strict but lcp ez; 0 = inactive; -1 = k==1; -2 = k==n
  out["ESGS"][:case] = zeros(Int64, nr, nk) # case id: 0: inactive, 1: active
  out["PLCP"][:case] = zeros(Int64, nr, nk) # case id: 0: inactive, 1: active
  out["GRID"][:case] = zeros(Int64, nr, nk) # case id: 0: inactive, 1: active
  out["GURO"][:case] = zeros(Int64, nr, nk) # case id: 0: inactive, 1: active
  
  # bestfeas
  out["EIPS"][:bestfeas] = zeros(Bool, nr, nk) # does 1 satisfy Mk(xbar) <= r and have lowest objval
  out["ESGS"][:bestfeas] = zeros(Bool, nr, nk) # does 2 satisfy Mk(xbar) <= r and have lowest objval
  out["PLCP"][:bestfeas] = zeros(Bool, nr, nk) # does 3 satisfy Mk(xbar) <= r and have lowest objval
  out["GRID"][:bestfeas] = zeros(Bool, nr, nk) # does 4 satisfy Mk(xbar) <= r and have lowest objval
  out["GURO"][:bestfeas] = zeros(Bool, nr, nk) # does 5 satisfy Mk(xbar) <= r and have lowest objval

  # other
  out["ALL"][:obj_1234] = Array{Vector{Int64}}(undef, nr, nk) # best objective value
  out["ALL"][:active] = zeros(Bool, nr, nk) # is Mk(x0) > r?
  out["ALL"][:x0sort] = zeros(Float64, n) # x0sort
  out["ALL"][:x0] = zeros(Float64, n) # x0
  out["ALL"][:t_psort] = fill(+Inf, nk) # partial sort time
  out["ALL"][:sig_time] = fill(+Inf, nk) # permutation time
  out["ALL"][:sig] = zeros(Float64, n) # permutation

  
  # solution
  out["EIPS"][:xstar] = zeros(Float64, n)
  out["ESGS"][:xstar] = zeros(Float64, n)
  out["PLCP"][:xstar] = zeros(Float64, n)
  out["GRID"][:xstar] = zeros(Float64, n)
  out["GURO"][:xstar] = zeros(Float64, n)
  return out
  end


function time_test(n::Integer, nrep::Integer,
  rlevel::Vector, klevel::Vector,
  datapath::String,
  distribution::String="Uniform",
  maxn_gurobi::Integer=1000_000,
  maxn_grid::Integer=100_000,
  expers=nothing,
  rep_offset=0,
  alg1::Bool=true, alg2::Bool=true, alg3::Bool=true, alg4::Bool=true, alg5::Bool=true, 
  )
  # setup
  nr = length(rlevel)
  nk = length(klevel)
  alg_names = ["EIPS", "ESGS", "PLCP", "GRID", "GURO"]
  out = initialize_output(n, nr, nk)

  x0prepop = false # prepopulate xbarsort .= x0sort
  hist = false # don't record history
  verb = false # no console output
  

  # test output
  println(datapath * "/$(distribution)/$n/")
  mkpath(datapath * "/$(distribution)/$n/")
  writedlm(datapath*"/$(distribution)/$(n)/rlevel.csv", rlevel)
  writedlm(datapath*"/$(distribution)/$(n)/klevel.csv", klevel)
  writedlm(datapath*"/$(distribution)/$(n)/maxn_gurobi.csv", maxn_gurobi)
  writedlm(datapath*"/$(distribution)/$(n)/maxn_grid.csv", maxn_grid)

  # sim_repi(repi::Int64, expers) = begin # modify `out` in parallel
  for repi in 1:nrep
    println("rep progress: $(round((repi-rep_offset) / nrep, digits=2))")
    total = length(klevel) * length(rlevel)
    remaining = total
    completed = 0
    printpct = 5
    flag = 0
    
    # generate x0
    if distribution == "Normal"
      x0 = randn(n);
    elseif distribution == "Uniform"
      x0 = rand(n);
    end

    # GC.enable(false)

    out["ALL"][:x0] .= x0;
    # writedlm(datapath*"/$(distribution)/$(n)/x0_$repi.csv", out["ALL"][:x0])
    xstar = similar(x0);
    out["ALL"][:t_sort] = @elapsed begin
      x0sort = sort(x0, rev=true, alg=QuickSort)
    end
    out["ALL"][:sig_time] = @elapsed begin
      sig = sortperm(x0, rev=true, alg=QuickSort)
    end
    out["ALL"][:sig] = sig
    out["ALL"][:x0sort] .= x0sort;
    writedlm(datapath*"/$(distribution)/$(n)/t_sort_$repi.csv", out["ALL"][:t_sort])

    # GC.enable(true)
  
    # GC.enable(false)
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
        # GC.enable(false)
        r = Float64(rlevel[ri]) * tks
        if !isnothing(expers)
          if (rlevel[ri], klevel[ki]) ∉ expers
          # println("  --> skipping: $((rlevel[ri], klevel[ki]))")
          continue
          end
        end
        # println("$((rlevel[ri], klevel[ki]))")
        out["ALL"][:active][ri,ki] = sum(x0sort[1:k]) > r

        # GC.enable(false)

        # EIPS
        xstar = similar(x0);
        debug = false;
        @views flag, Ite, t = project_topksum_EIPS_experiment!(
          xstar, x0, r, k, false, debug, hist
        )
        D = out["EIPS"]
        D[:case][ri,ki] = flag
        D[:xstar] .= xstar
        D[:name] = "EIPS"
        D[:t_init][ri,ki] = t[1]
        D[:t_run][ri,ki] = t[2]
        D[:t_primal][ri,ki] = t[3]
        D[:tks][ri,ki] = max(sum(sort(D[:xstar], rev=true)[1:k]) - r, 0.0) #! is unsorted
        # D[:tks][ri,ki] = max(sum(sort(xstar, rev=true)[1:k]) - r, 0.0)
        D[:obj][ri,ki] = 0.5 * norm(xstar .- x0, 2)^2
        D[:nit_init][ri,ki] = Ite[1]
        D[:nit_p][ri,ki] = Ite[2]
        D[:nit_e][ri,ki] = Ite[3]
        D[:nit_total][ri,ki] = sum(Ite)
        if alg1
          writeout_piv(D, n, repi, datapath, distribution)
        end
        
        # GC.gc()
        # GC.enable(true)

        # GC.enable(false)

        # ESGS
        xstar = similar(x0sort);
        @views sp, k0k1, nit, t, case, lam = project_topksum_esgs_experiment!(
          xstar, x0sort, x0, r, k, sum(x0sort[1:k]) > r, x0prepop, verb, hist
        )
        time_sort = out["ALL"][:active][ri,ki] ? out["ALL"][:t_sort] + out["ALL"][:t_psort][ki] : out["ALL"][:t_psort][ki];
        D = out["ESGS"]
        D[:xstar] .= xstar
        D[:name] = "ESGS"
        # D[:k0][ri,ki] = k0k1[1]
        # D[:k1][ri,ki] = k0k1[2]
        D[:t_init][ri,ki] = t[1] + time_sort;
        D[:t_run][ri,ki] = t[2]
        D[:t_primal][ri,ki] = t[3] #+ time_convert
        D[:tks][ri,ki] = max(sum(sort(xstar, rev=true)[1:k]) - r, 0.0) #! is unsorted
        # D[:ord][ri,ki] = maximum(max.(diff(D[:xstar]), 0.0))
        # D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
        D[:obj][ri,ki] = 0.5 * norm(xstar .- x0, 2)^2
        D[:nit_total][ri,ki] = nit
        D[:case][ri,ki] = case
        if alg2
          writeout_piv(D, n, repi, datapath, distribution)
        end

        # GC.gc()
        # GC.enable(true)

        # GC.enable(false)

        # PLCP need to record permutation
        xbarsort = similar(x0sort)
        @views sp, ab, nit, t, case = project_topksum_plcp_experiment!(
          xbarsort, out["ALL"][:sig], x0, r, k, out["ALL"][:active][ri,ki], x0prepop, verb, hist
        )
        time_sort = out["ALL"][:active][ri,ki] ? out["ALL"][:sig_time] + out["ALL"][:t_psort][ki] : out["ALL"][:t_psort][ki];
        D = out["PLCP"]
        D[:xstar] .= xbarsort
        # k0k1 = get_k0k1(sort(D[:xbarsort], rev=true), k) #! soln not sorted but still can recover indices
        a = ab[1]
        b = ab[2]
        if a >=0 && b > 0
        a = max(a-1, 0)
        b = min(b+1, n)
        elseif case == 2
        a = k-1
        b = k
        end
        D[:name] = "PLCP"
        # D[:k0][ri,ki] = a
        # D[:k1][ri,ki] = b
        D[:t_init][ri,ki] = t[1] + time_sort
        D[:t_run][ri,ki] = t[2]
        D[:t_primal][ri,ki] = t[3]
        D[:tks][ri,ki] = max(sum(D[:xstar][1:k]) - r, 0.0) #! is sorted
        # D[:ord][ri,ki] = maximum(max.(diff(D[:xbarsort]), 0.0))
        # D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
        D[:obj][ri,ki] = 0.5 * norm(D[:xstar] .- out["ALL"][:x0], 2)^2
        D[:nit_total][ri,ki] = nit
        D[:case][ri,ki] = case
        if alg3
          writeout_piv(D, n, repi, datapath, distribution)
        end

        # GC.gc()
        # GC.enable(true)

        # GC.enable(false)

        # GRID # no need to record permutation
        if n < maxn_grid # && [ri, ki] ∉ specific_case
          xbarsort = similar(x0sort)
          @views sp, k0k1, nit, t, case, grid_solved = project_topksum_grid_experiment!(
            xbarsort, x0, x0sort, r, k, out["ALL"][:active][ri,ki], x0prepop, verb, hist, 3_000
          )
          time_sort = out["ALL"][:active][ri,ki] ? out["ALL"][:t_sort] + out["ALL"][:t_psort][ki] : out["ALL"][:t_psort][ki];
          D = out["GRID"]
          D[:xstar] .= xbarsort
          D[:name] = "GRID"
          # D[:k0][ri,ki] = k0k1[1]
          # D[:k1][ri,ki] = k0k1[2]
          D[:t_init][ri,ki] = t[1] + time_sort
          D[:t_run][ri,ki] = t[2]
          D[:t_primal][ri,ki] = t[3]
          D[:tks][ri,ki] = max(sum(sort(D[:xstar], rev=true)[1:k]) - r, 0.0) #! not necessarily sorted
          # D[:ord][ri,ki] = maximum(max.(diff(D[:xbarsort]), 0.0))
          # D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
          D[:obj][ri,ki] = 0.5 * norm(D[:xstar] .- out["ALL"][:x0], 2)^2
          D[:nit_total][ri,ki] = nit
          D[:case][ri,ki] = case
          if alg4
            writeout_piv(D, n, repi, datapath, distribution)
          end
          if !grid_solved
            D[:xstar] .= out["EIPS"][:xstar]
          end
        else
          D = out["GRID"]
          D[:xstar] .= out["EIPS"][:xstar]
          D[:name] = "GRID"
        end

        # GC.gc()
        # GC.enable(true)

        # GC.enable(false)

        # GURO need to record permutation
        if n < maxn_gurobi
          xbarsort .= x0sort
          @views case, nit, t, guro_solved = project_topksum_guro_experiment!(
            xbarsort, out["ALL"][:sig], x0, r, k, out["ALL"][:active][ri,ki]
          )
          time_sort = out["ALL"][:active][ri,ki] ? out["ALL"][:sig_time] + out["ALL"][:t_psort][ki] : out["ALL"][:t_psort][ki];
          D = out["GURO"]
          D[:name] = "GURO"
          D[:xstar] .= xbarsort
          # k0k1 = get_k0k1(sort(D[:xbarsort], rev=true), k)
          # D[:k0][ri,ki] = k0k1[1]
          # D[:k1][ri,ki] = k0k1[2]
          D[:t_init][ri,ki] = t[1] + time_sort
          D[:t_run][ri,ki] = t[2]
          D[:t_primal][ri,ki] = t[3]
          D[:tks][ri,ki] = max(sum(sort(D[:xstar], rev=true)[1:k]) - r, 0.0) #! not necessarily sorted
          # kappa = res[:kappa]
          # _k_ = res[:_k_]
          # D[:ord][ri,ki] = max(
          #   max(0.0, maximum(res[:x][_k_] .- res[:x][kappa])),
          #   max(0.0, maximum(res[:x][setdiff(collect(1:n), kappa)] .- res[:x][_k_])),
          # )
          # D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
          D[:obj][ri,ki] = 0.5 * norm(D[:xstar] .- out["ALL"][:x0], 2)^2
          D[:nit_total][ri,ki] = nit
          # D[:nactivecon][ri,ki] = res[:nactivecon]
          if alg5
            writeout_piv(D, n, repi, datapath, distribution)
          end
          if !guro_solved
            D[:xstar] .= out["EIPS"][:xstar]
          end
        else
          D = out["GURO"]
          D[:xstar] .= out["EIPS"][:xstar]
          D[:name] = "GURO"
        end

        # GC.gc()
        # GC.enable(true)

        # check
        diff_12 = norm(out["EIPS"][:xstar] .- out["ESGS"][:xstar], Inf)
        diff_13 = norm(out["EIPS"][:xstar] .- out["PLCP"][:xstar], Inf)
        diff_14 = norm(out["EIPS"][:xstar] .- out["GRID"][:xstar], Inf)
        diff_15 = norm(out["EIPS"][:xstar] .- out["GURO"][:xstar], Inf)
        diffs = [diff_12, diff_13, diff_14, diff_15]
        tol = 1e-3
        if argmax(diffs) == 1 && maximum(diffs) > tol
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in EIPS v ESGS)? tks infeas = $(
            round(max(0, sum(out["ESGS"][:xstar][1:k])-r),digits=5)
          ), diff_12 = $(
            round(diff_12,digits=8)
          )")
        elseif argmax(diffs) == 2 && maximum(diffs) > tol && n <= maxn_grid
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in EIPS v PLCP)? tks infeas = $(
            round(max(0, sum(out["PLCP"][:xstar][1:k])-r),digits=5)
          ), diff_13 = $(
            round(diff_13,digits=8)
          )")
        elseif argmax(diffs) == 3 && maximum(diffs) > tol && n <= maxn_gurobi
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in EIPS v GRID)? tks infeas = $(
            round(max(0, sum(out["GRID"][:xstar][1:k])-r),digits=5)
          ), diff_14 = $(
            round(diff_14,digits=8)
          )")
        elseif argmax(diffs) == 4 && maximum(diffs) > tol && n <= maxn_gurobi
          println("(n, repi, r, k) = $((n, repi, r, k))")
          @warn("numerical error in EIPS v GURO)? tks infeas = $(
            round(max(0, sum(out["GURO"][:xstar][1:k])-r),digits=5)
          ), diff_15 = $(
            round(diff_15,digits=8)
        )")
        end


      end
      # GC.enable(true)
    end
    # GC.enable(true)
  end
end

function writeout_piv(D::Dict, n::Integer, repi::Integer, datapath::String, distribution::String)
  if isequal(D[:name], "EIPS")
    id = 1
  elseif isequal(D[:name], "ESGS")
    id = 2
  elseif isequal(D[:name], "PLCP")
    id = 3
  elseif isequal(D[:name], "GRID")
    id = 4
  elseif isequal(D[:name], "GURO")
    id = 5
  else
    throw()
  end
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_case_$repi.csv", D[:case])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_total_$repi.csv", D[:t_init] + D[:t_run] + D[:t_primal])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_run_$repi.csv", D[:t_run])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_primal_$repi.csv", D[:t_primal])
  # writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_tks_$repi.csv", D[:tks])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_t_init_$repi.csv", D[:t_init])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_obj_$repi.csv", D[:obj])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_p_$repi.csv", D[:nit_p])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_e_$repi.csv", D[:nit_e])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_init_$repi.csv", D[:nit_init])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_nit_total_$repi.csv", D[:nit_total])
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_case_$repi.csv", D[:case])
end




using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings, Parsers
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")
include(PROJPATH * "/experiment/data_process.jl")
global const DATAPATH = PROJPATH * "/time_compare"

nlevel = 10 .^(collect(1:5));
# nlevel = [10^3]
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
ri = [2; 10; 11; 13; 14; 15];
ki = [1; 6; 10];
klevel_sub = klevel[ki];
rlevel_sub = rlevel[ri];

expers = [(x, y) for x in rlevel_sub for y in klevel_sub];
repeat = 100;
rep_offset = 0;
distribution = "Uniform";
for ni in eachindex(nlevel)
  local n = nlevel[ni]
  Random.seed!(n)
  println("===================")
  println("   n=$n  ")
  println("===================")
  time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, 1_000_000, 1_000_000, expers)
  # time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, 1, 1, expers)
end

out = load_results(DATAPATH, 1, 100, 1:5, 10^6, 10^7)
process_load_result(out, DATAPATH, 1:5, 10^6, 10^7)





#################################################################################################################################################

## deal with the case for GRID where n = 10^5
nlevel = [10^1; 10^5];
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10]

ri = [2; 10; 11; 13; 14; 15]
ki = [1; 6; 10]
klevel_sub = klevel[ki]
rlevel_sub = rlevel[ri]
expers = [(x, y) for x in rlevel_sub for y in klevel_sub]
repeat = 2;
rep_offset = 0;
distribution = "Uniform";
for ni in eachindex(nlevel)
  local n = nlevel[ni]
  Random.seed!(n)
  println("===================")
  println("   n=$n  ")
  println("===================")
  time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, 100_000_000, 100_000_000, expers, 0, false, false, false, true, false)
end


#################################################################################################################################################

# deal with the case for GRID and GURO where n = 10^6, 10^7
nlevel = [10^1; 10^6; 10^7];
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10]
ri = [2; 10; 11; 13; 14; 15]
ki = [1; 6; 10]
klevel_sub = klevel[ki]
rlevel_sub = rlevel[ri]
expers = [(x, y) for x in rlevel_sub for y in klevel_sub]
repeat = 2;
rep_offset = 0;
distribution = "Uniform";
for ni in eachindex(nlevel)
  local n = nlevel[ni]
  Random.seed!(n)
  println("===================")
  println("   n=$n  ")
  println("===================")
  time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, 100_000_000, 100_000_000, expers, 0, false, false, false, true, true)
end



# #################################################################################################################################################

# # 处理当n=10^7的case 
# # cannot use garbage collector
# nlevel = [10^1; 10^7];
# rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
# klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
# # klevel_sub = [1//10^4; 2//10; 4//10; 8//10] !!!!!!!!!!!!!
# # rlevel_sub = [1//10; 5//10; 9//10; 1//1; 110//100; 150//100] !!!!!!!!!!!!!!!
# # ri = [13;14;15]
# ri = [2; 10; 11; 13; 14; 15];
# ki = [1; 6; 10];
# klevel_sub = klevel[ki];
# rlevel_sub = rlevel[ri];
# expers = [(x, y) for x in rlevel_sub for y in klevel_sub];
# repeat = 2;
# rep_offset = 0;
# distribution = "Uniform";
# for ni in eachindex(nlevel)
#   local n = nlevel[ni]
#   Random.seed!(n)
#   println("===================")
#   println("   n=$n  ")
#   println("===================")
#   time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, 1, 1000_000_000, expers, 0, false, false, false, true, false)
# end



# #################################################################################################################################################

# # 处理一些不稳定的case
# # nlevel = 10 .^(collect(1:4));
# nlevel = [10^1; 10^5];
# rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
# klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
# # klevel_sub = [1//10^4; 2//10; 4//10; 8//10] !!!!!!!!!!!!!
# # rlevel_sub = [1//10; 5//10; 9//10; 1//1; 110//100; 150//100] !!!!!!!!!!!!!!!
# # ri = [13;14;15]
# ri = [2; 10; 11; 13; 14; 15];
# ki = [1; 6; 10];
# klevel_sub = klevel[ki];
# rlevel_sub = rlevel[ri];
# expers = [(x, y) for x in rlevel_sub for y in klevel_sub];
# repeat = 100;
# rep_offset = 0;
# distribution = "Uniform";
# for ni in eachindex(nlevel)
#   local n = nlevel[ni]
#   Random.seed!(n)
#   println("===================")
#   println("   n=$n  ")
#   println("===================")
#   time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, 10^6, 1, expers, 0, false, true, true, false, true)
# end

# out = load_results(DATAPATH, 1, 100, 10^8, 10^8)

# process_load_result(out, DATAPATH, 10^8, 10^8)


#################################################################################################################################################


# global const DATAPATH = PROJPATH * "/time_compare_most"
nlevel = [10^1; 10^7];
# nlevel = [10^3]
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
repeat = 100;
rep_offset = 0;
distribution = "Uniform";
for ni in eachindex(nlevel)
  local n = nlevel[ni]
  Random.seed!(n)
  println("===================")
  println("   n=$n  ")
  println("===================")
  time_test(n, repeat, rlevel, klevel, DATAPATH, distribution, 1, 1)
end

root_dir = DATAPATH * "/Uniform";
rm(root_dir * "/10", recursive=true)

out = load_results(DATAPATH, 1, 100, 1:5)

process_load_result(out, DATAPATH, 1:5)