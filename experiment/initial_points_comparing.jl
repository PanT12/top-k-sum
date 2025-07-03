function initialize_output(n, nr, nk)
  alg_names = ["without_init", "with_init"]
  out = Dict()
  out["without_init"] = Dict()
  out["with_init"] = Dict()
  out["ALL"] = Dict()
  
  # t_init
  out["without_init"][:t_init] = fill(+Inf, nr, nk)
  out["with_init"][:t_init] = fill(+Inf, nr, nk)
  
  # t_run
  out["without_init"][:t_run] = fill(+Inf, nr, nk)
  out["with_init"][:t_run] = fill(+Inf, nr, nk)
  
  # t_primal: primal recovery
  out["without_init"][:t_primal] = fill(+Inf, nr, nk)
  out["with_init"][:t_primal] = fill(+Inf, nr, nk)
  
  # obj
  out["without_init"][:obj] = fill(+Inf, nr, nk) # obj value 1
  out["with_init"][:obj] = fill(+Inf, nr, nk) # obj value 1
  
  # nit
  out["without_init"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 1
  out["with_init"][:nit_init] = zeros(Int64, nr, nk) # number of init iterations 2
  out["without_init"][:nit_p] = zeros(Int64, nr, nk) # number of iterations 1
  out["with_init"][:nit_p] = zeros(Int64, nr, nk) # number of iterations 2
  out["without_init"][:nit_e] = zeros(Int64, nr, nk) # number of iterations 1
  out["with_init"][:nit_e] = zeros(Int64, nr, nk) # number of iterations 2
  out["without_init"][:nit_total] = zeros(Int64, nr, nk) # total number of iterations 1
  out["with_init"][:nit_total] = zeros(Int64, nr, nk) # total number of iterations 2


  
  # case
  out["without_init"][:case] = zeros(Int64, nr, nk) # case id: 1 = strict, 2 = strict but lcp ez; 0 = inactive; -1 = k==1; -2 = k==n
  out["with_init"][:case] = zeros(Int64, nr, nk) # case id: 0: inactive, 1: active
  
  # bestfeas
  out["without_init"][:bestfeas] = zeros(Bool, nr, nk) # does 1 satisfy Mk(xbar) <= r and have lowest objval
  out["with_init"][:bestfeas] = zeros(Bool, nr, nk) # does 1 satisfy Mk(xbar) <= r and have lowest objval

  # other
  out["ALL"][:obj_1234] = Array{Vector{Int64}}(undef, nr, nk) # best objective value
  out["ALL"][:active] = zeros(Bool, nr, nk) # is Mk(x0) > r?
  out["ALL"][:x0sort] = zeros(Float64, n) # x0sort
  out["ALL"][:x0] = zeros(Float64, n) # x0
  out["ALL"][:t_psort] = fill(+Inf, nk) # partial sort time

  
  # solution
  out["without_init"][:xstar] = zeros(Float64, n)
  out["with_init"][:xstar] = zeros(Float64, n)
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

  x0prepop = false # prepopulate xbarsort .= x0sort
  hist = false # don't record history
  verb = false # no console output


  # test output
  println(datapath * "/$(distribution)/$n/")
  mkpath(datapath * "/$(distribution)/$n/")
  writedlm(datapath*"/$(distribution)/$(n)/rlevel.csv", rlevel)
  writedlm(datapath*"/$(distribution)/$(n)/klevel.csv", klevel)
  

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
      
      # if !isnothing(expers)
      #   GC.enable(false)
      # end

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
  
        # ESGS
        xstar_esgs = similar(x0sort);
        @views sp, k0k1, nit, t, case, lam = project_topksum_esgs_experiment!(
          xstar_esgs, x0sort, x0, r, k, sum(x0sort[1:k]) > r, x0prepop, verb, hist
        )
        # time_sort = out["ALL"][:active][ri,ki] ? out["ALL"][:t_sort] + out["ALL"][:t_psort][ki] : out["ALL"][:t_psort][ki];
        # D = out["ESGS"]
        # D[:xstar] .= xstar
        # D[:name] = "ESGS"
        # D[:k0][ri,ki] = k0k1[1]
        # D[:k1][ri,ki] = k0k1[2]
        # D[:t_init][ri,ki] = t[1] + time_sort;
        # D[:t_run][ri,ki] = t[2]
        # D[:t_primal][ri,ki] = t[3] #+ time_convert
        # D[:tks][ri,ki] = max(sum(sort(xstar, rev=true)[1:k]) - r, 0.0) #! is unsorted
        # # D[:ord][ri,ki] = maximum(max.(diff(D[:xstar]), 0.0))
        # # D[:infeas][ri,ki] = max(D[:tks][ri,ki], D[:ord][ri,ki])
        # D[:obj][ri,ki] = 0.5 * norm(xstar .- x0, 2)^2
        # D[:nit][ri,ki] = nit
        # D[:case][ri,ki] = case
        # writeout_piv(D, n, repi, datapath, distribution)

        # EIPS # without_init
        xstar = similar(x0);
        debug = false;
        @views flag, Ite, t = project_topksum_EIPS_experiment!(
          xstar, x0, r, k, true, debug, hist
        )
        D = out["without_init"]
        D[:case][ri,ki] = flag
        D[:xstar] .= xstar
        D[:name] = "without_init"
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


        # EIPS # with_init
        xstar = similar(x0);
        debug = false;
        @views flag, Ite, t = project_topksum_EIPS_experiment!(
          xstar, x0, r, k, false, debug, hist
        )
        D = out["with_init"]
        D[:case][ri,ki] = flag
        D[:xstar] .= xstar
        D[:name] = "with_init"
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


        # check
        tol = 1e-8;
        diff_12 = norm(out["without_init"][:xstar] .- xstar_esgs, Inf)
        diff_13 = norm(out["with_init"][:xstar] .- xstar_esgs, Inf)
        diffs = [diff_12, diff_13];
        if maximum(diffs) > tol
          if argmax(diffs) == 1
          # if diff_12 > tol
            println("(n, repi, rlevel, klevel) = $((n, repi, ri, ki))")
            @warn("numerical error in ESGS v without_init)? tks infeas = $(
              round(max(0, sum(sort(out["without_init"][:xstar], rev=true)[1:k])-r),digits=5)
            ), diff_12 = $(
              round(diff_12,digits=8)
            )")
          else
            println("(n, repi, rlevel, klevel) = $((n, repi, ri, ki))")
            @warn("numerical error in ESGS v with_init)? tks infeas = $(
              round(max(0, sum(sort(out["with_init"][:xstar], rev=true)[1:k])-r),digits=5)
            ), diff_13 = $(
              round(diff_13,digits=8)
            )")            
          end
        end

      end
      # if !isnothing(expers)
      #   GC.enable(true)
      # end

    end
  end
end

function writeout_piv(D::Dict, n::Integer, repi::Integer, datapath::String, distribution::String)
  if isequal(D[:name], "without_init")
    id = 1
  elseif isequal(D[:name], "with_init")
    id = 2
  else
    throw()
  end
  writedlm(datapath*"/$(distribution)/$(n)/out_$(id)_flag_$repi.csv", D[:case])
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

# load the result and load the result
# input: 
#   - datapath
#   - rep_lo: lower bound of nrep
#   - rep_up: upper bound of nrep
#   - distribution
# function load_results(datapath::String, rep_lo::Integer, rep_hi::Integer, 
#   maxn_grid::Integer=10^5, maxn_gurobi::Integer=10^5
# )
#   out = Dict()
#   distributions = readdir(datapath)
#   distributions = distributions[.!occursin.(".xlsx", distributions) .* .!occursin.(".DS_Store", distributions)]
#   # distributions = ["Uniform"] ###############
#   for di in eachindex(distributions)
#     distribution = distributions[di]
#     out[distribution] = Dict()
#     println("============================")
#     println("distribution = $(distribution)")
#     expers = readdir(datapath * "/" * distribution)
#     expers = expers[.!occursin.(".csv", expers) .* .!occursin.(".tex", expers) .* .!occursin.(".DS_Store", expers)]
#     # ns = [parse(Int64, match(r"(\d+)",x).match) for x in expers]
#     # expers = ["100000"] ####################
#     ns = [parse(Int64, x) for x in expers]
#     exper_mask = ones(Bool, length(expers))
#     for i in eachindex(expers)
#       exper = expers[i]
#       files = readdir(datapath * "/" * distribution * "/" * exper * "/")
#       # files = files[.!occursin.(".DS_Store", files)]
#       rep_mask = [x[1] for x in split.(files, ".csv")]
#       rep_mask = [x[end] for x in split.(rep_mask, "_")]
#       rep_mask = [(isnothing(x) ? false : (parse(Int64, x.match) >= rep_lo && parse(Int64, x.match) <= rep_hi)) for x in match.(r"(\d+)", rep_mask)]
#       if sum(rep_mask) == 0
#         exper_mask[i] = false
#       end
#     end

#     for (exper,n) in zip(expers[exper_mask],ns[exper_mask])
#       println("============================")
#       println("exper = $exper")
#       out[distribution][exper] = Dict()
#       files = readdir(datapath * "/" * distribution * "/" * exper * "/")
#       rep_mask = [x[1] for x in split.(files, ".csv")]
#       rep_mask = [x[end] for x in split.(rep_mask, "_")]
#       rep_mask = [(isnothing(x) ? true : (parse(Int64, x.match) >= rep_lo && parse(Int64, x.match) <= rep_hi)) for x in match.(r"(\d+)", rep_mask)]
#       files = files[rep_mask]
#       # active = files[occursin.("_active", files)]
#       # bestfeas = files[occursin.("bestfeas", files)]
#       case = files[occursin.("case", files)];
#       nit_init = files[occursin.("_nit_init", files)];
#       nit_p = files[occursin.("_nit_p", files)];
#       nit_e = files[occursin.("_nit_e", files)];
#       nit_total = files[occursin.("_nit_total", files)];
#       # nactivecon = files[occursin.("nactivecon", files)]
#       obj = files[occursin.("obj", files)];
#       # ord = files[occursin.("ord", files)]
#       t_init = files[occursin.("_t_init", files)]
#       t_sort = files[occursin.("t_sort", files)]
#       t_run = files[occursin.("t_run", files)]
#       t_primal = files[occursin.("t_primal", files)]
#       # t_psort = files[occursin.("t_psort", files)];
#       t_total = files[occursin.("_t_total", files)];
#       outnames = ["case", "nit_init", "nit_p", "nit_e", "nit_total", "obj", "t_total", "t_run", "t_primal", "t_init", "t_sort"];
#       # outfiles = [active, bestfeas, case, infeas, k0, k1, tks, nit, nactivecon, obj, ord, t_init, t_run, t_primal, t_psort]
#       outfiles = [case, nit_init, nit_p, nit_e, nit_total, obj, t_total, t_run, t_primal, t_init, t_sort];
#       klevel = files[occursin.("klevel", files)];
#       rlevel = files[occursin.("rlevel", files)];
#       # maxn_grid = files[occursin.("maxn_grid", files)]
#       # maxn_gurobi = files[occursin.("maxn_gurobi", files)]
#       nrep = maximum([parse(Int64,x.match) for x in match.(r"\d+(?=.csv)",t_total)]) - minimum([parse(Int64,x.match) for x in match.(r"\d+(?=.csv)",t_total)]) + 1
#       k_list = readdlm(datapath * "/" * distribution * "/" * exper * "/" * klevel[1])
#       nk = length(k_list)
#       r_list = readdlm(datapath * "/" * distribution * "/" * exper * "/" * rlevel[1])
#       nr = length(r_list)
#       dims = (nrep,nr,nk)
#       # out[exper]["rlevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in r_list])
#       out[distribution][exper]["rlevel"] = vec(r_list)
#       # out[distribution][exper]["klevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in k_list])
#       out[distribution][exper]["klevel"] = vec(k_list)
#       for algid in 1:2
#         out[distribution][exper][algid] = Dict()
#         for (on,of) in zip(outnames,outfiles)
#           println(on)
#           if on == "t_sort" || on == "t_psort"
#             continue          
#           else
#             x = [readdlm(datapath * "/" * distribution * "/" * exper * "/" * f) for f in of[occursin.("out_$(algid)_", of)]]
#             out[distribution][exper][algid][on] = y = hvncat(dims, true, vcat(x...)...) # y[i,:,:]==x[i]
#           end
#           @simd for i in 1:nrep
#             @assert(y[i,:,:]==x[i])
#           end
#         end
#       end
#       x = [readdlm(datapath * "/" * distribution * "/" * exper * "/" * f) for f in t_sort];
#       out[distribution][exper]["sort_time"] = vec(vcat(x...));
#       # x = [readdlm(datapath * "/" * distribution * "/" * exper * "/" * f) for f in t_psort];
#       # out[distribution][exper]["partial_sort_time"] = hvncat((nrep, length(k_list)), true, vcat(x...)...);
#       println("-> $exper loaded")
#     end
#   end
#   return out
# end

# function process_load_result(out::Dict, datapath::String)
#   # datapath = datapath * "/" * distribution * "/"
#   distributions = readdir(datapath)
#   distributions = distributions[.!occursin.(".xlsx", distributions) .* .!occursin.(".DS_Store", distributions)]
#   # distributions = ["Uniform"]
#   for di in eachindex(distributions)
#     distribution = distributions[di]
#     expers = readdir(datapath * "/" * distribution)
#     expers = expers[.!occursin.(".csv", expers) .* .!occursin.(".tex", expers) .* .!occursin.(".DS_Store", expers)]
#     for i in eachindex(expers)
#       n = expers[i]

#       result = out[distribution][n];
#       println(datapath * "/" * distribution * "/" * n * "/process/")
#       mkpath(datapath * "/" * distribution * "/" * n * "/process/");
#       # mean_t_psort = mapslices(mean, result["partial_sort_time"], dims=1)
#       # writedlm(datapath * "/" * distribution * "/" * n * "/process/partial_sort_time_mean.csv", mean_t_psort)
#       mean_t_sort = mapslices(mean, result["sort_time"], dims=1)
#       writedlm(datapath * "/" * distribution * "/" * n * "/process/sort_time_mean.csv", mean_t_sort)
#       std_t_sort = mapslices(std, result["sort_time"], dims=1)
#       writedlm(datapath * "/" * distribution * "/" * n * "/process/sort_time_std.csv", std_t_sort)
#       for algo = 1:2 #!!!
#         for key in keys(result[algo])
#           mean_key = mapslices(mean, result[algo][key], dims=1)
#           mean_key = dropdims(mean_key, dims=1)
#           writedlm(datapath * "/" * distribution * "/" * n * "/process/$(algo)_$(key)_mean.csv", mean_key)
#           if key == "t_total" || key == "t_init" || key == "t_run" || key == "t_primal"
#             mean_key = mapslices(mean, result[algo][key], dims=1)
#             mean_key = dropdims(mean_key, dims=1)
#             writedlm(datapath * "/" * distribution * "/" * n * "/process/$(algo)_$(key)_mean.csv", mean_key)

#             std_t_total = mapslices(std, result[algo][key], dims=1) 
#             std_t_total = dropdims(std_t_total, dims=1)
#             writedlm(datapath * "/" * distribution * "/" * n * "/process/$(algo)_$(key)_std.csv", std_t_total)
#           end
#         end
#       end

#       # for key in keys(result[1])
#       #   if key == "t_total" || key == "obj"
#       #     a = round.(readdlm(datapath * "/" * distribution * "/" * n * "/process/1_$(key)_mean.csv"), digits = 3)
#       #     b = round.(readdlm(datapath * "/" * distribution * "/" * n * "/process/2_$(key)_mean.csv"), digits = 3)
#       #     c = [string(a[i,j], "(", b[i,j], ")") for i in 1:size(a, 1), j in 1:size(a, 2)]
#       #     writedlm(datapath * "/" * distribution * "/" * n * "/process/$(key)_compare.csv", c)
#       #   end
#       # end
#       # p1 = heatmap(algo_CIPS_vs_ESGS, title="CIPS v.s. ESGS$(n)$(distribution)", xlabel="k level", ylabel="r level", colorbar=true)
#       # p2 = heatmap(algo_BIPS_vs_ESGS, title="BIPS v.s. ESGS$(n)$(distribution)", xlabel="k level", ylabel="r level", colorbar=true)
#       # p3 = heatmap(algo_CIPS_vs_BIPS, title="CIPS v.s. BIPS$(n)$(distribution)", xlabel="k level", ylabel="r level", colorbar=true)
#       # savefig(p1, datapath * "/" * distribution * "/" * n * "/process/CIPS v.s. ESGS.pdf")  # 将 "heatmap.png" 替换为你想要的文件名
#       # savefig(p2, datapath * "/" * distribution * "/" * n * "/process/BIPS v.s. ESGS.pdf")
#       # savefig(p3, datapath * "/" * distribution * "/" * n * "/process/CIPS v.s. BIPS.pdf")
#     end
#   end

# end


using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics, JuMP, Gurobi
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")
# include(PROJPATH * "/experiment/helper_experiment.jl")
global const DATAPATH = PROJPATH * "/initial_point_selecting"
# global const DATAPATH = PROJPATH * "/most_case"
# global const DATAPATH = PROJPATH * "/initial_point_selecting_7most"
# global const DATAPATH = PROJPATH * "/debug"




# debug size
# nlevel = [10^1, 10^2, 10^3, 10^4]

nlevel = 10 .^(collect(1:3));
# nlevel = [1;4;10;40;100;400].*10^4
# nlevel = [10^5]
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];

r_k_pair = [(2, 6), (2, 9), (10, 6), (10, 7)];

expers = [(rlevel[r], klevel[k]) for (r, k) in r_k_pair];

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
out = load_results(DATAPATH, 1, 100, 1:2)

process_load_result(out, DATAPATH, 1:2)



####################################################################################
# most cases

nlevel = [10^1; 10^3];
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]];
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
  time_test(n, repeat, rlevel, klevel, DATAPATH, distribution)
end

root_dir = DATAPATH * "/Uniform";
rm(root_dir * "/10", recursive=true)

out = load_results(DATAPATH, 1, 100, 1:2)

process_load_result(out, DATAPATH, 1:2)