#=================================================
load the result
input:
  - datapath
  - rep_lo: lower bound of nrep
  - rep_up: upper bound of nrep
output:
  - out
=================================================#
function load_results(datapath::String, rep_lo::Integer, rep_hi::Integer, algos,
  maxn_grid::Integer=10^5, maxn_gurobi::Integer=10^5
)
  out = Dict()
  distributions = readdir(datapath)
  distributions = distributions[.!occursin.(".xlsx", distributions) .* .!occursin.(".DS_Store", distributions)]
  # distributions = ["Uniform"] ###############
  for di in eachindex(distributions)
    distribution = distributions[di]
    out[distribution] = Dict()
    println("============================")
    println("distribution = $(distribution)")
    expers = readdir(datapath * "/" * distribution)
    expers = expers[.!occursin.(".csv", expers) .* .!occursin.(".tex", expers) .* .!occursin.(".DS_Store", expers)]
    ns = [parse(Int64, x) for x in expers]
    exper_mask = ones(Bool, length(expers))
    for i in eachindex(expers)
      exper = expers[i]
      files = readdir(datapath * "/" * distribution * "/" * exper * "/")
      # files = files[.!occursin.(".DS_Store", files)]
      rep_mask = [x[1] for x in split.(files, ".csv")]
      rep_mask = [x[end] for x in split.(rep_mask, "_")]
      rep_mask = [(isnothing(x) ? false : (parse(Int64, x.match) >= rep_lo && parse(Int64, x.match) <= rep_hi)) for x in match.(r"(\d+)", rep_mask)]
      if sum(rep_mask) == 0
        exper_mask[i] = false
      end
    end

    for (exper,n) in zip(expers[exper_mask],ns[exper_mask])
      println("============================")
      println("exper = $exper")
      out[distribution][exper] = Dict()
      files = readdir(datapath * "/" * distribution * "/" * exper * "/")
      rep_mask = [x[1] for x in split.(files, ".csv")]
      rep_mask = [x[end] for x in split.(rep_mask, "_")]
      rep_mask = [(isnothing(x) ? true : (parse(Int64, x.match) >= rep_lo && parse(Int64, x.match) <= rep_hi)) for x in match.(r"(\d+)", rep_mask)]
      files = files[rep_mask]
      # active = files[occursin.("_active", files)]
      # bestfeas = files[occursin.("bestfeas", files)]
      case = files[occursin.("case", files)];
      nit_init = files[occursin.("_nit_init", files)];
      nit_p = files[occursin.("_nit_p", files)];
      nit_e = files[occursin.("_nit_e", files)];
      nit_total = files[occursin.("_nit_total", files)];
      # nactivecon = files[occursin.("nactivecon", files)]
      obj = files[occursin.("obj", files)];
      # ord = files[occursin.("ord", files)]
      t_init = files[occursin.("_t_init", files)]
      t_sort = files[occursin.("t_sort", files)]
      t_run = files[occursin.("_t_run", files)]
      t_primal = files[occursin.("_t_primal", files)]
      t_psort = files[occursin.("t_psort", files)];
      t_total = files[occursin.("_t_total", files)];
      outnames = ["case", "nit_init", "nit_p", "nit_e", "nit_total", "obj", "t_total", "t_run", "t_primal", "t_init", "t_sort", "t_psort"];
      outfiles = [case, nit_init, nit_p, nit_e, nit_total, obj, t_total, t_run, t_primal, t_init, t_sort, t_psort];
      klevel = files[occursin.("klevel", files)];
      rlevel = files[occursin.("rlevel", files)];
      # maxn_grid = files[occursin.("maxn_grid", files)]
      # maxn_gurobi = files[occursin.("maxn_gurobi", files)]
      nrep = maximum([parse(Int64,x.match) for x in match.(r"\d+(?=.csv)",t_total)]) - minimum([parse(Int64,x.match) for x in match.(r"\d+(?=.csv)",t_total)]) + 1
      k_list = readdlm(datapath * "/" * distribution * "/" * exper * "/" * klevel[1])
      nk = length(k_list)
      r_list = readdlm(datapath * "/" * distribution * "/" * exper * "/" * rlevel[1])
      nr = length(r_list)
      dims = (nrep,nr,nk)
      # out[exper]["rlevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in r_list])
      out[distribution][exper]["rlevel"] = vec(r_list)
      # out[distribution][exper]["klevel"] = vec([[r[1]//r[2] for r in eachcol(parse.(Int64,split(x, "//")))][1] for x in k_list])
      out[distribution][exper]["klevel"] = vec(k_list)
      for algid in algos
        if algid == 4 && maxn_grid < n
          continue
        end
        if (algid == 5) && maxn_gurobi < n
          continue
        end
        out[distribution][exper][algid] = Dict()
        for (on,of) in zip(outnames,outfiles)
          println(on)
          if on == "t_sort" || on == "t_psort"
            continue          
          else
            x = [readdlm(datapath * "/" * distribution * "/" * exper * "/" * f) for f in of[occursin.("out_$(algid)_", of)]]
            if length(vcat(x...)) != dims[1] * dims[2] *dims[3]
              special = true
              out[distribution][exper][algid][on] = y = hvncat((Int64(length(vcat(x...))/dims[2]/dims[3]), dims[2], dims[3]), true, vcat(x...)...)
            else
              special = false
              out[distribution][exper][algid][on] = y = hvncat(dims, true, vcat(x...)...) # y[i,:,:]==x[i]
            end
          end
          if special
            @simd for i in 1:Int64(length(vcat(x...))/dims[2]/dims[3])
              @assert(y[i,:,:]==x[i])
            end
          else
            @simd for i in 1:nrep
              @assert(y[i,:,:]==x[i])
            end
          end
        end
      end
      x = [readdlm(datapath * "/" * distribution * "/" * exper * "/" * f) for f in t_sort];
      out[distribution][exper]["sort_time"] = vec(vcat(x...));
      x = [readdlm(datapath * "/" * distribution * "/" * exper * "/" * f) for f in t_psort];
      out[distribution][exper]["partial_sort_time"] = hvncat((nrep, length(k_list)), true, vcat(x...)...);
      println("-> $exper loaded")
    end
  end
  return out
end


#=================================================
process the result
input:
  - out: the loaded result from the function load_results
  - datapath
  - algos: the alg data
output:
  - calculate the mean and the standard error
=================================================#
function process_load_result(out::Dict, datapath::String, algos, maxn_grid::Integer=10^5, maxn_gurobi::Integer=10^5)
  # datapath = datapath * "/" * distribution * "/"
  distributions = readdir(datapath)
  distributions = distributions[.!occursin.(".xlsx", distributions) .* .!occursin.(".DS_Store", distributions)]
  # distributions = ["Uniform"]
  for di in eachindex(distributions)
    distribution = distributions[di]
    expers = readdir(datapath * "/" * distribution)
    expers = expers[.!occursin.(".csv", expers) .* .!occursin.(".tex", expers) .* .!occursin.(".DS_Store", expers)]
    for i in eachindex(expers)
      n = expers[i]

      result = out[distribution][n];
      println(datapath * "/" * distribution * "/" * n * "/process/")
      mkpath(datapath * "/" * distribution * "/" * n * "/process/");
      mean_t_psort = mapslices(mean, result["partial_sort_time"], dims=1)
      writedlm(datapath * "/" * distribution * "/" * n * "/process/partial_sort_time_mean.csv", mean_t_psort)
      std_t_psort = mapslices(std, result["partial_sort_time"], dims=1)
      writedlm(datapath * "/" * distribution * "/" * n * "/process/partial_sort_time_std.csv", std_t_psort)

      mean_t_sort = mapslices(mean, result["sort_time"], dims=1)
      writedlm(datapath * "/" * distribution * "/" * n * "/process/sort_time_mean.csv", mean_t_sort)
      std_t_sort = mapslices(std, result["sort_time"], dims=1)
      writedlm(datapath * "/" * distribution * "/" * n * "/process/sort_time_std.csv", std_t_sort)
      for algo in algos #!!!
        if algo == 4 && maxn_grid < Parsers.parse(Int64, n)
          continue
        end
        if (algo == 5) && maxn_gurobi < Parsers.parse(Int64, n)
          continue
        end
        for key in keys(result[algo])
          mean_key = mapslices(mean, result[algo][key], dims=1)
          mean_key = dropdims(mean_key, dims=1)
          writedlm(datapath * "/" * distribution * "/" * n * "/process/$(algo)_$(key)_mean.csv", mean_key)
          if key == "t_total" || key == "t_init"
            mean_key = mapslices(mean, result[algo][key], dims=1)
            mean_key = dropdims(mean_key, dims=1)
            writedlm(datapath * "/" * distribution * "/" * n * "/process/$(algo)_$(key)_mean.csv", mean_key)

            std_t_total = mapslices(std, result[algo][key], dims=1) 
            std_t_total = dropdims(std_t_total, dims=1)
            writedlm(datapath * "/" * distribution * "/" * n * "/process/$(algo)_$(key)_std.csv", std_t_total)
          end
        end
      end
    end
  end

end



#=================================================
load the processed result used to plot
input:
  - datapath
  - algid
  - r_k_pair: the configuration of r and k we care
output:
  - n
  - data: processed data
  - out: raw data
=================================================#
function process_data_plot(datapath::String, algid::Vector{Int}, r_k_pair::Vector{Tuple{Int, Int}})
  root_dir = datapath * "/Uniform"
  files = readdir(root_dir);
  files = files[.!occursin.(".DS_Store", files)];
  n = sort([Parsers.parse(Int, sample_size) for sample_size in files])
  # if 10 in n
  #   n = n[2:end]
  # end
  files = [string(size) for size in n]
  out = Dict()
  sort_time = [];
  psort_time = [];

  file_name = ["total_mean", "total_std", "init_mean", "run_mean", "primal_mean"]
  file_name2 = ["nit_total_mean", "case_mean"]

  for folder in files
    # if folder == "10"
    #   continue
    # end
    out[folder] = Dict()
    process_path = joinpath(root_dir, folder, "process");
    sort_time_path = joinpath(process_path, "sort_time_mean.csv")
    sort_val = readdlm(sort_time_path)
    push!(sort_time, sort_val)
    psort_time_path = joinpath(process_path, "partial_sort_time_mean.csv")
    psort_val = readdlm(psort_time_path)
    push!(psort_time, psort_val)  

    for id in algid
      out[folder][id] = Dict()
      for data in file_name
        file = joinpath(process_path, string(id) * "_t_" * data * ".csv");
        out[folder][id][data] = readdlm(file)
      end
      for data in file_name2
        file = joinpath(process_path, string(id) * "_" * data * ".csv");
        out[folder][id][data] = readdlm(file)
      end

    end
  end
  
  data = Dict()
  data["sort_time"] = sort_time
  data["psort_time"] = psort_time

  for id in algid
    data[id] = Dict()
    for d in [file_name; file_name2]
      data[id][d] = Dict()
      for pair in r_k_pair
        rid, kid = pair
        val = []
        for folder in files
          # if folder == "10"
          #   continue
          # end
          push!(val, out[folder][id][d][rid, kid])
          data[id][d][pair] = val
        end
      end
    end
  end

  return n, data, out
end


# used to plot
my_color = Dict(
  "blackblue" => RGB(55/255, 103/255, 149/255), 
  "midblue" => RGB(082/255, 143/255, 173/255),
  "lightblue" => RGB(114/255, 188/255, 213/255),
  "yellow" => RGB(255/255, 208/255, 111/255),
  "red" => RGB(231/255, 098/255, 084/255),
  "orange"=> RGB(239/255, 138/255, 071/255),
  "green" => RGB(033/255, 158/255, 188/255),
  "blackyellow" => RGB(131/255, 064/255, 038/255)
);

my_color_list = [RGB(55/255, 103/255, 149/255), RGB(114/255, 188/255, 213/255), RGB(255/255, 208/255, 111/255), RGB(231/255, 098/255, 084/255)];

shape = [:circle, :diamond, :utriangle, :square, :star, :dtriangle, :pentagon, :cross, :xcross]