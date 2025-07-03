using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings, DataFrames
# using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics, Parsers, Plots, Measures, StatsPlots
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/src/base.jl")
# include(PROJPATH * "/experiment/helper_experiment.jl")
global const DATAPATH = PROJPATH * "/initial_point_selecting"


rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];

r_k_pair = [(2, 6), (2, 9), (10, 6), (10, 7)];

expers = [(rlevel[r], klevel[k]) for (r, k) in r_k_pair];
alg_name = ["without init", "with init"];
algid = [1, 2];

# function process_data_plot(datapath::String, algid::Vector{Int}, r_k_pair::Vector{Tuple{Int, Int}})
#   root_dir = datapath * "/Uniform"
#   files = readdir(root_dir);
#   files = files[.!occursin.(".DS_Store", files)];
#   n = sort([Parsers.parse(Int, sample_size) for sample_size in files])
#   # if 10 in n
#   #   n = n[2:end]
#   # end
#   files = [string(size) for size in n]
#   out = Dict()
#   sort_time = [];
#   psort_time = [];

#   file_name = ["total_mean", "total_std", "init_mean", "run_mean", "primal_mean"]
#   file_name2 = ["nit_total_mean", "case_mean"]

#   for folder in files
#     # if folder == "10"
#     #   continue
#     # end
#     out[folder] = Dict()
#     process_path = joinpath(root_dir, folder, "process");
#     sort_time_path = joinpath(process_path, "sort_time_mean.csv")
#     sort_val = readdlm(sort_time_path)
#     push!(sort_time, sort_val)
#     psort_time_path = joinpath(process_path, "partial_sort_time_mean.csv")
#     psort_val = readdlm(psort_time_path)
#     push!(psort_time, psort_val)  

#     for id in algid
#       out[folder][id] = Dict()
#       for data in file_name
#         file = joinpath(process_path, string(id) * "_t_" * data * ".csv");
#         out[folder][id][data] = readdlm(file)
#       end
#       for data in file_name2
#         file = joinpath(process_path, string(id) * "_" * data * ".csv");
#         out[folder][id][data] = readdlm(file)
#       end

#     end
#   end
  
#   data = Dict()
#   data["sort_time"] = sort_time
#   data["psort_time"] = psort_time

#   for id in algid
#     data[id] = Dict()
#     for d in [file_name; file_name2]
#       data[id][d] = Dict()
#       for pair in r_k_pair
#         rid, kid = pair
#         val = []
#         for folder in files
#           # if folder == "10"
#           #   continue
#           # end
#           push!(val, out[folder][id][d][rid, kid])
#           data[id][d][pair] = val
#         end
#       end
#     end
#   end

#   return n, data, out
# end

# my_color = Dict(
#   "blackblue" => RGB(55/255, 103/255, 149/255), 
#   "midblue" => RGB(082/255, 143/255, 173/255),
#   "lightblue" => RGB(114/255, 188/255, 213/255),
#   "yellow" => RGB(255/255, 208/255, 111/255),
#   "red" => RGB(231/255, 098/255, 084/255),
#   "orange"=> RGB(239/255, 138/255, 071/255),
#   "green" => RGB(033/255, 158/255, 188/255),
#   "blackyellow" => RGB(131/255, 064/255, 038/255)
# );

function initial_point_comparing_plot(pair, data, left::Bool=false, right::Bool=false, bottom::Bool=false)
  flag = Int(mean(data[2]["case_mean"][pair][3:end]))
  sort_time = vcat(data["sort_time"]...)
  p = plot(log10.(data[2]["total_mean"][pair]), marker=:circle, label=alg_name[2], color=my_color["blackblue"], 
    linewidth=3, markersize = 6, markerstrokewidth=0, legend=false
  )
  plot!(log10.(data[1]["total_mean"][pair]), marker=:utriangle, label=alg_name[1], color=my_color["red"], 
    linewidth=3, markersize = 6, markerstrokewidth=0
  )
  plot!(log10.(sort_time), line=:dash, label="QuickSort", color=:black, 
    linewidth=2, markersize = 6, markerstrokewidth=0
  )

  xlabel = ["10^{$(Int(i))}" for i in log10(n[1]):log10(n[end])]
  xticks = (1:length(n), xlabel)
  ylabel = ["10^{$(Int(i))}" for i in -5:-1]
  yticks = (-5:-1, ylabel)
  plot!(xticks=xticks, yticks=yticks, ytickfontsize=10, xtickfontsize=10)

  if left
    ylabel!("Time (sec)", fontfamily="Times New Roman", yguidefontsize = 15)
  end
  if bottom
    xlabel!("n", fontfamily="Times New Roman", xguidefontsize = 15)
  end

  yaxis2 = twinx();
  bar!(yaxis2, (1:6).-0.15, data[2]["nit_total_mean"][pair], 
    bar_width = 0.3, color=my_color["blackblue"], label = alg_name[2], 
    alpha = 0.3, ylims = (0, 25), legend=false
  )
  Plots.annotate!(yaxis2, [(i-0.15, height+0.5, text(string(round(height, digits=1)), 9, my_color["blackblue"], :center))
    for (i, height) in enumerate(data[2]["nit_total_mean"][pair])], 
    xticks=xticks
  )
  bar!(yaxis2, (1:6).+0.15, data[1]["nit_total_mean"][pair], 
    bar_width = 0.3, color=my_color["red"], label = alg_name[1],
    alpha = 0.3, ylims = (0, 20), legend=false
  )
  Plots.annotate!(yaxis2, [(i+0.15, height+0.5, text(string(round(height, digits=1)), 9, my_color["red"], :center))
    for (i, height) in enumerate(data[1]["nit_total_mean"][pair])], 
    xticks=xticks
  )
  plot!(yaxis2, yticks=[0, 10, 20], ytickfontsize=10, y_foreground_color_axis=my_color["blackyellow"], 
    y_guidefontcolor=my_color["blackyellow"],
    y_foreground_color_text=my_color["blackyellow"]
  )
  if right
    ylabel!(yaxis2, "Number of Iterations", fontfamily="Times New Roman", yguidefontsize=15)
  end

  rid, kid = pair
  taur, tauk = rlevel[rid], klevel[kid]

  title!("Flag = $(flag): \$ \\tau_r = $(numerator(taur)) / $(denominator(taur)), 
    \\tau_k = $(numerator(tauk)) / $(denominator(tauk)) \$", fontfamily="Times New Roman"
  )
  plot!(framestyle=:box)
  plot!(grid=true)

  display(plot!)
  return p
end



# initial_point_comparing_plot(r_k_pair[1], data, true, true, true)

n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);
# four subplot
figures = [];
for i in eachindex(r_k_pair)
  left, right, bottom = false, false, false
  pair = r_k_pair[i]
  if i == 1 || i == 3
    left = true
  end
  if i == 3 || i == 4
    bottom = true
  end
  if i == 2 || i == 4
    right = true
  end
  p = initial_point_comparing_plot(pair, data, left, right, bottom)
  push!(figures, p)
end

# p = plot(figures..., layout=(2,2), size=(1200, 600), 
#   left_margin = [5mm 0mm], right_margin = [0mm 5mm],
#   bottom_margin = [5mm 0mm]
# )

# Plots.savefig(p, "initial_point_comparing.pdf");


p1=bar([0, 0]', grid = false, showaxis=false, legend_column = -1, 
  label=["with init" "without init"], fontfamily = "Times New Roman", 
  legend=:outertop, bottom_margin = -100mm, size = (1000, 70),
  color = [my_color["blackblue"] my_color["red"]], alpha=0.3, 
  left_margin = -100mm, legendfontsize = 12
)



p0=plot([0, 0, 0]',grid=false, showaxis=false,
  label=["with init" "without init" "quicksort"], 
  color=[my_color["blackblue"] my_color["red"] :black],
  shape = [:circle :utriangle :none], line = [:solid :solid :dash],
  legend=:outertop, legend_column = -1, size = (1000, 70),
  legendfontsize = 12, markerstrokewidth=0,
  fontfamily = "Times New Roman", markersize = 6
)

lay = @layout([[a b]; [c d]; [e{0.1h} f{0.1h}]]);
p = plot(figures..., p0, p1, layout=lay, size=(1200, 650), 
  bottommargin=2mm, leftmargin=[5mm 0mm], right_margin = [0mm 5mm]
)


# p = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 400), leftmargin=4mm, bottommargin=4mm)
Plots.savefig(p, "plot/initial_point_comparing.pdf")





function distributed(data, pair, left::Bool=false)
  flag = Int(mean(data[2]["case_mean"][pair][3:end]))
  init = [data[i]["init_mean"][pair][end] for i in 2:-1:1]
  pivot = [data[i]["run_mean"][pair][end] for i in 2:-1:1]
  exact = [data[i]["primal_mean"][pair][end] for i in 2:-1:1]

  df = DataFrame(x=["with init", "without init"], init=init, pivot=pivot, exact=exact);
  p = groupedbar(["with init", "without init"], 
    hcat(df.exact, df.pivot, df.init), bar_position=:stack, 
    ylims = (0,0.4), bar_width=0.3,
    label=["exact" "pivot" "init"], legend=false,
    color = [my_color["lightblue"] my_color["midblue"] my_color["blackblue"]],
    ytickfontsize=10, xtickfontsize=10
  )
  rid, kid = pair
  taur, tauk = rlevel[rid], klevel[kid]

  title!("Flag = $(flag): \$ \\tau_r = $(numerator(taur)) / $(denominator(taur)), 
    \\tau_k = $(numerator(tauk)) / $(denominator(tauk)) \$", 
    fontfamily="Times New Roman", titlefont=11)
  if left
    ylabel!("Time (sec)", fontfamily="Times New Roman")
  end
  yticks = [0.1, 0.2, 0.3, 0.4, 0.5]
  yticklabels = [string(size) for size in yticks]
  plot!(yticks=yticks,yticklabels=yticklabels)
  plot!(framestyle=:box)
  return p
end

figures = [];
for i in eachindex(r_k_pair)
  left, right = false, false
  if i == 1
    left = true
  end
  p = distributed(data, r_k_pair[i], left)
  push!(figures, p)
end

# legend
p1=bar([0, 0, 0]', grid = false, showaxis = false, legend_column = -1, 
  label = ["init" "pivot" "exact"], fontfamily = "Times New Roman", 
  legend = :outertop, bottom_margin = -100mm, size = (1000, 70),
  color = [my_color["blackblue"] my_color["midblue"] my_color["lightblue"]], 
  left_margin = -100mm, legendfontsize = 10
)

lay = @layout([[a b c d]; e{0.1h}]);
p = plot(figures...,p1, layout = lay, size = (1000,300),
  left_margin = [5mm 0mm 0mm 0mm 2mm], top_margin = [3mm 0mm],
  bottom_margin = [-1mm 0mm]
)

savefig(p, "plot/time_distributed.pdf")





# heatmap
alg_name = ["without init", "with init"];
algid = [1, 2];
n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);


data1 = readdlm(joinpath(DATAPATH, "Uniform", "10000000", "process", "1_t_total_mean.csv"))
data2 = readdlm(joinpath(DATAPATH, "Uniform", "10000000", "process", "2_t_total_mean.csv"))
ratio = data1 ./ data2
# diff = out[string(n[end])][1]["total_mean"] - out[string(n[end])][2]["total_mean"];
ratio = out[string(n[end])][1]["total_mean"] ./ out[string(n[end])][2]["total_mean"];
h1 = begin
  heatmap(ratio, 
    xlabel="\$\\tau_k\$", 
    ylabel="\$\\tau_r\$", 
    colorbar=true, 
    clims=(1, 3), 
    color=cgrad(:PuBu, [0.010, 0.990]),
    xticks=(1:length(klevel), string.(float.(klevel))),  
    guidefont=font(15),
    xrotation=45,
    yticks=(1:length(rlevel), string.(float.(rlevel))),       
    yflip=true,
    grid=true,
  )
  title!("without init vs with init at \$n=10^7\$", 
    fontfamily="Times New Roman", titlefont=18
  )
  ann = [(j, i, text(string(round(ratio[i, j], digits=1)), 10, :white))
    for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] >= 3
  ]
  ann2 = [(j, i, text(string(round(ratio[i, j], digits=1)), 10, :black))
    for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] < 1
  ]
  annotate!(ann)
  annotate!(ann2)
  plot!(yguidefontsize=20, ytickfontsize=14, xguidefontsize=20, xtickfontsize=14)
end



flag = out[string(n[end])][2]["case_mean"];
data_flag = similar(flag);
for i in 1:size(flag, 1)
  for j in 1:size(flag, 2)
    a = ceil(flag[i, j])
    if abs(a - flag[i, j]) < 1e-10
      data_flag[i, j] = Int(a)
    else
      data_flag[i, j] = round(flag[i, j], digits=2)
    end
  end
end
h2 = begin
  heatmap(data_flag, 
    xlabel="\$\\tau_k\$", 
    ylabel="\$\\tau_r\$", 
    colorbar=true, 
    # clims=(minimum(ratio), 8), 
    colormap=reverse(cgrad(:OrRd)),
    xticks=(1:length(klevel), string.(float.(klevel))),  
    guidefont=font(15),
    xrotation=45,
    yticks=(1:length(rlevel), string.(float.(rlevel))),       
    yflip=true,
    grid=true,
  )
  title!("Mean Flag at \$n=10^7\$", 
    fontfamily="Times New Roman", titlefont=18
  )
  ann1 = [(j, i, isinteger(flag[i, j]) ? text(string(Int(flag[i, j])), :white) : string(round(flag[i, j], digits=2)), 11) 
    for i in 1:length(rlevel) for j in 1:length(klevel) if flag[i, j] == -1.0
  ]
  annotate!(ann1)
  ann2 = [(j, i, isinteger(flag[i, j]) ? string(Int(flag[i, j])) : string(round(flag[i, j], digits=2)), 11) 
    for i in 1:length(rlevel) for j in 1:length(klevel) if flag[i, j] != -1.0
  ]
  annotate!(ann2)
  # annotate!(ann2)
  plot!(yguidefontsize=20, ytickfontsize=14, xguidefontsize=20, xtickfontsize=14)
end


h = plot(h1, h2, lay = (1,2), size = (1600, 800), 
  left_margin = [7mm 4mm], bottom_margin = [10mm 0mm], 
  top_margin = 3mm
)

savefig(h, "plot/initial_point_compare_heatmap.pdf")