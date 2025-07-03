using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings, Parsers, StatsPlots, Printf
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics, Measures
using SparseMatricesCSR, SparseArrays, Statistics
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")
include(PROJPATH * "/src/base.jl")
global const DATAPATH = PROJPATH * "/time_compare"

# nlevel = [collect(1:10).*10^5; [1;5;10;50;100;500;1000].*10^4]
# nlevel = collect(1:10).*10^5
nlevel = [10^4; 5*10^4; 10^5; 5*10^5; 10^6; 5*10^6; 10^7];
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
ri = [2; 10; 11; 13; 14; 15];
ki = [1; 6; 10];
klevel_sub = klevel[ki];
rlevel_sub = rlevel[ri];

expers = [(x, y) for x in rlevel_sub for y in klevel_sub];
r_k_pair = [(x, y) for x in ri for y in ki];

alg_name = ["EIPS", "ESGS", "PLCP", "GRID", "GURO"];
algid = [1,2,3,4,5];

# function process_data_plot(datapath::String)
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

#   file_name = ["total_mean", "total_std", "init_mean", "init_std", "run_mean", "primal_mean"]
#   file_name2 = ["nit_mean", "case_mean"]

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
#         # if rlevel[rid] >= 1
#         #   continue
#         # end
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
#   return n, out, data
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



function time_plot_active(data, item, pair, bottom::Bool=false, left::Bool=false)

  rid, kid = pair
  p = plot(n, data[1][item][pair], marker=:circle, label=alg_name[1], color=my_color["blackblue"], 
    yscale=:log10, xscale=:log10, linewidth=3, markersize = 6, markerstrokewidth=0, legend=false
  )
  plot!(n, data[2][item][pair], marker=:square, label=alg_name[2], color=my_color["red"],linewidth=3, markersize = 6, markerstrokewidth=0)
  plot!(n, data[3][item][pair], marker=:dtriangle, label=alg_name[3], color=my_color["lightblue"],linewidth=3, markersize = 6, markerstrokewidth=0)
  plot!(n, data[4][item][pair], marker=:utriangle, label=alg_name[4], color=my_color["yellow"],linewidth=3, markersize = 6, markerstrokewidth=0)
  plot!(n, data[5][item][pair], marker=:pentagon, label=alg_name[5], color=my_color["green"],linewidth=3, markersize = 6, markerstrokewidth=0)
  if left
    ylabel!("time (sec)", fontfamily="Times New Roman")
  end
  xticks = n
  xticklabels = [string(size) for size in n]
  yticks = [10^(-7), 10^(-5), 10^(-3), 10^(-1), 10^(1), 10^(3)]
  yticklabels = [string(size) for size in yticks]
  plot!(xticks=xticks, xticklabels=xticklabels, yticks=yticks,yticklabels=yticklabels)
  plot!(yguidefontsize=17, ytickfontsize=12, xguidefontsize=16, xtickfontsize=12)


  plot!(grid=true)
  if bottom
    xlabel!("n", fontfamily="Times New Roman")
  end
  rid, kid = pair
  taur, tauk = rlevel[rid], klevel[kid]

  title!("Experiment: \$ \\tau_r = $(numerator(taur)) / $(denominator(taur)), \\tau_k = $(numerator(tauk)) / $(denominator(tauk)) \$", 
    fontfamily="Times New Roman", titlefont=16
  )

  plot!(framestyle=:box)
  return p
end

function time_plot_nonactive(data, item, ri, kid, left::Bool=false)

  p = plot(n, data[2][item][(ri[4], kid)], marker=:utriangle, 
    label=alg_name[2], yscale=:log10, xscale=:log10,
    color=my_color["yellow"],linewidth=3, markersize = 6, markerstrokewidth=0
  )
  tau_r = rlevel[ri[4]];
  plot!(n, data[1][item][(ri[4], kid)], marker=:circle, 
    label="$(alg_name[1]): \$ \\tau_r = $(numerator(tau_r)) / $(denominator(tau_r)) \$", 
    color=my_color["blackblue"], 
    linewidth=3, markersize = 6, markerstrokewidth=0
  )
  tau_r = rlevel[ri[5]];
  plot!(n, data[1][item][(ri[5], kid)], marker=:square, 
    label="$(alg_name[1]): \$ \\tau_r = $(numerator(tau_r)) / $(denominator(tau_r)) \$", 
    color=my_color["green"], 
    linewidth=3, markersize = 6, markerstrokewidth=0
  )
  tau_r = rlevel[ri[6]];
  plot!(n, data[1][item][(ri[6], kid)], marker=:dtriangle, 
    label="$(alg_name[1]): \$ \\tau_r = $(numerator(tau_r)) / $(denominator(tau_r)) \$", 
    color=my_color["red"], 
    linewidth=3, markersize = 6, markerstrokewidth=0
  )
  if left
    ylabel!("time (sec)", fontfamily="Times New Roman")
  end
  xticks = n
  xticklabels = [string(size) for size in n]
  yticks = [10^(-6), 10^(-5), 10^(-4), 10^(-3), 10^(-2)]
  yticklabels = [string(size) for size in yticks]
  plot!(xticks=xticks, xticklabels=xticklabels, yticks=yticks,yticklabels=yticklabels)

  plot!(grid=true)
  xlabel!("n", fontfamily="Times New Roman")
  tauk = klevel[kid]
  title!("Experiment: \$ \\tau_k = $(numerator(tauk)) / $(denominator(tauk)) \$", 
    fontfamily="Times New Roman", titlefont = 17
  )

  plot!(framestyle=:box)
  plot!(legend=:topleft, legendfont="Times New Roman", legendfontsize = 14)
  plot!(yguidefontsize=17, ytickfontsize=12, xguidefontsize=16, xtickfontsize=12)
  return p
end



n, out, data = process_data_plot(DATAPATH, algid, r_k_pair);

figures_nonactive = [];
for d in eachindex(ki)
  left, bottom = false, false
  if d == 1
    left = true
  elseif d == 2
    bottom = true
  end
  kid = ki[d]
  p = time_plot_nonactive(data, "init_mean", ri, kid, left);
  push!(figures_nonactive, p)
end

p = plot(figures_nonactive..., layout = (1,3), size = (1500, 400),
  top_margin = [4mm 0mm], left_margin = [9mm 0mm 0mm],
  bottom_margin = [10mm 0mm]
)

savefig(p, "plot/time_compare_nonactive.pdf")




r_k_pair_active = [(rid, kid) for (rid, kid) in r_k_pair if rlevel[rid] < 1];

left_show_label = [1; 4; 7];
bottom_show_label = [7; 8; 9];
figures_active = [];
for i in eachindex(r_k_pair_active)
  bottom, left = false, false
  pair = r_k_pair_active[i]
  println(pair)
  if i in left_show_label
    left = true
  end
  if i in bottom_show_label
    bottom = true
  end
  p = time_plot_active(data, "total_mean", pair, bottom, left)
  push!(figures_active, p)
end
lay = @layout([[a b c]; [d e f]; [g h i]; j{0.1h}])

plot(figures_active..., layout = (3,3), size = (1500,1000),
  left_margin = [7mm 1mm 1mm], margin=-1mm
)

p0=plot([0, 0, 0, 0, 0]',grid=false, marker=:circle, showaxis=false,
  label=[alg_name[1] alg_name[2] alg_name[3] alg_name[4] alg_name[5]], 
  shape=[:circle :square :dtriangle :utriangle :pentagon], 
  color=[my_color["blackblue"] my_color["red"] my_color["lightblue"] my_color["yellow"] my_color["green"]],
  legend=:outertop, legend_column = -1, size = (1000, 70), markerstrokewidth=0,
  legendfont=13, fontfamily="Times New Roman"
)
lay = @layout([[a b c]; [d e f]; [g h i]; j{0.1h}])
p = plot(figures_active...,p0, layout = lay, size = (1500,1000),
  left_margin = [7mm 1mm 1mm], bottom_margin = [1mm 2mm 3mm]
)

savefig(p, "plot/time_compare_active.pdf")







# heatmap
alg_name = ["EIPS", "ESGS", "PLCP"];
algid = [1,2,3];
n, out, data = process_data_plot(DATAPATH, algid, r_k_pair);

function time_compare_heatmap(algid)
  c = cgrad(:PuBu, [0.010, 0.990])
  # diff = out[string(n[end])][algid]["total_mean"] - out[string(n[end])][1]["total_mean"];
  ratio = out[string(n[end])][algid]["total_mean"] ./ out[string(n[end])][1]["total_mean"];
  h = begin
    heatmap(ratio, 
      xlabel="\$\\tau_k\$", 
      ylabel="\$\\tau_r\$", 
      colorbar=true, 
      clims=(1, 10), 
      # colormap=:PuBu,
      color=c,
      xticks=(1:length(klevel), string.(float.(klevel))),  
      guidefont=font(15),
      xrotation=45,
      yticks=(1:length(rlevel), string.(float.(rlevel))),       
      yflip=true,
      grid=true,
    )
    title!("$(alg_name[algid]) vs EIPS at \$n=10^7\$", 
      fontfamily="Times New Roman", titlefont=18
    )
    ann = [(j, i, text(string(round(ratio[i, j], digits=1)), 10, :white))
      for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] >= 10
    ]
    ann2 = [(j, i, text(string(round(ratio[i, j], digits=1)), 10, :black))
      for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] <= 1
    ]
    annotate!(ann)
    annotate!(ann2)
    plot!(yguidefontsize=20, ytickfontsize=14, xguidefontsize=20, xtickfontsize=14)
  end
end

figures = [];
for i in 2:3
  h = time_compare_heatmap(i)
  push!(figures, h)
end
h = plot(figures..., lay = (1,2), size = (1600, 800), 
  left_margin = [7mm 2mm], bottom_margin = [10mm 0mm], 
  top_margin = 3mm
)

savefig(h, "plot/time_compare_heatmap.pdf")




# latex table
algid = 3;
pair = (11, 10);
values = data[algid]["total_mean"][pair];
values2 = data[algid]["total_std"][pair];


formatted_values = map(value -> @sprintf("%.2e", value), values);
formatted_values = map(value -> replace(value, "e+0" => "e+", "e-0" => "e-"), formatted_values);

formatted_values2 = map(value -> @sprintf("%.2e", value), values2);
formatted_values2 = map(value -> replace(value, "e+0" => "e+", "e-0" => "e-"), formatted_values2);

for i in 3:length(formatted_values)
  print(formatted_values[i], " (", formatted_values2[i], ") & ")
end
println()