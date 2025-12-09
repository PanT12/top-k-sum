using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings, DataFrames
# using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics, Parsers, Plots, Measures, StatsPlots
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/src/base.jl")
include(PROJPATH * "/experiment/data_process.jl")
global const DATAPATH = PROJPATH * "/initial_point_selecting"


rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];

r_k_pair = [(2, 6), (2, 9), (10, 6), (10, 7)];

expers = [(rlevel[r], klevel[k]) for (r, k) in r_k_pair];
alg_name = ["ConInit", "AdaInit"];
algid = [1, 2];


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

function initial_point_comparing_plot(pair, data, left::Bool=false, right::Bool=false, bottom::Bool=false)
  flag = Int(mean(data[2]["case_mean"][pair][3:end]))
  sort_time = vcat(data["sort_time"]...)
  p = plot(log10.(data[2]["total_mean"][pair]), marker=:circle, label=alg_name[2], color=my_color["blackblue"], 
    linewidth=3, markersize = 6, markerstrokewidth=0, legend=true,
    legendfontsize=12, 
    legend_background_color = RGBA(1, 1, 1, 0.8)
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
    alpha = 0.3, ylims = (0, 25), legend=:bottomright, 
    legendfontsize=12, legendfontcolor=my_color["blackyellow"],
    legend_background_color = RGBA(1, 1, 1, 0.8)
  )
  Plots.annotate!(yaxis2, [(i-0.15, height+0.5, text(string(round(height, digits=1)), 10, my_color["blackblue"], :center))
    for (i, height) in enumerate(data[2]["nit_total_mean"][pair])], 
    xticks=xticks
  )
  bar!(yaxis2, (1:6).+0.15, data[1]["nit_total_mean"][pair], 
    bar_width = 0.3, color=my_color["red"], label = alg_name[1],
    alpha = 0.3, ylims = (0, 20)
  )
  Plots.annotate!(yaxis2, [(i+0.15, height+0.5, text(string(round(height, digits=1)), 10, my_color["red"], :center))
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

  title!("Flag = $(flag): \$\\tau_r\$ = $(numerator(taur)) / $(denominator(taur)), \$\\tau_k\$ = $(numerator(tauk)) / $(denominator(tauk))", 
    fontfamily="Times New Roman"
  )
  plot!(framestyle=:box)
  plot!(grid=true)

  display(plot!)
  return p
end


n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);

p = initial_point_comparing_plot(r_k_pair[3], data, true, true, true)
Plots.savefig(p, "plot/initial_point_comparing_flag2.pdf");

p = initial_point_comparing_plot(r_k_pair[4], data, true, true, true)
Plots.savefig(p, "plot/initial_point_comparing_flag3.pdf");




function distributed(data, pair)
  flag = Int(mean(data[2]["case_mean"][pair][3:end]))
  init = [data[i]["init_mean"][pair][end] for i in 2:-1:1]
  pivot = [data[i]["run_mean"][pair][end] for i in 2:-1:1]
  exact = [data[i]["primal_mean"][pair][end] for i in 2:-1:1]

  df = DataFrame(x=["AdaInit", "ConInit"], init=init, pivot=pivot, exact=exact);
  p = groupedbar(["AdaInit", "ConInit"], 
    hcat(df.exact, df.pivot, df.init), bar_position=:stack, 
    ylims = (0,0.4), bar_width=0.3,
    label=["exact" "pivot" "init"],
    color = [my_color["lightblue"] my_color["midblue"] my_color["blackblue"]],
    ytickfontsize=25, xtickfontsize=30,
    legendfontsize=30, legend=(0.5, 0.9), left_margin=2mm,
    size=(800, 700)
  )
  rid, kid = pair
  taur, tauk = rlevel[rid], klevel[kid]

  title!("Flag = $(flag): \$ \\tau_r \$ = $(numerator(taur)) / $(denominator(taur)), \$ \\tau_k \$ = $(numerator(tauk)) / $(denominator(tauk))", 
    fontfamily="Times New Roman", titlefont=30
  )
  ylabel!("Time (sec)", fontfamily="Times New Roman", yguidefontsize=30)
  yticks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
  yticklabels = [string(size) for size in yticks]
  plot!(yticks=yticks,yticklabels=yticklabels)
  plot!(framestyle=:box)
  return p
end

p = distributed(data, r_k_pair[2])
Plots.savefig(p, "plot/time_distributed_flag1.pdf");

p = distributed(data, r_k_pair[3])
Plots.savefig(p, "plot/time_distributed_flag2.pdf");

p = distributed(data, r_k_pair[4])
Plots.savefig(p, "plot/time_distributed_flag3.pdf");





# heatmap
alg_name = ["ConInit", "AdaInit"];
algid = [1, 2];
n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);

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
    colorbar=false, 
    clims=(-1, 4), 
    colormap=(cgrad(:OrRd)),
    xticks=(1:length(klevel), string.(float.(klevel))),  
    guidefont=font(15),
    xrotation=45,
    yticks=(1:length(rlevel), string.(float.(rlevel))),       
    yflip=true,
    grid=true,
    size = (700, 800), 
    left_margin = 4mm
  )
  title!("Mean Flag at \$n=10^7\$", 
    fontfamily="Times New Roman", titlefont=20
  )
  ann1 = [(j, i, isinteger(flag[i, j]) ? text(string(Int(flag[i, j])), :balck) : string(round(flag[i, j], digits=2)), 11) 
    for i in 1:length(rlevel) for j in 1:length(klevel) if flag[i, j] == -1.0 || flag[i, j] == 0.0 || flag[i, j] == 1.0
  ]
  annotate!(ann1)
  ann2 = [(j, i, isinteger(flag[i, j]) ? text(string(Int(flag[i, j])), :white) : string(round(flag[i, j], digits=2)), 11) 
    for i in 1:length(rlevel) for j in 1:length(klevel) if flag[i, j] == 2.0 || flag[i, j] == 3.0
  ]
  annotate!(ann2)
  # annotate!(ann2)
  plot!(yguidefontsize=25, ytickfontsize=14, xguidefontsize=25, xtickfontsize=14)
end

savefig(h2, "plot/initial_point_compare_heatmap_meanflag.pdf");



data1 = readdlm(joinpath(DATAPATH, "Uniform", "10000000", "process", "1_t_total_mean.csv"))
data2 = readdlm(joinpath(DATAPATH, "Uniform", "10000000", "process", "2_t_total_mean.csv"))
ratio = data1 ./ data2
h1 = begin
  heatmap(ratio, 
    xlabel="\$\\tau_k\$", 
    ylabel="\$\\tau_r\$", 
    colorbar=true, 
    clims=(0.9, 2), 
    color=cgrad(:PuBu, [0.10, 0.90]),
    xticks=(1:length(klevel), string.(float.(klevel))),  
    guidefont=font(15),
    xrotation=45,
    yticks=(1:length(rlevel), string.(float.(rlevel))),       
    yflip=true,
    grid=true,
    size = (700, 800)
  )
  title!("ConInit vs AdaInit at \$n=10^7\$", 
    fontfamily="Times New Roman", titlefont=20
  )
  ann = [(j, i, text(string(round(ratio[i, j], digits=1)), 12, :white))
    for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] >= 2
  ]
  ann2 = [(j, i, text(string(round(ratio[i, j], digits=1)), 12, :black))
    for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] < 1
  ]
  ann3 = [(j, i, text(string("+"), 10, :black))
    for i in 1:length(rlevel) for j in 1:length(klevel) if data_flag[i, j] == -1.0
  ]
  ann4 = [(j, i, text(string("+"), 10, :black))
    for i in 1:length(rlevel) for j in 1:length(klevel) if data_flag[i, j] == 1.0
  ]
  ann5 = [(j, i, text(string("+"), 10, :black))
    for i in 1:length(rlevel) for j in 1:length(klevel) if data_flag[i, j] == 0.0
  ]
  annotate!(ann)
  annotate!(ann2)
  annotate!(ann3)
  annotate!(ann4)
  annotate!(ann5)
  plot!(yguidefontsize=25, ytickfontsize=14, xguidefontsize=25, xtickfontsize=14)

end

savefig(h1, "plot/initial_point_compare_heatmap_ratio.pdf");