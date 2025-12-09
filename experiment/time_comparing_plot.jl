using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings, Parsers, StatsPlots, Printf
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics, Measures
using SparseMatricesCSR, SparseArrays, Statistics
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")
include(PROJPATH * "/src/base.jl")
include(PROJPATH * "/experiment/data_process.jl")
global const DATAPATH = PROJPATH * "/time_compare"

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
    ylabel!("Time (sec)", fontfamily="Times New Roman", yguidefontsize=15)
  end
  xticks = n
  xticklabels = [string(size) for size in n]
  yticks = [10^(-7), 10^(-5), 10^(-3), 10^(-1), 10^(1), 10^(3)]
  yticklabels = [string(size) for size in yticks]
  plot!(xticks=xticks, xticklabels=xticklabels, yticks=yticks,yticklabels=yticklabels)
  plot!(yguidefontsize=17, ytickfontsize=13, xguidefontsize=16, xtickfontsize=13)


  plot!(grid=true)
  if bottom
    xlabel!("n", fontfamily="Times New Roman", xguidefontsize=15)
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
    ylabel!("Time (sec)", fontfamily="Times New Roman", yguidefontsize=15)
  end
  xticks = n
  xticklabels = [string(size) for size in n]
  yticks = [10^(-6), 10^(-5), 10^(-4), 10^(-3), 10^(-2)]
  yticklabels = [string(size) for size in yticks]
  plot!(xticks=xticks, xticklabels=xticklabels, yticks=yticks,yticklabels=yticklabels)

  plot!(grid=true)
  xlabel!("n", fontfamily="Times New Roman", xguidefontsize=15)
  tauk = klevel[kid]
  title!("Experiment: \$ \\tau_k = $(numerator(tauk)) / $(denominator(tauk)) \$", 
    fontfamily="Times New Roman", titlefont = 17
  )

  plot!(framestyle=:box)
  plot!(legend=:topleft, legendfont="Times New Roman", legendfontsize = 14)
  plot!(yguidefontsize=17, ytickfontsize=13, xguidefontsize=16, xtickfontsize=13)
  return p
end



n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);

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
  left_margin = [7mm 1mm 1mm], margin=-1mm, 
  top_margin=1mm, bottom_margin=5mm
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
  left_margin = [7mm 1mm 1mm], bottom_margin = [1mm 2mm 5mm]
)

savefig(p, "plot/time_compare_active.pdf")







# heatmap
alg_name = ["EIPS", "ESGS", "PLCP"];
algid = [1,2,3];
n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);

function time_compare_heatmap(algid)
  c = cgrad(:PuBu, [0.010, 0.990])
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
      size=(700,800)
    )
    title!("$(alg_name[algid]) vs EIPS at \$n=10^7\$", 
      fontfamily="Times New Roman", titlefont=20
    )
    ann = [(j, i, text(string(round(ratio[i, j], digits=1)), 10, :white))
      for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] >= 10
    ]
    ann2 = [(j, i, text(string(round(ratio[i, j], digits=1)), 10, :black))
      for i in 1:length(rlevel) for j in 1:length(klevel) if ratio[i, j] <= 1
    ]
    annotate!(ann)
    annotate!(ann2)
    plot!(yguidefontsize=22, ytickfontsize=14, xguidefontsize=22, xtickfontsize=14)
  end
end

h = time_compare_heatmap(2)
savefig(h, "plot/time_compare_heatmap_ESGS.pdf")

h = time_compare_heatmap(3)
savefig(h, "plot/time_compare_heatmap_PCLP.pdf")



# latex table
algid = 5;
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