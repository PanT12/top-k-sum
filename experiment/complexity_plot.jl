using LinearAlgebra, Random, DelimitedFiles, Plots, LaTeXStrings, GLM, DataFrames
# using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics, Parsers, Plots, Measures
# using SparseMatricesCSR, SparseArrays
global const PROJPATH = pwd();
include(PROJPATH * "/src/base.jl")
include(PROJPATH * "/experiment/data_process.jl")
global const DATAPATH = PROJPATH * "/complexity"
# global const DATAPATH = PROJPATH * "/complexity_tau_r=1"


rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10];
ri = [2; 6; 10; 13; 15; 17];
ki = [1; 6; 8; 12];
klevel_sub = klevel[ki];
rlevel_sub = rlevel[ri];
expers = [(x, y) for x in rlevel_sub for y in klevel_sub];

r_k_pair = [(x, y) for x in ri for y in ki];

alg_name = ["with init", "without init"];
algid = [1];

# estimate the slope: if the slope ≈ 1: linear behavior
function estimate_slope(data::Dict, n::Array, ki::Int, ri::Int, algid)
  df = DataFrame(x = log10.(n), y = log10.([float(item) for item in data[algid]["total_mean"][(ri, ki)]]));
  model = lm(@formula(y ~ x), df);
  return round(coef(model)[2], digits=2)
end


# plot. 
#=================================================
Note: if r is smaller than the top-k-sum, we only plot the time spent in the initialization step;
      otherwise, we plot the whole time.
=================================================#
function complexity_plot(algid, ri, ki_list, n::Array, data::Dict, left::Bool=false, bottom::Bool=false)
  taur = rlevel[ri]
  key = taur >= 1.0 ? "init_mean" : "total_mean"
  println(key)
  i = 2
  slope = estimate_slope(data, n, ki_list[1], ri, algid)
  p = plot(n, data[algid][key][(ri, ki_list[1])], shape=shape[1], legendfont=10, label="$(slope)", 
    yscale=:log10, xscale=:log10, markersize=6, color=my_color_list[1], markerstrokewidth=0, 
    linewidth=3, legendfontsize=10
  ) 
  while i <= length(ki)
    slope = estimate_slope(data, n, ki_list[i], ri, algid)
    # plot!(n, data[ki[i]], marker=:circle, label="\$\\tau_k=\$$(klevel_sub[i])")
    plot!(n, data[algid][key][(ri, ki_list[i])], shape=shape[i], markersize=6, 
      label="$(slope)", color=my_color_list[i], markerstrokewidth=0, linewidth=3
    )
    i+=1
  end

  if left
    ylabel!("Time (sec)", fontfamily="Times New Roman", yguidefontsize=15)
  end
  if bottom
    xlabel!("n", fontfamily="Times New Roman", xguidefontsize=15)
  end

  title!("Experiment: \$ \\tau_r = $(numerator(taur)) / $(denominator(taur)) \$", fontfamily="Times New Roman")

  xticks = [10^4, 10^5, 10^6, 10^7]
  yticks = [10^(-4), 10^(-3), 10^(-2), 10^(-1), 10^(0)]
  plot!(xticks=xticks, yticks=yticks, ytickfontsize=12, xtickfontsize=12)

  plot!(framestyle=:box)
  plot!(grid=true)
  plot!(legend=:topleft, legendfont="Times New Roman", legendtitle="slope")

  display(plot!)
  return p
end




n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);

figures = [];
for i in eachindex(ri)
  rid = ri[i]
  left, bottom = false, false
  if i in [1, 4]
    left = true
  end
  if i in [4, 5, 6]
    bottom = true
  end
  p = complexity_plot(1, rid, ki, n, data, left, bottom)
  push!(figures, p)
end

plot(figures..., layout = (2, 3), size = (1200, 600),
  left_margin = [5mm 2mm 2mm], bottom_margin = [5mm 0mm]
)

# legend
p0=plot([0, 0, 0, 0]',grid=false, showaxis=false,
  label=["\$\\tau_k = $(numerator(klevel_sub[1])) / $(denominator(klevel_sub[1]))\$" "\$\\tau_k = $(numerator(klevel_sub[2])) / $(denominator(klevel_sub[2]))\$" "\$\\tau_k = $(numerator(klevel_sub[3])) / $(denominator(klevel_sub[3]))\$" "\$\\tau_k = $(numerator(klevel_sub[4])) / $(denominator(klevel_sub[4]))\$"], 
  shape=[shape[1] shape[2] shape[3] shape[4]], 
  color=[my_color_list[1] my_color_list[2] my_color_list[3] my_color_list[4]],
  legend=:outertop, legend_column = -1, size = (1000, 70), 
  legendfontsize = 12, markersize = 6, markerstrokewidth=0
)

lay = @layout([[a b c]; [d e f];g{0.1h}])
p = plot(figures..., p0, layout=lay, size=(1200, 650), 
  # bottommargin=2mm, leftmargin=5mm, padding=-1mm, margin=-2mm
  left_margin=[7mm -2mm -2mm], bottommargin=[2mm 2mm 10mm]
)

savefig(p, "plot/complexity_plot.pdf")










# plot the special case: τ_r = 1, that is r is equal to the top-k-sum
global const DATAPATH = PROJPATH * "/complexity_tau_r=1"
# special case tau_r = 1;
rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]; [100; 101; 110; 120; 150; 200] .// 100];
klevel = [[1; 10; 100; 500] .// 10^4; [1;2;3;4;5;6;7;8;9] .// 10]

# = 1 part 
ri = [13]
ki = [1; 2; 3; 4; 5; 6; 7; 8; 9]

klevel_sub = klevel[ki];
rlevel_sub = rlevel[ri];
expers = [(x, y) for x in rlevel_sub for y in klevel_sub];

r_k_pair = [(x, y) for x in ri for y in ki];

n, data, out = process_data_plot(DATAPATH, algid, r_k_pair);

function complexity_plot_r1(algid, ri, ki_list, n::Array, data::Dict, left::Bool=false, bottom::Bool=false)
  taur = rlevel[ri]
  key = taur >= 1.0 ? "init_mean" : "total_mean"
  println(key)
  i = 2
  slope = estimate_slope(data, n, ki_list[1], ri, algid)
  tauk = klevel[ki_list[1]]
  numer, denom = numerator(tauk), denominator(tauk)
  p = plot(n, data[algid][key][(ri, ki_list[1])], shape=shape[1], legendfont=10,
    yscale=:log10, xscale=:log10, markersize=6, markerstrokewidth=0, 
    linewidth=3, legendfontsize=12, label="\$\\tau_k=\$ $(numer)/$(denom): $(slope)",
    legendtitlefontsize=15, 
    left_margin=2mm, bottom_margin=2mm
  ) 
  while i <= length(ki)
    slope = estimate_slope(data, n, ki_list[i], ri, algid)
    tauk = klevel[ki_list[i]]
    numer, denom = numerator(tauk), denominator(tauk)
    plot!(n, data[algid][key][(ri, ki_list[i])], shape=shape[i], markersize=6, 
      markerstrokewidth=0, linewidth=3, label="\$\\tau_k=\$ $(numer)/$(denom): $(slope)"
    )
    i+=1
  end


  ylabel!("Time (sec)", fontfamily="Times New Roman", yguidefontsize=15)
  xlabel!("n", fontfamily="Times New Roman", xguidefontsize=15)

  title!("Experiment: \$ \\tau_r = $(numerator(taur)) / $(denominator(taur)) \$", fontfamily="Times New Roman")

  xticks = [10^6, 10^7, 10^8]
  yticks = [10^(-2.5), 10^(-2), 10^(-1.5), 10^(-1), 10^(-0.5)]
  plot!(xticks=xticks, yticks=yticks, ytickfontsize=10, xtickfontsize=10)

  plot!(framestyle=:box)
  plot!(grid=true)
  plot!(legend=:outertopright, legendfont="Times New Roman", legendtitle="slope")
  plot!(size = (700, 400))

  display(plot!)
  return p
end


p = complexity_plot_r1(1, ri[1], ki, n, data)

savefig(p, "plot/complexity_r=1_plot.pdf")