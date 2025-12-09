using LinearAlgebra, Random, DelimitedFiles
using Gurobi, JuMP, SparseMatricesCSR, SparseArrays, Statistics
using SparseMatricesCSR, SparseArrays, Statistics
using Plots, LaTeXStrings
global const PROJPATH = pwd();
include(PROJPATH * "/experiment/projection.jl")


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

n = 10^2;
Random.seed!(2)
x0 = rand(n);
x0sort = sort(x0, rev=true);
# rlevel = [[0;1;2;3;4;5;6;7;8;9] .//10; [99//100, 999//1000]];
rlevel = [collect(1:2:90) .// 100; collect(91:1:99) .// 100];
k = max(2, Int64(ceil(n*0.4)))
tks = sum(x0sort[1:k]);
flag_list = [];
case_list = [];
xstar = similar(x0);
for ri in 1:length(rlevel)
  r = tks * float(rlevel[ri]) # 1
  @views case, flag, Ite, t = project_topksum_EIPS_experiment!(
    xstar, x0, r, k, false, false, false
  )
  push!(flag_list, flag);
  push!(case_list, case);
end

# plot flag distribution
p1 = plot(rlevel, flag_list, marker=:circle, line=:dash,
  xlabel=L"r / T_{(k)}(a)", ylabel="Value", markersize=3,
  title="Flag and case distributions",
  fontfamily="Times New Roman", titlefont=14, label="Flag", 
  legend=:bottomright, legendfontsize=10,
  xguidefontsize=15, yguidefontsize=15,
  markercolor=my_color["blackblue"],
)
plot!(ytickfontsize=10, xtickfontsize=10)
plot!(rlevel, case_list, marker=:utriangle, line=:dot, label="Case",
  markercolor=my_color["red"]
)
plot!(framestyle=:box)




rlevel = collect(950:1:990) .// 1000;
flag_list = [];
case_list = [];
xstar = similar(x0);
for ri in 1:length(rlevel)
  r = tks * float(rlevel[ri]) # 1
  @views case, flag, Ite, t = project_topksum_EIPS_experiment!(
    xstar, x0, r, k, false, false, false
  )
  push!(flag_list, flag);
  push!(case_list, case);
end

# plot flag distribution
p2 = plot(rlevel, flag_list, marker=:circle, line=:dash,
  xlabel=L"r / T_{(k)}(a)", ylabel="Value", markersize=3, 
  title="Flag and case distributions",
  fontfamily="Times New Roman", titlefont=14, label="Flag", 
  legend=:bottomright, legendfontsize=10,
  xguidefontsize=15, yguidefontsize=15,
  markercolor=my_color["blackblue"]
)
plot!(ytickfontsize=10, xtickfontsize=10)
plot!(rlevel, case_list, marker=:utriangle, line=:dot, label="Case",
  markercolor=my_color["red"]
)
plot!(framestyle=:box)


lay = @layout([a b]);
p = plot(p1, p2, layout=lay, size=(1200, 400), 
  bottommargin=7mm, leftmargin=[10mm 10mm], right_margin = [0mm 5mm]
)

savefig(p, "plot/flag_case_distribution.pdf")