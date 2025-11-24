using Gurobi, JuMP

function default_proj_options()
  options = Dict()
  options[:maxtime] = 3000 # seconds
  options[:method] = 2 # barrier
  options[:maxiter] = 10_000
  options[:fea_tol] = 1e-9
  options[:opt_tol] = 1e-9
  options[:bar_tol] = 1e-9
  options[:nthreads] = 8
  return options
end

function project_topksum_guro!(
  xbarsort::Vector, sig::Vector, x0::Vector, r::Real, k::Integer, active::Bool,
  options::Dict=default_proj_options()
)
    
  # initialization
  time_init = @elapsed begin
    n = length(x0)
    if !active
      time_primal = @elapsed xbarsort .= x0
      case = 0;
      solved=true
      return case, 0, (0.0, 0.0, time_primal), solved
    end

    # model = Model(Gurobi.Optimizer)
    # sig = sortperm(x0, rev=true, alg=QuickSort);
    x0sort = x0[sig];
    outputFlag = 0
    model = Model(optimizer_with_attributes(()->Gurobi.Optimizer(), "OutputFlag" => outputFlag, "Threads" => 0)); 
    MOI.set(model, MOI.Silent(), true);
    # variables
    @variable(model, x[1:n])

    # objective function
    @objective(model, Min, sum((x[i] - x0sort[i])^2 for i in 1:n))

    # Define constraints
    indk = [i <= k ? 1.0 : 0.0 for i in 1:n]
    # top-k-sum constraint
    @constraint(model, sum(indk[i] * x[i] for i in 1:n) <= r)

    for i in 1:k-1
      @constraint(model, x[i] >= x[k])
    end

    for j in k+1:n
      @constraint(model, x[j] <= x[k])
    end

    # Set parameters
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "Method", options[:method])
    set_optimizer_attribute(model, "BarIterLimit", options[:maxiter])
    set_optimizer_attribute(model, "FeasibilityTol", options[:fea_tol])
    set_optimizer_attribute(model, "OptimalityTol", options[:opt_tol])
    set_optimizer_attribute(model, "BarConvTol", options[:bar_tol])
    set_optimizer_attribute(model, "Threads", options[:nthreads])
  end

  # Solve the model
  solvetime = @elapsed optimize!(model);
  time_run = solve_time(model)

  # result
  status = termination_status(model)

  solved = false
  time_primal = @elapsed begin
    if status == MOI.OPTIMAL
      solved=true
      xbarsort[sig] .= value.(x)
    else
      @warn("Not founded.")
    end
  end

  return 1, 0, (time_init, time_run, time_primal), solved

end