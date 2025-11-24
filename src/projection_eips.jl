include("base.jl")
function project_topksum_EIPS!(
  xstar::AbstractVector{Tfr}, a::AbstractVector{Tfr}, 
  r::Tfr, k::Ti, not_select::Bool=false, debug::Bool=false, hist::Bool=false,
)where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}

  ## Test active
  ## initialization
  len = length(a);
  tol::Tfr = Tfr(1e-8)
  
  if k == 1
    flag = -2
    xstar .= min.(r, a)
    return flag, (0, 0, 0)
  end
  if k == len
    term = (r - max(r, sum(a)))/len
    flag = -2
    time_exact = @elapsed xstar .= a .+ term
    return flag, (0, 0, 0)
  end

  # initialization step
  Ite_init, u, unew, flag = initialization!(a, r, k, hist)

  if flag == 0
    xstar .= a
    return flag, (Ite_init, 0, 0)
  end

  flag, diffL, vlengthL, uC, uL, uR = select_initial_points(
    a, r, k, Ite_init, unew, u, flag, not_select, hist
  )

  if flag == 10
    if uR.val != Inf
      g_subseq = view(a, uR.g .<= a)
      f_subseq = view(a, 1:len)
    else
      g_subseq = view(a, 1:len)
      f_subseq = view(a, falses(len))
    end
  elseif flag == -1
    solution_form!(xstar, uR.val, uL.val, a)
    return flag, (Ite_init, 0, 0)
  elseif flag == 1
    g_subseq = f_subseq = view(a, a.>=uR.g-tol)
    f_subseq = view(f_subseq, f_subseq.>=unew)
  elseif flag == 2
    g_subseq = view(a, 1:len)
    f_subseq = view(a, a .>= uL.val)
  elseif flag == 3
    f_subseq = g_subseq = view(a, a .< uR.val)
    g_subseq = view(g_subseq, uR.g-tol .<= g_subseq)
  end

  uL, vlengthL, diffL, g_subseq, Ite_pivot = pivot(
    a, r, k, unew, uL, uR, uC, f_subseq, g_subseq, diffL, vlengthL, not_select, debug, hist
  )

  Ite_exact, uL = exact(
    a, r, k, uL, g_subseq, diffL, vlengthL, debug, hist
  )
  solution_form!(xstar, uL.val, uL.f, a)
  return flag, (Ite_init, Ite_pivot, Ite_exact)
end
