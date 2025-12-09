#=================================================
EIPS: detailed timing of Algorithm EIPS
input:
  - xstar: preallocated output vector
  - a: input vector
  - r: scalar budget parameter
  - k: integers
  - hist: binary for history
output:
  - iter_num: (Ite_init, Ite_pivot, Ite_exact)
=================================================#
include("base.jl")
function project_topksum_EIPS!(
  xstar::AbstractVector{Tfr}, a::AbstractVector{Tfr}, 
  r::Tfr, k::Ti, hist::Bool=false,
)where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}

  len = length(a);
  tol::Tfr = Tfr(1e-8)
  
  if k == 1
    flag = -2
    xstar .= min.(r, a)
    if hist
      return (0, 0, 0)
    end
  end
  if k == len
    term = (r - max(r, sum(a)))/len
    flag = -2
    xstar .= a .+ term
    if hist
      return (0, 0, 0)
    end
  end

  Ite_init, u, unew, flag, max_a = initialization!(a, r, k, hist)

  if flag == 0
    xstar .= a
    return (Ite_init, 0, 0)
  end

  flag, diffL, vlengthL, uC, uL, uR = select_initial_points(
    a, r, k, max_a, Ite_init, unew, u, flag, false, hist
  )


  if flag == -1
    solution_form!(xstar, uR.val, uL.val, a)
    if hist
      return (Ite_init, 0, 0)
    end
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

  uL, vlengthL, diffL, g_subseq, Ite_pivot = pivot(a, r, k, unew, uL, uR, uC, f_subseq, g_subseq, diffL, vlengthL, false, false, hist)

  Ite_exact, uL = exact(a, r, k, uL, g_subseq, diffL, vlengthL, false, hist)
  solution_form!(xstar, uL.val, uL.f, a)
  if hist
    return (Ite_init, Ite_pivot, Ite_exact)
  end
end