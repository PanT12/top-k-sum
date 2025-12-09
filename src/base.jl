mutable struct point
  val::Float64 # value
  m::Int64 # the number of elements larger than value
  s::Float64 # the sum of these elements
  f::Float64 # F_r(u)
  g::Float64 # G_r(u)
end



# return the number and the sum of the elements in I(u, upper)
function num_and_sum(
    a::AbstractVector{Tfr}, lower::Tfr, upper::Tfr=Inf
)where {Tfr<:Union{AbstractFloat,Rational}}
  m = 0
  sum_element = 0.0
  @inbounds @simd for xi in a
    if upper > xi >= lower
      m += 1
      sum_element += xi
    end
  end
  return m, sum_element
end


# return the solution
function solution_form!(
  xstar::AbstractVector{Tfr}, ustar::Tfr, lstar::Tfr, a::AbstractVector{Tfr}
)where {Tfr<:Union{AbstractFloat,Rational}}
  xstar .= a
  λstar = ustar - lstar
  @inbounds @simd for i in eachindex(a)
    val = a[i]
    if val > ustar
      xstar[i] = val - λstar
    elseif lstar <= val <= ustar
      xstar[i] = lstar
    else
      xstar[i] = val
    end
  end
end

#=================================================
G-searching with refinement
input:
  - a: input vector
  - f: F_r(u)
  - u: current u
  - H: k-m
  - init: apply whether filter! or filter. true: filter; false: !filter
  - hist: whether print some details
output:
  - ρ: G_r(u)
  - vlength+add_number: the number of elements in I(G_r(u), u)
=================================================#

function project_topksum_g_searching(
  a::AbstractVector{Tfr},
  f::Tfr, u::Tfr, H::Int,
  init::Bool = true,
  hist::Bool = false
) where {Tfr<:Union{AbstractFloat,Rational}}

  if H == 0
    return maximum(@views a[a.<u]), 1
  end

  # initialization
  solved = false
  vlength, sumv1 = num_and_sum(a, f, u)
  
  if hist
    # println("init vlength is $(vlength), init sum is $(sumv1). time used for calculating sum is $(time).")
  end
  sumv = 0.0;
  add_number = 0;
  right = false; # determine left or right
  @inbounds if vlength < H+1
    right = true # then it must be in the right, we need to add some elements such that vlength = H+1
    if init
      a = filter(x -> x < f, a)
    else
      filter!(x -> x < f, a)
    end

    add_number = H + 1 - vlength;
    sumv = 0.0
    @inbounds @simd for i in 1:add_number
      sumv += a[i]
    end
  end
  ρ = (sumv1 + sumv - H*u)/(vlength - H + add_number) # initial lower bound
  if hist
    # println("init ρ is $(ρ), sum is $(sumv1), length is $(vlength).")
  end

  if !right
    # we need to determine left or right
    if f > ρ # on the right
      if init
        a = filter(x -> x < f, a)
      else
        filter!(x -> x < f, a) # need add elements and remove elements
      end
      right = true
    else # on the left
      if init
        a = filter(x -> u > x >= f, a)
      else
        filter!(x -> u > x >= f, a) # only need remove elements 
      end
      add_number = vlength
      vlength = 0;
    end
  end

  # add elements
  if right
    @inbounds @fastmath for i = (add_number+1):length(a)
      if a[i] > ρ
        a[add_number += 1] = a[i]
        ρ += (a[i] - ρ) / (vlength - H + add_number)
      end
    end
  end
  if hist
    # println("1. rho is $(ρ) and add_number is $(add_number).")
  end
  time3 = @elapsed begin
  # remove elements
  while true
    vlengthold = add_number
    add_number = 0
    @fastmath for i = 1:vlengthold
      if a[i] > ρ - 1e-8
        @inbounds a[add_number += 1] = a[i]
      else
        bottom = (vlengthold - i + add_number + vlength - H);
        # if bottom == 0
        #   ρ = a[i]
        #   add_number = 1
        #   solved = true
        #   break
        # end
        ρ += (ρ - a[i]) / bottom
      end
      # println("rho is $(ρ) and add_number is $(add_number).")
    end
    if add_number == vlengthold
      solved = true
      break
    end
  end
  end

  if hist
    # println("time1 is $(time1), time2 is $(time2), time3 is $(time3).")
  end

  if solved
    return ρ, vlength+add_number
  end

end

#=================================================
Initialization step
input:
  - a: input vector
  - r: scalar budget parameter
  - k: integer
  - hist: whether print some details
output:
  - Ite_init: The number of iteration in the initialization step
  - u: u^i (>= a_k)
  - unew: u^{i+1} (< a_k)
  - flag: 0: Top-k-sum is less than r; 
          -1: Top-k-sum is larger than r. a_{k+1} <= unew < a_k
          -2: Top-k-sum is larger than r.
=================================================#
function initialization!(
  a::AbstractVector{Tfr}, r::Tfr, k::Ti,
  hist::Bool=false, 
)where {Tfr<:Union{AbstractFloat,Rational}, Ti<:Integer}
  
  tol::Tfr = 1e-8 # stopping criteria
  l_fold = unew = r/k;
  max_a = uold = maximum(a)
  Ite_init = 1;
  flag = -2;
  sum_a_r1_old::Tfr = 0.0
  m_old::Ti = 0 
  if hist
    println("-----------------Initialization step begins-----------------")
  end
  while true
    m, sum_a_r1 = num_and_sum(a, unew)
    l_fnew = (r - sum_a_r1 + m*unew)/k

    if unew - l_fnew < tol
      flag = 0;
      if hist
        println("Flag is 0. Top-k-sum is less than r.")
      end
      break
    end

    if m == k
      flag = -1;
      if hist
        println("Flag is -1. Top-k-sum is larger than r. Need to test case 3.")
      end
      break
    end

    if m > k
      if hist
        println("Top-k-sum is larger than r.")
      end
      break
    end

    Ite_init += 1;
    uold = unew;
    l_fold = l_fnew
    m_old = m;
    sum_a_r1_old = sum_a_r1

    unew = (uold*m - k*l_fold)/(m-k);

    if hist
      println("Ite $(Ite_init): u is $(uold), F_r(u) is $(l_fold).")
    end
  end
  u = point(uold, m_old, sum_a_r1_old, l_fold, Inf)
  return Ite_init, u, unew, flag, max_a
end




#=================================================
select initial points for searching u_r^*
input:
  - a: input vector
  - r: scalar budget parameter
  - k: integer
  - Ite_init: The number of iteration in the initialization step
  - unew: u^{i+1} (< a_k)
  - uold: u^i (>= a_k)
  - flag: the value of flag obtained in the initialization step
  - coninit: select initial points consistently
  - hist: whether print some details
output:
  - flag: -1, 1, 2, 3, 10
  - diffL: G_r(uL) - F_r(uL) (>0) if uL(= uold) exists
  - vlengthL: the number of elements in I(G_r(uL), uL) if uL(= uold) exists
  - uC
  - uL
  - uR
=================================================#
function select_initial_points(
  a::AbstractVector{Tfr}, r::Tfr, k::Ti, max_a, Ite_init::Ti, unew::Tfr,
  uold::point, flag::Ti,
  coninit::Bool=false, hist::Bool=false, 
)where {Tfr<:Union{AbstractFloat,Rational}, Ti<:Integer}

  vlengthL = 0
  diffL = Inf
  tol = 1e-10
  uL, uR = point(-Inf, 0, 0.0, 0.0, 0.0), point(Inf, 0, 0.0, 0.0, 0.0)

  if flag == -1
    ak = minimum(a[a.>=unew])
    num_ak, sum_akf = num_and_sum(a, ak)
    Fr_ak = (r - sum_akf + num_ak * ak)/k
    num_Fr_ak, sum_Fr_ak = num_and_sum(a, Fr_ak)
    compare_result = (r - k * Fr_ak) - (sum_Fr_ak - num_Fr_ak * Fr_ak) + k * (ak - Fr_ak);    
    if abs(compare_result) < tol
      uL.val = Fr_ak; uR.val = ak # uL.val denotes l*, uR.val denotes u*
      return flag, diffL, vlengthL, Inf, uL, uR
    end
  end

  if coninit
    l_fold = r/k;
    num_f, sum_f = num_and_sum(a, l_fold)
    compare_result = - (sum_f - num_f * l_fold) + k * (max_a - l_fold)
    if compare_result < 0
      uL.val = l_fold; # uL.val denotes l*, uR.val denotes u*
      flag = -1;
      return flag, diffL, vlengthL, Inf, uL, uR
    end
    l_gold, vlength = project_topksum_g_searching(a, l_fold, max_a + 1e-8, k, true, hist) 
    uR = point(max_a + 1e-8, 0, 0.0, l_fold, l_gold)
    slope = -k/(vlength-k);
    uC = (uR.f - uR.g)/slope + uR.val
    if hist
        println("-----------------Pivot step begins-----------------")
        println("Iteration 0: g is $(l_gold), f is $(l_fold) with difference f - g is $(l_fold - l_gold). Current u is $(max_a)")
        println("-----------------------------------------------")
    end
    flag = 10;
    return flag, diffL, vlengthL, uC, uL, uR
  end


  num_f, sum_f = num_and_sum(a, uold.f)
  compare_result = (r - k * uold.f) - (sum_f - num_f * uold.f) + k * (uold.val - uold.f)
  if Ite_init == 1
    if compare_result < 0
      uL.val = r/k; # uL.val denotes l*, uR.val denotes u*
      flag = -1;
      return flag, diffL, vlengthL, Inf, uL, uR
    end
    flag = 1;
    l_gold, vlength = project_topksum_g_searching(a, uold.f, uold.val + 1e-8, k, true, hist) 
    uR = point(uold.val + 1e-8, uold.m, uold.s, uold.f, l_gold)
    slope = -k/(vlength-k);
    uC = (uR.f - uR.g)/slope + uR.val
  else
    # control whether separate
    if compare_result < -tol # left
      flag = 2;
      l_gold, vlength = project_topksum_g_searching(a, uold.f, uold.val, k-uold.m, true, hist);
      uL = point(uold.val, uold.m, uold.s, uold.f, l_gold)
      uC = ((k*(uL.g - uL.f))/uL.m + 2*uL.val)/2;
      vlengthL = vlength
      diffL = uL.g - uL.f;
    else # right
      flag = 3;
      l_gold, _ = project_topksum_g_searching(a, uold.f, uold.val, k-uold.m, true, hist) 
      uR = point(uold.val, uold.m, uold.s, uold.f, l_gold)
      intp = (k*(uR.g-uR.f))/uR.m + uR.val
      uC = (uold.val + max(intp, unew))/2;
    end
  end
  if hist
      println("Flag is $(flag).")
      println("-----------------Pivot step begins-----------------")
      println("Iteration 0: l_g is $(l_gold), l_f is $(uold.f) with difference l_f - l_g is $(uold.f - l_gold). Current u is $(uold.val)")
      println("-----------------------------------------------")
  end

  return flag, diffL, vlengthL, uC, uL, uR
end



#=================================================
Pivot step
input:
  - a: input vector
  - r: scalar budget parameter
  - k: integer
  - lb_in_init: lower bound (u^{i+1}) in the initialization step.
  - uL
  - uR
  - uC
  - f_subseq: a_f: the filtered sequence used to calculate F_r(u)
  - g_subseq: a_g: the filtered sequence used to calculate G_r(u)
  - diffL: G_r(uL) - F_r(uL) (>0) if uL(= uold) exists in the last step
  - vlengthL: the number of elements in I(G_r(uL), uL) if uL(= uold) exists in the last step
  - coninit: select initial points consistently
  - debug: check the calculation of G_r(u)
  - hist: whether print some details
output:
  - u_L
  - diffL: G_r(u_L) - F_r(u_L) (>0) if u_L(= uold) exists
  - vlengthL: the number of elements in I(G_r(u_L), u_L) if u_L(= uold) exists
  - g_subseq: the filtered sequence used to calculate G_r(u)
  - Ite: the number of iteration in the pivot step
=================================================#
# function pivot0(
#   a::AbstractVector{Tfr}, r::Tfr, k::Ti, lb_in_init::Tfr,
#   uL::point, uR::point, uC::Tfr,
#   f_subseq::AbstractVector{Tfr}, g_subseq::AbstractVector{Tfr},
#   diffL::Tfr, vlengthL::Ti, coninit::Bool, 
#   debug::Bool=false, hist::Bool=false
# )where {Tfr<:Union{AbstractFloat,Rational}, Ti<:Integer}

#   tol = 1e-8
#   m=0;
#   Ite = 1;
#   solved = false;
#   sum_a_r1 = +Inf
#   sum_r1 = 0; 
#   while Ite<=100
#     sum_r1, sum_a_r1 = num_and_sum(f_subseq, uC)
#     m = sum_r1 + uR.m # calculate I(a, unew)
#     if m > k || (!coninit && uC <= lb_in_init)  # useless point
#       if hist
#         println("Skip $(uC) and take uC=$((uC + uR.val) / 2).")
#       end
#       uC = (uC + uR.val) / 2
#       Ite+=1
#       continue
#     end
#     # f searching 
#     sum_larger_u = sum_a_r1 + uR.s;
#     uCf = (r - sum_larger_u + m * uC)/k
#     # g searching
#     uCg, vlength = project_topksum_g_searching(collect(g_subseq), uCf, uC, k-m, false, hist)
#     if debug
#       error = (k-sum(a.>=uC))*(uC-uCg) - sum(max.(a[a.<uC].-uCg, 0))
#       if abs(error) > 1e-5
#         @warn("debuging. error is $(error).")
#       end
#     end
#     if hist
#       println("Iteration $(Ite): m is $(m), g is $(uCg), f is $(uCf) with difference f - g is $(uCf - uCg). Current u is $(uC), vlength is $(vlength)")
#     end

#     diff = uCg - uCf;
#     if abs(diff) < tol
#       uL = point(uC, m, sum_larger_u, uCf, uCg)
#       diffL = diff;
#       solved = true;
#       break
#     end
#     # filter g_subseq or f_subseq
#     if diff < 0
#       g_subseq = view(g_subseq, uCg-tol .<= g_subseq .< uC)
#       f_subseq = view(f_subseq, f_subseq .< uC)
#     else
#       # left hand side
#       f_subseq = view(f_subseq, f_subseq .>= uC)
#       vlengthL = vlength
#       diffL = diff;
#     end

#     # the case updating uC
#     if hist
#       println("uR is $(uR) and uL is $(uL)")
#     end
#     if uR.val == +Inf
#       if diff < 0
#         uR = point(uC, m, sum_larger_u, uCf, uCg)
#       else
#         uLB = deepcopy(uL)
#         uL = point(uC, m, sum_larger_u, uCf, uCg)
#       end
#     elseif uL.val == -Inf
#       if diff > 0
#         uL = point(uC, m, sum_larger_u, uCf, uCg)
#       else
#         uRB = deepcopy(uR)
#         uR = point(uC, m, sum_larger_u, uCf, uCg)
#       end
#     else
#       if hist
#         println("m - m_uL is $(m - uL.m). m - m_uR is $(m - uR.m).")
#       end

#       uCold = uC
#       if m == uL.m
#         if diff <= 0
#           break
#         else
#           uL = point(uC, m, sum_larger_u, uCf, uCg)
#           uC = minimum(@views f_subseq[f_subseq .> uC])
#           uCf = uCf + m/k * (uC - uCold)
#           if hist
#             println("jump to a larger point $(uC) and ucf is $(uCf).")
#           end
#           uCg, vlength = project_topksum_g_searching(collect(g_subseq), uCf, uC, k-m, false, hist)
#           if debug
#             error = (k-sum(a.>=uC))*(uC-uC) - sum(max.(a[a.<uC].-uC, 0))
#             if abs(error) > 1e-5
#               @warn("debuging. error is $(error).")
#             end
#           end
#           diff = uCg - uCf
#           if diff <= 0
#             g_subseq = view(g_subseq, uCg-tol .<= g_subseq .< uC)
#             break
#           else
#             uL = point(uC, m, sum_larger_u, uCf, uCg)
#           end
#         end
#       elseif m == uR.m
#         if diff >= 0
#           uL = point(uC, m, sum_larger_u, uCf, uCg)
#           vlengthL = vlength;
#           break
#         else
#           uC = maximum(@views f_subseq[f_subseq .< uC])
#           if hist
#             println("jump to a smaller point $(uC).")
#           end
#           uCf = uCf + m/k * (uC - uCold)
#           uCg, vlength = project_topksum_g_searching(collect(g_subseq), uCf, uC, k-m-1, false, hist)
#           if debug
#             error = (k-sum(a.>uC))*(uC-uC) - sum(max.(a[a.<=uC].-uC, 0))
#             if abs(error) > 1e-5
#               @warn("debuging. error is $(error).")
#             end
#           end
#           diff = uCg - uCf
#           if diff >= 0
#             uL = point(uC, m, sum_larger_u, uCf, uCg) # note: take the right-derivative
#             vlengthL = vlength + 1; diffL = diff; # note: take the right-derivative
#             break
#           else
#             uR = point(uC, m+1, sum_larger_u + uC, uCf, uCg)
#             f_subseq = view(f_subseq, f_subseq .< uC)
#           end
#         end
#       else
#         if diff < 0
#           uR = point(uC, m, sum_larger_u, uCf, uCg)
#         else
#           uL = point(uC, m, sum_larger_u, uCf, uCg)
#         end
#       end
#     end

#     if hist
#       println("-----------------------------------------------")
#     end

#     # generate unew
#     Ite += 1;
#     if uR.val == +Inf
#       uC = generate_newpoint(uLB, uL)
#     elseif uL.val == -Inf
#       uC = generate_newpoint(uRB, uR)
#     else
#       uC = generate_newpoint(uR, uL)
#     end 

#   end

#   return uL, vlengthL, diffL, g_subseq, Ite

# end

function pivot(
  a::AbstractVector{Tfr}, r::Tfr, k::Ti, lb_in_init::Tfr,
  uL::point, uR::point, uC::Tfr,
  f_subseq::AbstractVector{Tfr}, g_subseq::AbstractVector{Tfr},
  diffL::Tfr, vlengthL::Ti, coninit::Bool, 
  debug::Bool=false, hist::Bool=false
)where {Tfr<:Union{AbstractFloat,Rational}, Ti<:Integer}

  tol = 1e-8
  m=0;
  Ite = 1;
  solved = false;
  sum_a_r1 = +Inf
  sum_r1 = 0; 

  while Ite<=100
    sum_r1, sum_a_r1 = num_and_sum(f_subseq, uC)
    m = sum_r1 + uR.m # calculate I(a, unew)
    if m > k || (uC <= lb_in_init)  # useless point
      if hist
        println("Skip $(uC) and take uC=$((uC + uR.val) / 2).")
      end
      uC = (uC + uR.val) / 2
      Ite+=1
      continue
    end
    # f searching 
    sum_larger_u = sum_a_r1 + uR.s;
    uCf = (r - sum_larger_u + m * uC)/k
    # g searching
    uCg, vlength = project_topksum_g_searching(collect(g_subseq), uCf, uC, k-m, false, hist)
    if debug
      error = (k-sum(a.>=uC))*(uC-uCg) - sum(max.(a[a.<uC].-uCg, 0))
      if abs(error) > 1e-5
        @warn("debuging. error is $(error).")
      end
    end
    if hist
      println("Iteration $(Ite): m is $(m), g is $(uCg), f is $(uCf) with difference f - g is $(uCf - uCg). current u is $(uC).")
    end

    diff = uCg - uCf;
    if abs(diff) < tol
      uL = point(uC, m, sum_larger_u, uCf, uCg)
      diffL = diff;
      solved = true;
      break
    end
    # filter g_subseq or f_subseq
    if diff < 0
      g_subseq = view(g_subseq, uCg-tol .<= g_subseq .< uC)
      f_subseq = view(f_subseq, f_subseq .< uC)
    else
      # left hand side
      f_subseq = view(f_subseq, f_subseq .>= uC)
      vlengthL = vlength
      diffL = diff;
    end

    if hist
      println("uR is $(uR) and uL is $(uL). m - m_uL is $(m - uL.m). m - m_uR is $(m - uR.m).")
    end

    uCold = uC
    if m == uL.m
      if diff <= 0
        break
      else
        uL = point(uC, m, sum_larger_u, uCf, uCg)
        uC = minimum(@views f_subseq[f_subseq .> uC])
        uCf = uCf + m/k * (uC - uCold)
        if hist
          println("jump to a larger point $(uC) and ucf is $(uCf).")
        end
        uCg, vlength = project_topksum_g_searching(collect(g_subseq), uCf, uC, k-m, false, hist)
        if debug
          error = (k-sum(a.>=uC))*(uC-uC) - sum(max.(a[a.<uC].-uC, 0))
          if abs(error) > 1e-5
            @warn("debuging. error is $(error).")
          end
        end
        diff = uCg - uCf
        if diff <= 0
          g_subseq = view(g_subseq, uCg-tol .<= g_subseq .< uC)
          break
        end
      end
    elseif m == uR.m
      if diff >= 0
        uL = point(uC, m, sum_larger_u, uCf, uCg)
        vlengthL = vlength;
        break
      else
        uC = maximum(@views f_subseq[f_subseq .< uC])
        if hist
          println("jump to a smaller point $(uC).")
        end
        uCf = uCf + m/k * (uC - uCold)
        cur_uC_num = sum(uC .== f_subseq)
        uCg, vlength = project_topksum_g_searching(collect(g_subseq), uCf, uC, k-m-cur_uC_num, false, hist)
        if debug
          error = (k-sum(a.>uC))*(uC-uC) - sum(max.(a[a.<=uC].-uC, 0))
          if abs(error) > 1e-5
            @warn("debuging. error is $(error).")
          end
        end
        diff = uCg - uCf
        if diff >= 0
          uL = point(uC, m, sum_larger_u, uCf, uCg) # note: take the right-derivative
          vlengthL = vlength + cur_uC_num; diffL = diff; # note: take the right-derivative
          break
        else
          m += cur_uC_num; sum_larger_u += cur_uC_num * uC; 
          f_subseq = view(f_subseq, f_subseq .< uC)
        end
      end
    end

    if uR.val == +Inf
      if diff < 0
        uR = point(uC, m, sum_larger_u, uCf, uCg)
      else
        uLB = deepcopy(uL)
        uL = point(uC, m, sum_larger_u, uCf, uCg)
      end
    elseif uL.val == -Inf
      if diff > 0
        uL = point(uC, m, sum_larger_u, uCf, uCg)
      else
        uRB = deepcopy(uR)
        uR = point(uC, m, sum_larger_u, uCf, uCg)
      end
    else
      if diff < 0
        uR = point(uC, m, sum_larger_u, uCf, uCg)
      else
        uL = point(uC, m, sum_larger_u, uCf, uCg)
      end
    end

    if hist
      println("-----------------------------------------------")
    end

    # generate unew
    Ite += 1;
    if uR.val == +Inf
      uC = generate_newpoint(uLB, uL)
    elseif uL.val == -Inf
      uC = generate_newpoint(uRB, uR)
    else
      uC = generate_newpoint(uR, uL)
    end 

  end

  return uL, vlengthL, diffL, g_subseq, Ite

end

# updating rule of uC
function generate_newpoint(x::point, y::point)
  return (y.val * (x.g - x.f) - x.val * (y.g - y.f)) / ((x.g - x.f) - (y.g - y.f))
end


#=================================================
Exact step
input:
  - a: input vector
  - r: scalar budget parameter
  - k: integer
  - uL
  - g_subseq: a_g: the filtered sequence used to calculate G_r(u)
  - diffL: G_r(uL) - F_r(uL) (>0) if uL(= uold) exists in the last step
  - vlengthL: the number of elements in I(G_r(uL), uL) if uL(= uold) exists in the last step
  - debug: check the calculation of G_r(u)
  - hist: whether print some details
output:
  - u_L: u_r^*
  - Ite: the number of iteration in the exact step
=================================================#
function exact(
  a::AbstractVector{Tfr}, r::Tfr, k::Ti, 
  uL::point, g_subseq::AbstractVector{Tfr},
  diffL::Tfr, vlengthL::Ti,
  debug::Bool=false, hist::Bool=false
)where {Tfr<:Union{AbstractFloat,Rational}, Ti<:Integer}

  if hist
    println("-----------------Exact step begins-----------------")
    println("uL is $(uL), diff is $(diffL), m is $(uL.m)")
  end
  m = uL.m
  tol = 1e-8
  Ite = 0;
  while abs(diffL) > tol
    Ite += 1
    uLoldval = uL.val
    s = -m/k + (m-k)/(vlengthL - (k-m))
    uL.val = -diffL/s + uLoldval
    # println("s is $(s), vlengthL is $(vlengthL)")
    uL.f = uL.f + m/k * (uL.val - uLoldval)
    # println("uL.val is $(uL.val), uL.f is $(uL.f), g_subseq is $(g_subseq).")
    uL.g, vlengthL = project_topksum_g_searching(g_subseq, uL.f, uL.val, k-m, true, hist)
    if hist
      println("Iteration $(Ite): m is $(m), g is $(uL.g), f is $(uL.f) with difference f - g is $(uL.f - uL.g). Current u is $(uL.val)")
      println("-----------------------------------------------")
    end
    if debug
      error = (k-sum(a.>=uL.val))*(uL.val-uL.g) - sum(max.(a[a.<uL.val].-uL.g, 0))
      if abs(error) > 1e-5
        @warn("debuging. error is $(error).")
      end
    end
    diffL = uL.g - uL.f
    # if Ite > 10
    #   break
    # end
  end

  # result
  if m == k
    uL.val = minimum(@views a[a.>=uL.val])
    uL.m, sum_a_r1 = num_and_sum(a, uL.val)
    uL.f = (r - sum_a_r1 + uL.m * uL.val)/k
    if hist
      println("Fake point! It should be $(uL.val), $(uL.f).")
    end
  end
  return Ite, uL
end

