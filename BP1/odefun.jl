const year_seconds = 31556926
const sim_years = 300

# using OrdinaryDiffEq
# using DiffEqCallbacks
using DifferentialEquations
using Printf
using Plots

function odefun(dψV, ψδ, p, t)
  reject_step = p.reject_step
  Vp = p.Vp
  lop = p.lop
  F = p.F
  u = p.u
  Δτ = p.Δτ
  ge = p.ge
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τz0 = p.τz0
  RSDc = p.RSDc
  RSf0 = p.RSf0
  LFtoB = p.LFtoB


  if reject_step[1]
    return
  end
  δNp = div(length(ψδ), 2)
  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[ δNp .+ (1:δNp) ]


  bc_Dirichlet = (lf, x, y, e) -> (2-lf)*(δ ./ 2) + (lf-1)*fill(t * Vp/2, size(x))
  bc_Neumann   = (lf, x, y, nx, ny, e) -> zeros(size(x))


  # solve for displacements everywhere in domain
  e = 1
  locbcarray_mod!(ge, lop[e], LFtoB, bc_Dirichlet, bc_Neumann,(e))
  u[:] = F \ ge


  # set up rates of change for  state and slip
  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[δNp .+ (1:δNp) ]

  dψ .= 0
  V  .= 0


  # Update the fault data
  Δτ .= 0
  lf1 = 1


  Δτ .= -μshear .* computetraction_mod(lop[1], lf1, u, δ)

  for n = 1:δNp
    ψn = ψ[n]
    an = RSa[n]

    τn = Δτ[n] + τz0
    if isnan(τn)
      println("τ reject")
      reject_step[1] = true
      return
    end

    VR = abs(τn / η)
    VL = -VR
    Vn = V[n]
    obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
    (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                 atolx = 1e-9, rtolx = 1e-9)

           
    if isnan(Vn) || iter < 0
      println("V reject")
      reject_step[1] = true
      return
          #error()
    end
    V[n] = Vn  # save this value

    dψ[n] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0)
    if !isfinite(dψ[n])
      println("ψ reject")
      dψ[n] = 0
      reject_step[1] = true
      return
    end
  end

  nothing
end


struct ODE_results
   t_list::Array{Any,1}
   V_list::Array{Any,1}
   δ_list::Array{Any,1}
end

function saveslip(ψδ,t,i,ODEresults,p,base_name="",tdump=100)
    Vmax = 0.0
    # t_list = []
    # V_list = []
    if isdefined(i,:fsallast)
      δNp = div(length(ψδ),2)
      dψV = i.fsallast
      ψ = @view dψV[1:δNp]
      V = @view dψV[δNp .+ (1:δNp)]
      Vmax = maximum(abs.(extrema(V)))
      δ = @view ψδ[δNp .+ (1:δNp)]
      # @show δ
      push!(ODEresults.t_list,t)
      push!(ODEresults.V_list,Vmax)
      push!(ODEresults.δ_list,copy(δ))
      stations = Integer.(range(1,stop=δNp,length=81))
      # @show stations
      # @show typeof(ODEresults.δ_list)
      # if (length(ODEresults.t_list) == 10)
      if (t == sim_years * year_seconds)
        open("$(base_name)V.dat","w") do f
          write(f,"t V ")
          for i in stations
            # write(f,"$(stations[i])")
            write(f,"$((i-1)/(δNp-1)*40000) ")
          end
          write(f,"\n")
        for n = 1:length(ODEresults.t_list)
        # for n = 1:3
          # write(f, "$(ODEresults.t_list[n]) $(ODEresults.V_list[n]) \n")
          # write(f, "$(ODEresults.t_list[n]) $(log10(ODEresults.V_list[n])) $(ODEresults.δ_list[n])")
          write(f, "$(ODEresults.t_list[n]) $(log10(ODEresults.V_list[n]))")
          # for i in 1:div(length(V)+1,2)
          #   write(f,"$(ODEresults.δ_list[n][2*i-1])")
          # end

          for i in Integer.(range(1,81,length=81))
              write(f," $(ODEresults.δ_list[n][i])")
          end
          write(f,"\n")
        end
        # write(f,ODEresults.δ_list)
        @show δ
        @show ODEresults.δ_list
      end
      # @show ODEresults.δ_list
    end
    end
    Vmax
end