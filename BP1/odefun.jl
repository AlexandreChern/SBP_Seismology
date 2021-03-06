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
  # τn = p.τn
  τ = p.τ

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
    τ[n] = τn
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


function setupfaultstations(locations)
  T = eltype(locations)
  @assert size(locations,2) == 2
end


struct ODE_results
   t_list::Array{Any,1}
   V_list::Array{Any,1}
   δ_list::Array{Any,1}
  #  station_1::Array{Any,1}
  #  station_2::Array{Any,1}
  #  station_3::Array{Any,1}
  #  station_4::Array{Any,1}
  #  station_5::Array{Any,1}
  #  station_6::Array{Any,1}
  #  station_7::Array{Any,1}
  #  station_8::Array{Any,1}
  #  station_9::Array{Any,1}
  #  station_10::Array{Any,1}
  #  station_11::Array{Any,1}
  stations::Dict{Int64,Array{Any,1}}
end

function saveslip(ψδ,t,i,ODEresults,yf,stations,station_indices,p,base_name="",tdump=100)
    Vmax = 0.0
    # t_list = []
    # V_list = []
    if isdefined(i,:fsallast)
      δNp = div(length(ψδ),2)
      dψV = i.fsallast
      dψ = @view dψV[1:δNp]
      V = @view dψV[δNp .+ (1:δNp)]
      Vmax = maximum(abs.(extrema(V)))
      δ = @view ψδ[δNp .+ (1:δNp)]
      ψ = @view ψδ[1:δNp]
      # @show δ
      push!(ODEresults.t_list,t)
      push!(ODEresults.V_list,Vmax)
      push!(ODEresults.δ_list,copy(δ))
      # push!(ODEresults.station_1,V[station_indices[1]])
      # push!(ODEresults.station_2,V[station_indices[2]])
      # push!(ODEresults.station_3,V[station_indices[3]])
      # push!(ODEresults.station_4,V[station_indices[4]])
      # push!(ODEresults.station_5,V[station_indices[5]])
      # push!(ODEresults.station_6,V[station_indices[6]])
      # push!(ODEresults.station_7,V[station_indices[7]])
      # push!(ODEresults.station_8,V[station_indices[8]])
      # push!(ODEresults.station_9,V[station_indices[9]])
      # push!(ODEresults.station_10,V[station_indices[10]])
      # push!(ODEresults.station_1,V[station_indices[11]])
      # push!(ODEresults.stations[1],V[station_indices[1]])
      # push!(ODEresults.stations[2],V[station_indices[2]])
      for i =1:11
        # @show p.τ[station_indices[i]]
        push!(ODEresults.stations[i],[V[station_indices[i]],p.τ[station_indices[i]],p.RSDc * exp((ψ[station_indices[i]] - p.RSf0) / p.RSb) /
        p.RSV0])
      end
      # stations = Integer.(range(1,stop=δNp,length=81))
      # @show stations
      # @show typeof(ODEresults.δ_list)
      # if (length(ODEresults.t_list) == 10)
      if (t == sim_years * year_seconds)
        @show station_indices
        for i in range(1,stop=length(station_indices))
          station_name = Symbol("station_",i)
          # @show station_name
          station_id = station_indices[i]
          open("$(base_name)$(stations[i]).dat","w") do f
            write(f,"t slip slip_rate shear_stress state\n")
            for n = 1:length(ODEresults.t_list)
              # write(f,"$(ODEresults.t_list[n]) $(ODEresults.δ_list[n][station_id]) $(log10(abs(ODEresults.(eval(station_name[i])))))\n")
              write(f,"$(ODEresults.t_list[n]) $(ODEresults.δ_list[n][station_id]) $(log10(abs(ODEresults.stations[i][n][1]))) $(ODEresults.stations[i][n][2]) $(log10(ODEresults.stations[i][n][3]))\n")
            end
            # for n = 1:length(ODE_results.t_list)
            #   write(f,"$(ODEresults.t_list[n]) $(ODEresults.δ_list[n][station_id]")
            # end
          end
        end
      end
      if (t == sim_years * year_seconds)
        open("$(base_name)V.dat","w") do f
          write(f,"z \n")
          write(f,"t Slip_rate Slip")
          write(f,"\n") 
          write(f,"0.0 0.0 ")
          # for i in stations
          for i in 1:δNp
            # write(f,"$(stations[i])")
            # write(f,"$(Integer((i-1)*40000/(δNp-1))) ")
            write(f,"$(yf[i]) ")
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

          # for i in Integer.(range(1,81,length=81))
          for i in 1:δNp
              # write(f," $(ODEresults.δ_list[n][i])")
              write(f," $(ODEresults.δ_list[n][i])")
          end
          write(f,"\n")
        end
        # write(f,ODEresults.δ_list)
        # @show δ
        # @show ODEresults.δ_list
      end
      # @show ODEresults.δ_list
    end
    end
    Vmax
end