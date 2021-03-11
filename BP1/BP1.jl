include("./global_curved.jl")
include("odefun.jl")
using Plots

function main()
  sim_years = 3000.

  Vp = 1e-9 # plate rate
  ρ = 2.670
  cs = 3.464
  σn = 50
  RSamin = 0.01
  RSamax = 0.025
  RSb = 0.015
  RSDc = 0.008 # change it to be 0.008
  RSf0 = 0.6
  RSV0 = 1e-6
  RSVinit = 1e-9
  RSH1 = 15
  RSH2 = 18

  μshear = cs^2 * ρ
  η = μshear / (2 * cs)


  # N + 1 is the number of grid points in each dimension
  # N = 200
  N = 400
  δNp = N + 1

  # SBP interior order
  SBPp   = 2

  # mesh file side set type to actually boundary condition type
  bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
            BC_JUMP_INTERFACE]

  (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("../meshes/1_1_block.inp";
            bc_map=bc_map)

   # number of elements and faces
   (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
   @show (nelems, nfaces)


   # EToN0 is the base mesh size (e.g., before refinement)
   EToN0 = zeros(Int64, 2, nelems)
   EToN0[1, :] .= N
   EToN0[2, :] .= N

   # Exact solution
   Lx = 80
   Ly = 80


     # Set up the local grid dimensions
     Nr = EToN0[1, :]
     Ns = EToN0[2, :]

     # Dictionary to store the operators
     OPTYPE = typeof(locoperator(2, 8, 8))
     lop = Dict{Int64,OPTYPE}() # Be extra careful about the () here

    #  el_x = 10e12 # set it to be infinity to have even spread
    #  el_y = 10e12
    el_x = 10
    el_y = 10
    xt = (r,s) -> (el_x .* tan.(atan((Lx )/el_x).* (0.5*r .+ 0.5))  , el_x .* sec.(atan((Lx )/el_x).* (0.5*r .+ 0.5)).^2 * atan((Lx)/el_x) * 0.5 ,zeros(size(s)))
    yt = (r,s) -> (el_y .* tan.(atan((Ly )/el_y).* (0.5*s .+ 0.5))  , zeros(size(r)), el_y .* sec.(atan((Ly )/el_y).*(0.5*s .+ 0.5)) .^2 * atan((Ly )/el_y) * 0.5 )


     # create metrics
     metrics = create_metrics(SBPp, Nr[1], Ns[1], xt, yt) # not quite sure about this part

     # create local operator
     LFtoB = [BC_DIRICHLET,BC_DIRICHLET,BC_NEUMANN,BC_NEUMANN]
     lop = Dict{Int64,OPTYPE}() # This step to create a dict is essential
     lop[1] = locoperator(SBPp, Nr[1], Ns[1], metrics, LFtoB) # this function might not be correct

     # obtain M
     factorization = (x) -> cholesky(Symmetric(x))
     M = SBPLocalOperator1(lop, Nr[1], Ns[1], factorization)

     # obtain ge that stores boundary data
     ge = zeros((Nr[1]+1) * (Ns[1]+1))
     e = 1


     # boundary conditions
     bc_Dirichlet = (lf, x, y, e) -> zeros(size(x))
     bc_Neumann   = (lf, x, y, nx, ny, e) -> zeros(size(x))
     # set right-hand side and solve
     locbcarray_mod!(ge, lop[e], LFtoB, bc_Dirichlet, bc_Neumann,(e))
     u = M.F[e] \ ge

  Δτ = zeros(N+1)
  τ = zeros(N+1)

  # Assemble fault variables/data
  RSa = zeros(N+1)

  xf = lop[1].facecoord[1][1]
  yf = lop[1].facecoord[2][1]
  for n = 1:N+1
      RSa[n] = RSamin - (RSamin - RSamax) *
        min(1, max(0, (RSH1 - yf[n])/(RSH1 - RSH2)))
  end

  τz0 = σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                 exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
                                      RSamax)) + η * RSVinit


  θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τz0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)


  ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)


  ψδ = zeros(2δNp)
  ψδ[1:δNp] .= ψ0

  #δ = zeros(δNp)
  #pyplot()
  #display(plot(yf, δ))
  #sleep(1)
  # @show yf
  stations = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35]
  # @show stations
  # numstations = length(stations)
  # station_ind = zeros(Int64, numstations)
  @show size(yf)
  @show size(stations)

  function find_station_index(stations,grid_points)
    numstations = length(stations)
    @show numstations
    station_ind = zeros(numstations)
    @show station_ind
    for i in range(1,stop=numstations)
      @show argmin(abs.(grid_points .- stations[i]))
        station_ind[i] = argmin(abs.(grid_points .- stations[i]))
        # station_ind[i] = abs.(grid_points .- stations[i])[1]
    end
    return Integer.(station_ind)
  end  

  station_indices = find_station_index(stations,yf)
  @show station_indices
  # station_face = zeros(Int64, numstations)
  # for s = 1:numstations
  #   n = station_ind[s] = argmin(abs.(fault_y-stations[s]))
  #   f = station_face[s] = findlast((m) -> m <= n, FToδstarts)
  #   println((stations[s], fault_y[station_ind[s]], FToδstarts[f] <= n < FToδstarts[f+1]))
  # end
  # station_t = zeros(Float64, ceil(Int64, sim_years)*10)
  # station_V = zeros(Float64, numstations, ceil(Int64, sim_years)*10)
  # station_τ = zeros(Float64, numstations, ceil(Int64, sim_years)*10)
  # station_θ = zeros(Float64, numstations, ceil(Int64, sim_years)*10)
  # station_δ = zeros(Float64, numstations, ceil(Int64, sim_years)*10)

  # stations_locations = [0 0
  #                       0 -2.5
  #                       0 -5
  #                       0 -7.5
  #                       0 -10
  #                       0 -12.5
  #                       0 -15
  #                       0 -17.5
  #                       0 -20
  #                       0 -25
  #                       0 -30
  #                      ]
  # y_loc = stations_locations[:,2]
  # y_ind = 

  

  # set up parameters sent to right hand side
  odeparam = (reject_step = [false],
              Vp=Vp,
              lop=lop,
              F = M.F[e],
              u=u,
              Δτ = Δτ,
              τ = τ,
              ge = ge,
              μshear=μshear,
              RSa=RSa,
              RSb=RSb,
              σn=σn,
              η=η,
              RSV0=RSV0,
              τz0=τz0,
              # τn = τz0,
              RSDc=RSDc,
              RSf0=RSf0,
              LFtoB = LFtoB
             )

  dψV = zeros(2δNp)
  tspan = (0, sim_years * year_seconds)
  prob = ODEProblem(odefun, ψδ, tspan, odeparam)
  function stepcheck(_, p, _)
    if p.reject_step[1]
      p.reject_step[1] = false
      println("reject")
      return true
    end
    return false
  end

  ODEresults = ODE_results([],[],[],Dict(i=>[] for i=1:11))

  cb = SavingCallback((ψδ,t,i)->saveslip(ψδ,t,i,ODEresults,yf,stations,station_indices,odeparam,"BP1_",10*year_seconds),SavedValues(Float64,Float64))

  sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, dt=year_seconds,
              atol = 1e-5, rtol = 1e-3, save_everystep=true,
              internalnorm=(x, _)->norm(x, Inf),callback=cb)


  return (sol, yf)
end


function plot_slip(S, yf, stride_time)

  m = length(yf)
  no_time_steps = size(S.t)
  slip_final = S.u[end][end]

  for i = 1:stride_time:no_time_steps[1]

    slip_t = S.u[i][m+1:end] # slip at time t
    #pyplot()
    display(plot(slip_t, -yf, xtickfont=font(18),
    ytickfont=font(18),
    guidefont=font(18),
    legendfont=font(18), ylabel = "Depth (km)", xlabel = "Slip (m)", xlims = (0, slip_final)))
    sleep(0.1)
  end

  #nothing
end




(S, yf) = main()
# plot_slip(S, yf, 10)
