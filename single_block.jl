include("global_curved.jl")
let

    # setting parameters for FD Cycle
    sim_years = 10000

    Vp = 1e-9
    ρ = 2.670
    cs = 3.464
    σn = 50
    RSamin = 0.010
    RSamax = 0.025
    RSb = 0.015
    RSDc = 0.008
    RSf0 = 0.6
    RSV0 = 1e-6
    RSH1 = 15
    RSH2 = 18

    μshear = cs^2 * ρ
    η = μshear / (2 * cs)



    # SBP interior order
    SBPp   = 6

    # number of levels
    num_of_levels = 6 # no mesh refinement right now

    # mesh file side set type to actually boundary condition type
    bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
            BC_JUMP_INTERFACE]
    
    (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/1_1_block.inp";
            bc_map=bc_map)
    
    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    @show (nelems, nfaces)

    plot_connectivity(verts, EToV)
    # This is the base mesh size in each dimension
    N1 = N0 = 16

     # EToN0 is the base mesh size (e.g., before refinement)
    EToN0 = zeros(Int64, 2, nelems)
    EToN0[1, :] .= N0
    EToN0[2, :] .= N1

    @assert typeof(EToV) == Array{Int,2} && size(EToV) == (4, nelems)
    @assert typeof(EToF) == Array{Int,2} && size(EToF) == (4, nelems)
    @assert maximum(maximum(EToF)) == nfaces


    # Exact solution
    # Lx = maximum(verts[1,:])
    # Ly = maximum(abs.(verts[2,:]))
    Lx = 80
    Ly = 80
    (kx, ky) = (2 * π / Lx, 4 * π / Ly)
    vex(x, y, e) = begin # I only have one element. no need for EToDomain[e]
        if EToDomain[e] == 1
            # return cos.(kx * x) .* cosh.(ky * y)
            # return x
            return cos.(kx * x) .* cos.(ky * y)
        # elseif EToDomain[e] == 2
        #     return 10 .+ cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_x(x, y, e) = begin
        if EToDomain[e] == 1
            # return -kx * sin.(kx * x) .* cosh.(ky * y)
            # return ones(size(x)[1])
            return - kx * sin.(kx * x) .* cos.(ky * y)
        # elseif EToDomain[e] == 2
        #     return -kx * sin.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_y(x, y, e) = begin
        if EToDomain[e] == 1
            # return ky * cos.(kx * x) .* sinh.(ky * y)
            # return zeros(size(x)[1])
            return - ky * cos.(kx * x) .* sin.(ky * y)
        # elseif EToDomain[e] == 2
        #     return ky * cos.(kx * x) .* sinh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_xx(x, y, e) = begin
        if EToDomain[e] == 1
            # return -kx^2 * cos.(kx * x) .* cosh.(ky * y)
            # return zeros(size(x)[1])
            return - kx^2 * cos.(kx * x) .* cos.(ky * y)
        # elseif EToDomain[e] == 2
        #     return -kx^2 * cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_xy(x, y, e) = begin
        if EToDomain[e] == 1
            # return -kx * ky * sin.(kx * x) .* sinh.(ky * y)
            # return zeros(size(x)[1])
            return kx * ky * sin.(kx * x) .* sin.(ky * y)
        # elseif EToDomain[e] == 2
        #     return -kx * ky * sin.(kx * x) .* sinh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_yy(x, y, e) = begin
        if EToDomain[e] == 1
            # return ky^2 * cos.(kx * x) .* cosh.(ky * y)
            # return zeros(size(x)[1])
            return - ky^2 * cos.(kx * x) .* cos.(ky * y)
        # elseif EToDomain[e] == 2
        #     return ky^2 * cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end


    ϵ = zeros(num_of_levels)
    @show verts
    for lvl = 1:length(ϵ)
        # Set up the local grid dimensions
        Nr = EToN0[1, :] * (2^(lvl - 1))
        Ns = EToN0[2, :] * (2^(lvl - 1))

        @show Nr
        @show Ns

        # Dictionary to store the operators
        OPTYPE = typeof(locoperator(2, 8, 8))
        lop = Dict{Int64,OPTYPE}() # Be extra careful about the () here

        (x1, x2, x3, x4) = verts[1,:]
        (y1, y2, y3, y4) = verts[2,:]

        # xt = ?
        # function xt(r,s)
        #     return r
        # end
        # xt = (r,s) -> (0.5 * r.^3 .+ r, 1.5 .* r .* r .+ 1, zeros(size(r)))
        el_x = 10 # set it to be infinity to have even spread
        el_y = 10
        # xt = (r,s) -> (el_x .* tan.(atan((Lx + 1)/el_x).* r)  .- 1, el_x .* sec.(atan((Lx + 1)/el_x).* r).^2 * atan((Lx + 1)/el_x) ,zeros(size(s)))
        xt = (r,s) -> (el_x .* tan.(atan((Lx )/el_x).* (0.5*r .+ 0.5))  , el_x .* sec.(atan((Lx )/el_x).* (0.5*r .+ 0.5)).^2 * atan((Lx)/el_x) * 0.5 ,zeros(size(s)))
        # yt = (r,s) -> (el_y .* tan.(atan((Ly + 1)/el_y).* s)  .- 1 , zeros(size(r)), el_y .* sec.(atan((Ly + 1)/el_y).*s) .^2 * atan((Ly + 1)/el_y) )
        yt = (r,s) -> (el_y .* tan.(atan((Ly )/el_y).* (0.5*s .+ 0.5))  , zeros(size(r)), el_y .* sec.(atan((Ly )/el_y).*(0.5*s .+ 0.5)) .^2 * atan((Ly )/el_y) * 0.5 )

        # xt = (r,s) -> (r, ones(size(r)), zeros(size(r))) # maybe the mapping is wrong
        # yt = ?
        # function yt(r,s)
        #     return s
        # end
        # yt = (r,s) -> (0.5 * s.^3 .+ s,zeros(size(s)), 1.5 .* s .* s .+ 1)
        
        # create metrics
        metrics = create_metrics(SBPp, Nr[1], Ns[1], xt, yt) # not quite sure about this part
        # println("Metrics created successfully")

        # create local operator
        LFtoB = [BC_DIRICHLET,BC_DIRICHLET,BC_NEUMANN,BC_NEUMANN]
        lop = Dict{Int64,OPTYPE}() # This step to create a dict is essential
        lop[1] = locoperator(SBPp, Nr[1], Ns[1], metrics, LFtoB) # this function might not be correct
        # for locoperator() call with more elements, the length of lop is the number of the elements 
        # for locoperator() call with a single element, the length of lop is the length of the tuple. (how many items for a single elemet)
        # println("lop created successfully")
        # @show typeof(lop)
        @show length(lop)

        # obtain M 
        factorization = (x) -> cholesky(Symmetric(x))
        M = SBPLocalOperator1(lop, Nr[1], Ns[1], factorization)
        # println("M created successfully")
        # M = locoperator(SBPp, Nr[1],Ns[1],metrics,LFtoB)

        # obtain g
        # source function
        ge = zeros((Nr[1]+1) * (Ns[1]+1))
        gδe = zeros((Nr[1]+1) * (Ns[1] + 1))
        e = 1
        δ = 1
        source = (x, y, e) -> (-vex_xx(x, y, e)  - vex_yy(x, y, e))
        # locbcarray to implement

        # boundary conditions
        bc_Dirichlet = (lf, x, y, e, δ) -> vex(x, y, e)
        bc_Neumann   = (lf, x, y, nx, ny, e, δ) -> (nx .* vex_x(x, y, e)
                                                + ny .* vex_y(x, y, e))
        in_jump      = (lf, x, y, e, δ) -> begin
            f = EToF[lf, e]
            if EToS[lf, e] == 1
                if EToO[lf, e]
                    return -δ[FToδstarts[f]:(FToδstarts[f + 1] - 1)]
                else
                    error("shouldn't get here")
                end
            else
                if EToO[lf, e]
                    return  δ[FToδstarts[f]:(FToδstarts[f + 1] - 1)]
                else
                    return  δ[(FToδstarts[f + 1] - 1):-1:FToδstarts[f]]
                end
            end
        end

        locbcarray!(ge, gδe, lop[e], LFtoB, bc_Dirichlet, bc_Neumann, in_jump,(e,δ))
        locsourcearray!((ge),source,lop[e],e)
        # @show g
        # @show M
        numerical_solution = M.F[e] \ ge

        odeparams = ( reject_step = [false],
            Vp = Vp,
            lop = lop,
            EToDomain = EToDomain,
            EToF = EToF,
            EToV = EToV,
            EToN0 = EToN0,
            FToB = FToB,
            FToδstarts = FToδstarts
        )

        # exact_solution = vex(lop[e].coord[1], lop[e].coord[2],e)
        # diff = numerical_solution - exact_solution[:]
        # diff_matrix = reshape(numerical_solution,Ns[1]+1,Nr[1]+1)
        # err = sqrt(diff' * lop[1].JH * diff)
        # println("Err: $err")
        # println("LogErr: $(log(2,err))")
    end
end
