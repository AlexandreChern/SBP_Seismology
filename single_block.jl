include("global_curved.jl")
let
    # SBP interior order
    SBPp   = 6

    # number of levels
    num_of_levels = 1 # no mesh refinement right now

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
    Lx = maximum(verts[1,:])
    Ly = maximum(abs.(verts[2,:]))
    (kx, ky) = (2 * π / Lx, 4 * π / Ly)
    vex(x, y, e) = begin
        if EToDomain[e] == 1
            return cos.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return 10 .+ cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_x(x, y, e) = begin
        if EToDomain[e] == 1
            return -kx * sin.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return -kx * sin.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_y(x, y, e) = begin
        if EToDomain[e] == 1
            return ky * cos.(kx * x) .* sinh.(ky * y)
        elseif EToDomain[e] == 2
            return ky * cos.(kx * x) .* sinh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_xx(x, y, e) = begin
        if EToDomain[e] == 1
            return -kx^2 * cos.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return -kx^2 * cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_xy(x, y, e) = begin
        if EToDomain[e] == 1
            return -kx * ky * sin.(kx * x) .* sinh.(ky * y)
        elseif EToDomain[e] == 2
            return -kx * ky * sin.(kx * x) .* sinh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_yy(x, y, e) = begin
        if EToDomain[e] == 1
            return ky^2 * cos.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return ky^2 * cos.(kx * x) .* cosh.(ky * y)
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

        # Dictionary to store the operators
        OPTYPE = typeof(locoperator(2, 8, 8))
        lop = Dict{Int64,OPTYPE}() # Be extra careful about the () here

        (x1, x2, x3, x4) = verts[1,:]
        (y1, y2, y3, y4) = verts[2,:]

        # xt = ?
        # xt(r,s) = r
        # yt = ?
        # yt(r,s) = s
        
        # create metrics
        metrics = create_metrics(SBPp, Nr, Ns, xt, yt) # not quite sure about this part

        # create local operator
        LFtoB = [BC_DIRICHLET,BC_DIRICHLET,BC_NEUMANN,BC_NEUMANN]
        lop = locoperator(SBPp, Nr, Ns, metrics, LFtoB) # this function might not be correct 

        # obtain M 
        factorization = (x) -> cholesky(Symmetric(x))
        M = SBPLocalOperator1(lop, Nr, Ns, factorization)

        # obtain g





    end
end
