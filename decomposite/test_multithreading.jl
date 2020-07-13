using .Threads
using Dates
using BenchmarkTools
using ArgParse
include("global_curved.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s  begin
        "--block_num", "-b"
            help = "an option with an argument"
            arg_type = Int
            default = 4
        "--level_num", "-l"
            help = "another option with an argument"
            arg_type = Int
            default = 4

        "--sbp_level", "-o"
            help = "SBP operators order"
            arg_type = Int
            default = 6
        # "--flag1"
        #     help = "an option without argument, i,e, a flag"
        #     action = :store_true
        # "arg1"
        #     help = "a positional argument"
        #     required = true
    end
    return parse_args(s)
end

# function main()
#     parsed_args = parse_commandline()
#     println("Parsed args:")
#     for (arg,val) in parsed_args
#         println(" $arg => $val")
#     end
# end

parsed_args = parse_commandline()

# @show parsed_args
#
block_num = parsed_args["block_num"]
level_num = parsed_args["level_num"]
SBP_lvl = parsed_args["sbp_level"]

# @show block_num
# @show level_num


let
    # number of blocks in each side
    # n_block = 8
    n_block = block_num
    # n_block = $block_
    # SBP interior order
    # SBPp   = 6
    SBPp = SBP_lvl
    # num_of_lvls = 4
    num_of_lvls = level_num
    @show n_block
    @show num_of_lvls
    @show SBPp
    current_time = now()
    string_time =  String((Symbol("_",Dates.month(current_time),'_',Dates.day(current_time),'_',Dates.hour(current_time),'_',Dates.minute(current_time))))
    input_file_name =  String((Symbol(n_block,"_",n_block,"_block.inp")))
    output_file_name =  String((Symbol(n_block,"_",n_block,"_p_",SBPp,string_time,"_output.txt")))
    fileio = open("results/" * output_file_name,"w")


    # mesh file side set type to actually boundary condition type
    bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,BC_JUMP_INTERFACE]
    (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("../meshes/" * input_file_name;bc_map = bc_map)
    # EToV defines the element by its vertices
    # EToF defines element by its four faces, in global face number
    # FToB defines whether face is Dirichlet (1), Neumann (2), interior jump (7)
    #      or just an interior interface (0)
    # EToDomain is 1 if element is inside circle; 2 otherwise

    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    @show (nelems, nfaces)

    # This is needed to fix up points that should be on the boundary of the
    # circle, but the mesh didn't quite put them there
    for v in 1:size(verts, 2)
        x, y = verts[1, v], verts[2, v]
        if abs(hypot(x,y) - 1) < 1e-5
            Q = atan(y, x)
            verts[1, v], verts[2, v] = cos(Q), sin(Q)
        end
    end


    # Plot the original connectivity before mesh warping
    plot_connectivity(verts, EToV)

    # This is the base mesh size in each dimension
    N1 = N0 = 16

    # EToN0 is the base mesh size (e.g., before refinement)
    EToN0 = zeros(Int64, 2, nelems)
    EToN0[1, :] .= N0
    EToN0[2, :] .= N1

    @assert typeof(EToV) == Array{Int, 2} && size(EToV) == (4, nelems)
    @assert typeof(EToF) == Array{Int, 2} && size(EToF) == (4, nelems)
    @assert maximum(maximum(EToF)) == nfaces

    # Determine secondary arrays
    # FToE : Unique Global Face to Element Number
    #        (the i'th column of this stores the element numbers that share the
    #        global face number i)
    # FToLF: Unique Global Face to Element local face number
    #        (the i'th column of this stores the element local face numbers that
    #        shares the global face number i)
    # EToO : Element to Unique Global Faces Orientation
    #        (the i'th column of this stores the whether the element and global
    #        face are oriented in the same way in physical memory or need to be
    #        rotated)
    # EToS : Element to Unique Global Face Side
    #        (the i'th column of this stores whether an element face is on the
    #        plus side or minus side of the global face)
    (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

    # Exact solution
    Lx = maximum(verts[1,:])
    Ly = maximum(abs.(verts[2,:]))
    (kx, ky) = (2*π / Lx, 4*π / Ly)
    vex(x,y,e) = begin
        if EToDomain[e] == 1
            return cos.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return 10 .+ cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_x(x,y,e) = begin
        if EToDomain[e] == 1
            return -kx * sin.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return -kx * sin.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_y(x,y,e) = begin
        if EToDomain[e] == 1
            return ky * cos.(kx * x) .* sinh.(ky * y)
        elseif EToDomain[e] == 2
            return ky * cos.(kx * x) .* sinh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_xx(x,y,e) = begin
        if EToDomain[e] == 1
            return -kx^2 * cos.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return -kx^2 * cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_xy(x,y,e) = begin
        if EToDomain[e] == 1
            return -kx * ky * sin.(kx * x) .* sinh.(ky * y)
        elseif EToDomain[e] == 2
            return -kx * ky * sin.(kx * x) .* sinh.(ky * y)
        else
            error("invalid block")
        end
    end
    vex_yy(x,y,e) = begin
        if EToDomain[e] == 1
            return ky^2 * cos.(kx * x) .* cosh.(ky * y)
        elseif EToDomain[e] == 2
            return ky^2 * cos.(kx * x) .* cosh.(ky * y)
        else
            error("invalid block")
        end
    end

    ϵ = zeros(num_of_lvls)
    ϵ_test = zeros(num_of_lvls) # this array store errors calculated from multiple run-time which is wrong and very big
    for lvl = 1:length(ϵ)
        start = time()
        # Set up the local grid dimensions
        Nr = EToN0[1, :] * (2^(lvl-1))
        Ns = EToN0[2, :] * (2^(lvl-1))

        #
        # Build the local volume operators
        #

        # Dictionary to store the operators
        OPTYPE = typeof(locoperator(2, 8, 8))
        lop = Dict{Int64, OPTYPE}()

        # Loop over blocks and create local operators
        @threads for e = 1:nelems
            # Get the element corners
            (x1, x2, x3, x4) = verts[1, EToV[:, e]]
            (y1, y2, y3, y4) = verts[2, EToV[:, e]]

            # Initialize the block transformations as transfinite between the corners
            ex = [(α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
            (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
            (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
            (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]
            exα = [(α) -> -x1 / 2 + x3 / 2,
            (α) -> -x2 / 2 + x4 / 2,
            (α) -> -x1 / 2 + x2 / 2,
            (α) -> -x3 / 2 + x4 / 2]
            ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
            (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
            (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
            (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]
            eyα = [(α) -> -y1 / 2 + y3 / 2,
            (α) -> -y2 / 2 + y4 / 2,
            (α) -> -y1 / 2 + y2 / 2,
            (α) -> -y3 / 2 + y4 / 2]

            # For blocks on the circle, put in the curved edge transform
            if FToB[EToF[1, e]] == BC_JUMP_INTERFACE
                error("curved face 1 not implemented yet")
            end
            if FToB[EToF[2, e]] == BC_JUMP_INTERFACE
                error("curved face 2 not implemented yet")
            end
            if FToB[EToF[3, e]] == BC_JUMP_INTERFACE
                Q1 = atan(y1, x1)
                Q2 = atan(y2, x2)
                if !(-π/2 < Q1 - Q2 < π/2)
                    Q2 -= sign(Q2) * 2 * π
                end
                ex[3] = (α) -> cos.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
                ey[3] = (α) -> sin.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
                β3 = (Q2 - Q1) / 2
                exα[3] = (α) -> -β3 .* sin.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
                eyα[3] = (α) -> +β3 .* cos.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
            end
            if FToB[EToF[4, e]] == BC_JUMP_INTERFACE
                Q3 = atan(y3, x3)
                Q4 = atan(y4, x4)
                if !(-π/2 < Q3 - Q4 < π/2)
                    error("curved face 4 angle correction not implemented yet")
                end
                ex[4] = (α) -> cos.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
                ey[4] = (α) -> sin.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
                β4 = (Q4 - Q3) / 2
                exα[4] = (α) -> -β4 .* sin.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
                eyα[4] = (α) -> +β4 .* cos.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
            end

            # Create the volume transform as the transfinite blending of the edge
            # transformations
            xt(r,s) = transfinite_blend(ex[1], ex[2], ex[3], ex[4],exα[1], exα[2], exα[3], exα[4],r, s)
            yt(r,s) = transfinite_blend(ey[1], ey[2], ey[3], ey[4],eyα[1], eyα[2], eyα[3], eyα[4],r, s)


            metrics = create_metrics(SBPp, Nr[e], Ns[e], xt, yt)

            # Build local operators
            lop[e] = locoperator(SBPp, Nr[e], Ns[e], metrics, FToB[EToF[:, e]])
        end
        m_list = [lop[e].M̃ for e in 1:nelems]
        println(length(m_list))
        println(length(unique(m_list)))


        # If this is the first mesh level plot the mesh
        # lvl == 1 && plot_blocks(lop)

        #
        # Do some assemble of the global volume operators
        #
        # (M, FbarT, D, vstarts, FToλstarts) =
        # LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
        #                    (x) -> cholesky(Symmetric(x)))

        M = SBPLocalOperator1_forming(lop,Nr,Ns)
        chol_fac = (x) -> cholesky(Symmetric(x))
        lu_fac = (x) -> lu(x)
        fac_method = chol_fac
        # factors = SBPLocalOperator1_factorization(lop,Nr,Ns,(x) -> cholesky(Symmetric(x)))
        factors = SBPLocalOperator1_factorization(lop,Nr,Ns,fac_method)

        vstarts = M.offset
        (FToλstarts, FbarT,D) = threaded_gloλoperator(lop,M.offset,FToB,FToE,FToLF,EToO,EToS,Nr,Ns)

        println()

        @show lvl
        write(fileio,"\n")
        write(fileio,"lvl = $lvl \n")

        println("N: $(2^(lvl+3)*n_block)")
        write(fileio, "N: $(2^(lvl+3)*n_block)\n")

        # locfactors = M.F
        # @show length(locfactors)
        locfactors = factors
        # @show length(factors)

        # Get a unique array indexes for the face to jumps map
        FToδstarts = bcstarts(FToB, FToE, FToLF, BC_JUMP_INTERFACE, Nr, Ns)

        # Compute the number of volume, trace (λ), and jump (δ) points
        VNp = vstarts[nelems+1]-1
        λNp = FToλstarts[nfaces+1]-1
        δNp = FToδstarts[nfaces+1]-1
        # @show (VNp, λNp, δNp)

        # Build the (sparse) λ matrix using the schur complement and factor
        start_assembleλmatrix = time()
        (B,elapsed) = assembleλmatrix_test(FToλstarts, vstarts, EToF, FToB, locfactors, D, FbarT)
        println("Time for direct solve in forming λ: $elapsed")
        write(fileio,"Time for direct solve in forming λ: $elapsed\n")

        elapsed_assembleλmatrix = time() - start_assembleλmatrix
        # Test code to show the size of each input variabls
        # @show Base.summarysize(FToλstarts)
        # @show Base.summarysize(vstarts)
        # @show Base.summarysize(EToF)
        # @show Base.summarysize(FToB)
        # @show Base.summarysize(locfactors)
        # @show Base.summarysize(D)
        # @show Base.summarysize(FbarT)
        # @show Base.summarysize(B)
        println(length(locfactors))
        println(length(unique(locfactors)))
        println("Time elapsed (assembleλmatrix) for lvl $lvl = $elapsed_assembleλmatrix")
        write(fileio,"Time elapsed (assembleλmatrix) for lvl $lvl = $elapsed_assembleλmatrix\n")
        #
        # rs1 = @benchmark assembleλmatrix_test($FToλstarts, $vstarts, $EToF, $FToB, $locfactors, $D, $FbarT)

        # println(Base.summarysize(B))
        # display(rs1)
        # 
        # rs2 = @benchmark assembleλmatrix($FToλstarts, $vstarts, $EToF, $FToB, $locfactors, $D, $FbarT)
        #
        # display(rs2)

        BF = cholesky(Symmetric(B))

        (bλ, λ, gδ) = (zeros(λNp), zeros(λNp), zeros(λNp))
        (Δ, u, g) = (zeros(VNp), zeros(VNp), zeros(VNp))
        δ = zeros(δNp)
        for f = 1:nfaces
            if FToB[f] == BC_JUMP_INTERFACE
                (e1, e2) = FToE[:, f]
                (lf1, lf2) = FToLF[:, f]
                (xf, yf) = lop[e1].facecoord
                @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
                vex(xf[lf1], yf[lf1], e2) - vex(xf[lf1], yf[lf1], e1)
            end
        end

        bc_Dirichlet = (lf, x, y, e, δ) -> vex(x, y, e)
        bc_Neumann = (lf, x, y, nx, ny, e, δ) -> (nx .* vex_x(x, y, e) + ny .* vex_y(x, y, e))
        in_jump  = (lf, x, y, e, δ) -> begin
        f = EToF[lf, e]
        if EToS[lf, e] == 1
            if EToO[lf, e]
                return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
            else
                error("shouldn't get here")
            end
        else
            if EToO[lf, e]
                return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
            else
                return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
            end
        end
    end

    @threads for e = 1:nelems
        gδe = ntuple(4) do lf
            f = EToF[lf, e]
            if EToO[lf, e]
                return @view gδ[FToλstarts[f]:(FToλstarts[f+1]-1)]
            else
                return  @view gδ[(FToλstarts[f+1]-1):-1:FToλstarts[f]]
            end
        end
        locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), gδe, lop[e],
        FToB[EToF[:,e]], bc_Dirichlet, bc_Neumann, in_jump, (e, δ))

        source = (x, y, e) -> (-vex_xx(x, y, e)  - vex_yy(x, y, e))
        locsourcearray!((@view g[vstarts[e]:vstarts[e+1]-1]), source, lop[e], e)
    end

    LocalToGLobalRHS!(bλ, g, gδ,  u, locfactors, FbarT, vstarts)
    start_gs = time()
    λ[:] = BF \ bλ
    elapsed_gs = time() - start_gs
    println("Time elapsed for the global solve is approximately $elapsed_gs")
    write(fileio,"Time elapsed for the whole code is approximately $elapsed_gs \n")

    u[:] = -FbarT' * λ
    u[:] .= g .+ u

    @threads for e = 1:nelems
        F = locfactors[e]
        (x, y) = lop[e].coord
        JH = lop[e].JH

        @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
        #=
        ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), F,
        (@view u[vstarts[e]:(vstarts[e+1]-1)]))
        =#

        @views Δ[vstarts[e]:(vstarts[e+1]-1)] = (u[vstarts[e]:(vstarts[e+1]-1)] -
        vex(x[:], y[:], e))
        ϵ[lvl] += Δ[vstarts[e]:(vstarts[e+1]-1)]' * JH * Δ[vstarts[e]:(vstarts[e+1]-1)]
    end
    ϵ[lvl] = sqrt(ϵ[lvl])
    @show (lvl, ϵ[lvl])
    elapsed = time() - start


    # starting timing with repeating times
    repeat_times = 10
    elapsed1 = elapsed2 = elapsed3 = 0.0

    for n = 1:repeat_times
        start1 = time()
        @threads for e=1:nelems
            F = locfactors[e]
            (x, y) = lop[e].coord
            JH = lop[e].JH
        end
        elapsed1 += time() - start1


        start2 = time()
        @threads for e = 1:nelems
            F = locfactors[e]
            (x, y) = lop[e].coord
            JH = lop[e].JH

            @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
            #=
            ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), F,
            (@view u[vstarts[e]:(vstarts[e+1]-1)]))
            =#
        end
        elapsed2 += time() - start2


        start3 = time()
        @sync @threads for e = 1:nelems
            F = locfactors[e]
            (x, y) = lop[e].coord
            JH = lop[e].JH
            @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
            #=
            ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), F,
            (@view u[vstarts[e]:(vstarts[e+1]-1)]))
            =#

            @views Δ[vstarts[e]:(vstarts[e+1]-1)] = (u[vstarts[e]:(vstarts[e+1]-1)] -
            vex(x[:], y[:], e))
            ϵ_test[lvl] += Δ[vstarts[e]:(vstarts[e+1]-1)]' * JH * Δ[vstarts[e]:(vstarts[e+1]-1)]
        end
        elapsed3 += time() - start3
    end
    println("Time elapsed for the whole code is approximately $elapsed")
    write(fileio,"Time elapsed for the whole code is approximately $elapsed \n")

    println("Time elapsed (reading matrices) for lvl $lvl = $(elapsed1/repeat_times)")
    write(fileio,"Time elapsed (reading matrices) for lvl $lvl = $(elapsed1/repeat_times)\n")

    println("Time elapsed (linear solve with reading matrices) for lvl $lvl = $(elapsed2/repeat_times)")
    write(fileio,"Time elapsed (linear Solve with reading matrices) for lvl $lvl = $(elapsed2/repeat_times)\n")

    println("Time elapsed (All three parts) for lvl $lvl = $(elapsed3/repeat_times)")
    write(fileio,"Time elapsed (All three parts) for lvl $lvl = $(elapsed3/repeat_times)\n")
    write(fileio,string((lvl,ϵ[lvl])) * "\n")

end
println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))
write(fileio,string((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2)))
nothing
close(fileio)
end
