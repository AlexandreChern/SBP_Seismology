include("global_curved.jl");

SBPp   = 6

bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
          BC_JUMP_INTERFACE]
(verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/square_circle.inp";
        bc_map = bc_map)
verts

EToV

EToF

FToB

EToDomain

(nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
@show (nelems, nfaces)

for v in 1:size(verts, 2)
  x, y = verts[1, v], verts[2, v]
  if abs(hypot(x,y) - 1) < 1e-5
    Q = atan(y, x)
    verts[1, v], verts[2, v] = cos(Q), sin(Q)
  end
end

verts

plot_connectivity(verts, EToV)

N1 = N0 = 16

EToN0 = zeros(Int64,2,nelems)

EToN0[1,:] .= N0

@. EToN0[1,:] = N0
@. EToN0[2,:] = N1

@assert typeof(EToV) == Array{Int, 2} && size(EToV) == (4,nelems)
@assert typeof(EToF) == Array{Int, 2} && size(EToF) == (4,nelems)
@assert maximum(maximum(EToF)) == nfaces

(FToE, FToLF, EToO, EToS) = connectivityarrays(EToV,EToF)

Lx = maximum(verts[1,:])
Ly = maximum(abs.(verts[2,:]))

(kx, ky) = (2*π / Lx, 4*π / Ly)



# Here it defines vex function
vex(x,y,e) = begin
    if EToDomain[e] == 1
        return cos.(kx * x) .* cosh.(ky * y)
    elseif EToDomain[e] == 2
        return 10 .+ cos.(kx * x) .* cosh(ky * y)
    else
        error("invalid block")
    end
end



# vex_x and vex_y define first derivative of vex functions respectively

vex_x(x,y,e) =  begin
    if ETODomain[e] == 1
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
        return ky * cos.(kx * x) .* sinh.*(ky * y)
    else
        error("invalid block")
    end
end

# vex_xx vex_xy, vex_yy define the second order derivatives of the vex functions respectively




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




ϵ = zeros(4)

lvl = 1

Nr = EToN0[1, :] * (2^(lvl - 1))
Ns = EToN0[2, :] * (2^(lvl - 1))

OPTYPE = typeof(locoperator(2,8,8))
lop = Dict{Int64, OPTYPE}()

for e = 1:nelems
    (x1, x2, x3, x4) = verts[1,EToV[:,e]]
    (y1, y2, y3, y4) = verts[2,EToV[:,e]]
end

e = 1

verts[1,EToV[:,1]]

verts

(x1, x2, x3, x4) = verts[1,EToV[:,1]]
(y1, y2, y3, y4) = verts[2,EToV[:,1]]
ex = [ (α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
        (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
        (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
        (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]

ex[1](3)

exα = [
    (α) -> -x1 / 2 + x3 / 2,
    (α) -> -x2 / 2 + x4 / 2,
    (α) -> -x1 / 2 + x2 / 2,
    (α) -> -x3 / 2 + x4 / 2
]


ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
      (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
      (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
      (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]

eyα = [(α) -> -y1 / 2 + y3 / 2,
       (α) -> -y2 / 2 + y4 / 2,
       (α) -> -y1 / 2 + y2 / 2,
       (α) -> -y3 / 2 + y4 / 2]

if FToB[EToF[1, e]] == BC_JUMP_INTERFACE
 error("curved face 1 not implemented yet")
end
if FToB[EToF[2, e]] == BC_JUMP_INTERFACE
 error("curved face 2 not implemented yet")
end

if FToB[EToF[3,e]] == BC_JUMP_INTERFACE
    Q1 = atan(y1,x1)
    Q2 = atan(y2,x2)
    if !(-π/2 < Q1 - Q2 < π/2)
        Q2 -= sign(Q2) * 2 * π
    end
    ex[3] = (α) -> cos.(Q1 * (1 .- α)/2 + Q2 * (1 .+ α) / 2)
    ey[3] = (α) -> sin.(Q1 * (1 .- α)/2 + Q2 * (1 .+ α) / 2)
    β3 = (Q2 - Q1)/2
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

xt(r,s) = transfinite_blend(ex[1], ex[2], ex[3], ex[4],
                            exα[1], exα[2], exα[3], exα[4],
                            r, s)

yt(r,s) = transfinite_blend(ey[1], ey[2], ey[3], ey[4],
                            eyα[1], eyα[2], eyα[3], eyα[4],
                            r, s)


metrics = create_metrics(SBPp, Nr[e], Ns[e], xt, yt)

lop_test = locoperator(SBPp, Nr[e], Ns[e], metrics, FToB[EToF[:,e]])


lop[e] = locoperator(SBPp, Nr[e], Ns[e], metrics, FToB[EToF[:,e]]);






lop_test
