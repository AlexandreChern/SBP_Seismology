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
