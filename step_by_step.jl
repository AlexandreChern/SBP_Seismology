include("global_curved.jl");

SBPp   = 6

bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
          BC_JUMP_INTERFACE]
(verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/square_circle.inp";
        bc_map = bc_map)
