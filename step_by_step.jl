include("global_curved.jl");

SBPp   = 6

# mesh file side set type to actually boundary condition type
bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
          BC_JUMP_INTERFACE]
(verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/square_circle.inp";
                                                   bc_map = bc_map)
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

ϵ = zeros(4)

length(ϵ)

lvl = 1

# Set up the local grid dimensions
Nr = EToN0[1, :] * (2^(lvl-1))
Ns = EToN0[2, :] * (2^(lvl-1))

#
# Build the local volume operators
#

# Dictionary to store the operators
OPTYPE = typeof(locoperator(2, 8, 8))
lop = Dict{Int64, OPTYPE}() # Be extra careful about the () here

nelems

e = 1

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
