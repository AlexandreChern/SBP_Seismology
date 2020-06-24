using SparseArrays
using LinearAlgebra
using .Threads

if (typeof(Base.find_package("UnicodePlots")) == Nothing)
    println("Package Unicode not installed")
    using Pkg;
    Pkg.add("UnicodePlots")
end
using UnicodePlots

include("diagonal_sbp.jl")

# flatten tuples to arrays
flatten_tuples(x) = reshape(collect(Iterators.flatten(x)), length(x[1]),
                            length(x))

⊗(A,B) = kron(A, B)     # Defines the syntax meaning of ⊗ as Kronecker Product
                        # Equivalent form can be ⊗ = kron, this will also work

const BC_DIRICHLET        = 1
const BC_NEUMANN          = 2
const BC_LOCKED_INTERFACE = 0
const BC_JUMP_INTERFACE   = 7

#{{{ Transfinite Blend
function transfinite_blend(α1, α2, α3, α4, α1s, α2s, α3r, α4r, r, s)
  # +---4---+
  # |       |
  # 1       2
  # |       |
  # +---3---+
  @assert [α1(-1) α2(-1) α1( 1) α2( 1)] ≈ [α3(-1) α3( 1) α4(-1) α4( 1)]


  x = (1 .+ r) .* α2(s)/2 + (1 .- r) .* α1(s)/2 +
      (1 .+ s) .* α4(r)/2 + (1 .- s) .* α3(r)/2 -
     ((1 .+ r) .* (1 .+ s) .* α2( 1) +
      (1 .- r) .* (1 .+ s) .* α1( 1) +
      (1 .+ r) .* (1 .- s) .* α2(-1) +
      (1 .- r) .* (1 .- s) .* α1(-1)) / 4

  xr =  α2(s)/2 - α1(s)/2 +
        (1 .+ s) .* α4r(r)/2 + (1 .- s) .* α3r(r)/2 -
      (+(1 .+ s) .* α2( 1) +
       -(1 .+ s) .* α1( 1) +
       +(1 .- s) .* α2(-1) +
       -(1 .- s) .* α1(-1)) / 4


  xs = (1 .+ r) .* α2s(s)/2 + (1 .- r) .* α1s(s)/2 +
       α4(r)/2 - α3(r)/2 -
      (+(1 .+ r) .* α2( 1) +
       +(1 .- r) .* α1( 1) +
       -(1 .+ r) .* α2(-1) +
       -(1 .- r) .* α1(-1)) / 4

  return (x, xr, xs)
end

function transfinite_blend(α1, α2, α3, α4, r, s, p)
  (Nrp, Nsp) = size(r)
  (Dr, _, _, _) = diagonal_sbp_D1(p, Nrp-1; xc = (-1,1))
  (Ds, _, _, _) = diagonal_sbp_D1(p, Nsp-1; xc = (-1,1))

  α2s(s) = α2(s) * Ds'
  α1s(s) = α1(s) * Ds'
  α4r(s) = Dr * α4(r)
  α3r(s) = Dr * α3(r)

  transfinite_blend(α1, α2, α3, α4, α1s, α2s, α3r, α4r, r, s)
end

function transfinite_blend(v1::T1, v2::T2, v3::T3, v4::T4, r, s
                          ) where {T1 <: Number, T2 <: Number,
                                   T3 <: Number, T4 <: Number}
  e1(α) = v1 * (1 .- α) / 2 + v3 * (1 .+ α) / 2
  e2(α) = v2 * (1 .- α) / 2 + v4 * (1 .+ α) / 2
  e3(α) = v1 * (1 .- α) / 2 + v2 * (1 .+ α) / 2
  e4(α) = v3 * (1 .- α) / 2 + v4 * (1 .+ α) / 2
  e1α(α) = -v1 / 2 + v3 / 2
  e2α(α) = -v2 / 2 + v4 / 2
  e3α(α) = -v1 / 2 + v2 / 2
  e4α(α) = -v3 / 2 + v4 / 2
  transfinite_blend(e1, e2, e3, e4, e1α, e2α, e3α, e4α, r, s)
end
#}}}

#{{{ connectivityarrays
function connectivityarrays(EToV, EToF)
  # number of elements
  nelems = size(EToV, 2)
  nfaces = maximum(maximum(EToF))

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side

  FToE  = zeros(Int64, 2, nfaces)
  FToLF = zeros(Int64, 2, nfaces)
  EToO  = Array{Bool,2}(undef, 4, nelems)
  EToS  = zeros(Int64, 4, nelems)

  # Local Face to Local Vertex map
  LFToLV = flatten_tuples(((1,3), (2, 4), (1,2), (3,4)))
  for e = 1:nelems
    for lf = 1:4
      gf = EToF[lf, e]
      if FToE[1, gf] == 0
        @assert FToLF[1, gf] == 0
        FToE[1, gf] = e
        FToLF[1, gf] = lf
        EToO[lf, e] = true
        EToS[lf, e] = 1
      else
        @assert FToE[2, gf] == 0
        @assert FToLF[2, gf] == 0
        FToE[2, gf] = e
        FToLF[2, gf] = lf
        EToS[lf, e] = 2

        ne = FToE[1, gf]
        nf = FToLF[1, gf]

        nv = EToV[LFToLV[:,nf], ne]
        lv = EToV[LFToLV[:,lf], e]
        if nv == lv
          EToO[lf, e] = true
        elseif nv[end:-1:1] == lv
          EToO[lf, e] = false
        else
          error("problem with connectivity")
        end
      end
    end
  end
  (FToE, FToLF, EToO, EToS)
end
#}}}

#{{{ locoperator
function create_metrics(pm, Nr, Ns,
                        xf=(r,s)->(r, ones(size(r)), zeros(size(r))),
                        yf=(r,s)->(s, zeros(size(s)), ones(size(s))))
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  # Derivative operators for the metric terms
  @assert pm <= 8
  pp = pm == 6 ? 8 : pm

  r = range(-1, stop=1, length=Nrp)
  s = range(-1, stop=1, length=Nsp)

  # Create the mesh
  r = ones(1, Nsp) ⊗ r
  s = s' ⊗ ones(Nrp)
  (x, xr, xs) = xf(r, s)
  (y, yr, ys) = yf(r, s)

  J = xr .* ys - xs .* yr
  @assert minimum(J) > 0

  rx =  ys ./ J
  sx = -yr ./ J
  ry = -xs ./ J
  sy =  xr ./ J

  # variable coefficient matrix components
  crr = J .* (rx .* rx + ry .* ry)
  crs = J .* (sx .* rx + sy .* ry)
  css = J .* (sx .* sx + sy .* sy)

  #
  # Block surface matrices
  #
  (xf1, yf1) = (view(x, 1, :), view(y, 1, :))
  nx1 = -ys[1, :]
  ny1 =  xs[1, :]
  sJ1 = hypot.(nx1, ny1)
  nx1 = nx1 ./ sJ1
  ny1 = ny1 ./ sJ1

  (xf2, yf2) = (view(x, Nrp, :), view(y, Nrp, :))
  nx2 =  ys[end, :]
  ny2 = -xs[end, :]
  sJ2 = hypot.(nx2, ny2)
  nx2 = nx2 ./ sJ2
  ny2 = ny2 ./ sJ2

  (xf3, yf3) = (view(x, :, 1), view(y, :, 1))
  nx3 =  yr[:, 1]
  ny3 = -xr[:, 1]
  sJ3 = hypot.(nx3, ny3)
  nx3 = nx3 ./ sJ3
  ny3 = ny3 ./ sJ3

  (xf4, yf4) = (view(x, :, Nsp), view(y, :, Nsp))
  nx4 = -yr[:, end]
  ny4 =  xr[:, end]
  sJ4 = hypot.(nx4, ny4)
  nx4 = nx4 ./ sJ4
  ny4 = ny4 ./ sJ4


  (coord = (x,y),
   facecoord = ((xf1, xf2, xf3, xf4), (yf1, yf2, yf3, yf4)),
   crr = crr, css = css, crs = crs,
   J=J,
   sJ = (sJ1, sJ2, sJ3, sJ4),
   nx = (nx1, nx2, nx3, nx4),
   ny = (ny1, ny2, ny3, ny4),
   rx = rx, ry = ry, sx = sx, sy = sy)
end

function locoperator(p, Nr, Ns, metrics=create_metrics(p,Nr,Ns),
                     LFToB = (BC_DIRICHLET, BC_DIRICHLET,
                              BC_DIRICHLET, BC_DIRICHLET);
                     τscale = 2,
                     crr = metrics.crr,
                     css = metrics.css,
                     crs = metrics.crs)
  csr = crs
  J = metrics.J

  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  # Derivative operators for the rest of the computation
  (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
  Qr = Hr * Dr
  QrT = sparse(transpose(Qr))

  (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))
  Qs = Hs * Ds
  QsT = sparse(transpose(Qs))

  # Identity matrices for the comuptation
  Ir = sparse(I, Nrp, Nrp)
  Is = sparse(I, Nsp, Nsp)

  #{{{ Set up the rr derivative matrix
  ISr0 = Array{Int64,1}(undef,0)
  JSr0 = Array{Int64,1}(undef,0)
  VSr0 = Array{Float64,1}(undef,0)
  ISrN = Array{Int64,1}(undef,0)
  JSrN = Array{Int64,1}(undef,0)
  VSrN = Array{Float64,1}(undef,0)

  (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Nr, rand(Nrp))
  IArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
  JArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
  VArr = Array{Float64,1}(undef,Nsp * length(Ae.nzval))
  stArr = 0

  ISr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
  JSr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
  VSr0 = Array{Float64,1}(undef,Nsp * length(S0e.nzval))
  stSr0 = 0

  ISrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
  JSrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
  VSrN = Array{Float64,1}(undef,Nsp * length(SNe.nzval))
  stSrN = 0
  for j = 1:Nsp
    rng = (j-1) * Nrp .+ (1:Nrp)
    (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Nr, crr[rng])
    (Ie, Je, Ve) = findnz(Ae)
    IArr[stArr .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JArr[stArr .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VArr[stArr .+ (1:length(Ve))] = Hs[j,j] * Ve
    stArr += length(Ve)

    (Ie, Je, Ve) = findnz(S0e)
    ISr0[stSr0 .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JSr0[stSr0 .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VSr0[stSr0 .+ (1:length(Ve))] =  Hs[j,j] * Ve
    stSr0 += length(Ve)

    (Ie, Je, Ve) = findnz(SNe)
    ISrN[stSrN .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JSrN[stSrN .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VSrN[stSrN .+ (1:length(Ve))] =  Hs[j,j] * Ve
    stSrN += length(Ve)
  end
  Ãrr = sparse(IArr[1:stArr], JArr[1:stArr], VArr[1:stArr], Np, Np)
  Sr0 = sparse(ISr0[1:stSr0], JSr0[1:stSr0], VSr0[1:stSr0], Np, Np)
  SrN = sparse(ISrN[1:stSrN], JSrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  Sr0T = sparse(JSr0[1:stSr0], ISr0[1:stSr0], VSr0[1:stSr0], Np, Np)
  SrNT = sparse(JSrN[1:stSrN], ISrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  #= affine mesh test
  # @assert Ãrr ≈ Ãrr'
  (D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, Nr)
  Ar = SN - S0 - Hr * D2
  @assert Ãrr ≈ Hs ⊗ Ar
  =#
  # @assert Sr0 ≈ ((sparse(Diagonal(crr[1   .+ Nrp*(0:Ns)])) * Hs) ⊗ S0)
  # @assert SrN ≈ ((sparse(Diagonal(crr[Nrp .+ Nrp*(0:Ns)])) * Hs) ⊗ SN)
  #}}}

  #{{{ Set up the ss derivative matrix
  (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, rand(Nsp))
  IAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
  JAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
  VAss = Array{Float64,1}(undef,Nrp * length(Ae.nzval))
  stAss = 0

  ISs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
  JSs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
  VSs0 = Array{Float64,1}(undef,Nrp * length(S0e.nzval))
  stSs0 = 0

  ISsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
  JSsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
  VSsN = Array{Float64,1}(undef,Nrp * length(SNe.nzval))
  stSsN = 0
  for i = 1:Nrp
    rng = i .+ Nrp * (0:Ns)
    (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, css[rng])
    R = Ae - Dr' * Hr * Diagonal(css[rng]) * Dr

    (Ie, Je, Ve) = findnz(Ae)
    IAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VAss[stAss .+ (1:length(Ve))] = Hr[i,i] * Ve
    stAss += length(Ve)

    (Ie, Je, Ve) = findnz(S0e)
    ISs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JSs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VSs0[stSs0 .+ (1:length(Ve))] = Hr[i,i] * Ve
    stSs0 += length(Ve)

    (Ie, Je, Ve) = findnz(SNe)
    ISsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JSsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VSsN[stSsN .+ (1:length(Ve))] = Hr[i,i] * Ve
    stSsN += length(Ve)
  end
  Ãss = sparse(IAss[1:stAss], JAss[1:stAss], VAss[1:stAss], Np, Np)
  Ss0 = sparse(ISs0[1:stSs0], JSs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsN = sparse(ISsN[1:stSsN], JSsN[1:stSsN], VSsN[1:stSsN], Np, Np)
  Ss0T = sparse(JSs0[1:stSs0], ISs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsNT = sparse(JSsN[1:stSsN], ISsN[1:stSsN], VSsN[1:stSsN], Np, Np)
  # @assert Ãss ≈ Ãss'
  #= affine mesh test
  (D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, Ns)
  As = SN - S0 - Hs * D2
  @assert Ãss ≈ As ⊗ Hr
  =#
  # @assert Ss0 ≈ (S0 ⊗ (Hr * sparse(Diagonal(css[1:Nrp]))))
  # @assert SsN ≈ (SN ⊗ (Hr * sparse(Diagonal(css[Nrp*Ns .+ (1:Nrp)]))))
  #}}}

  #{{{ Set up the sr and rs derivative matrices
  Ãsr = (QsT ⊗ Ir) * sparse(1:length(crs), 1:length(crs), view(crs, :)) * (Is ⊗ Qr)
  Ãrs = (Is ⊗ QrT) * sparse(1:length(csr), 1:length(csr), view(csr, :)) * (Qs ⊗ Ir)
  #}}}

  Ã = Ãrr + Ãss + Ãrs + Ãsr

  #
  # Boundary point matrices
  #
  Er0 = sparse([1], [1], [1], Nrp, Nrp)
  ErN = sparse([Nrp], [Nrp], [1], Nrp, Nrp)
  Es0 = sparse([1], [1], [1], Nsp, Nsp)
  EsN = sparse([Nsp], [Nsp], [1], Nsp, Nsp)

  er0 = sparse([1  ], [1], [1], Nrp, 1)
  erN = sparse([Nrp], [1], [1], Nrp, 1)
  es0 = sparse([1  ], [1], [1], Nsp, 1)
  esN = sparse([Nsp], [1], [1], Nsp, 1)

  er0T = sparse([1], [1  ], [1], 1, Nrp)
  erNT = sparse([1], [Nrp], [1], 1, Nrp)
  es0T = sparse([1], [1  ], [1], 1, Nsp)
  esNT = sparse([1], [Nsp], [1], 1, Nsp)

  #
  # Store coefficient matrices as matrices
  #
  crs0 = sparse(Diagonal(crs[1:Nrp]))
  crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
  csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
  csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

  #
  # Surface mass matrices
  #
  H1 = Hs
  H1I = HsI

  H2 = Hs
  H2I = HsI

  H3 = Hr
  H3I = HrI

  H4 = Hr
  H4I = HrI

  #
  # Penalty terms
  #
  if p == 2
    l = 2
    β = 0.363636363
    α = 1 / 2
  elseif p == 4
    l = 4
    β = 0.2505765857
    α = 17 / 48
  elseif p == 6
    l = 7
    β = 0.1878687080
    α = 13649 / 43200
  else
    error("unknown order")
  end

  ψmin = reshape((crr + css - sqrt.((crr - css).^2 + 4crs.^2)) / 2, Nrp, Nsp)
  @assert minimum(ψmin) > 0

  hr = 2 / Nr
  hs = 2 / Ns

  ψ1 = ψmin[  1, :]
  ψ2 = ψmin[Nrp, :]
  ψ3 = ψmin[:,   1]
  ψ4 = ψmin[:, Nsp]
  for k = 2:l
    ψ1 = min.(ψ1, ψmin[k, :])
    ψ2 = min.(ψ2, ψmin[Nrp+1-k, :])
    ψ3 = min.(ψ3, ψmin[:, k])
    ψ4 = min.(ψ4, ψmin[:, Nsp+1-k])
  end
  τ1 = (2τscale / hr) * (crr[  1, :].^2 / β + crs[  1, :].^2 / α) ./ ψ1
  τ2 = (2τscale / hr) * (crr[Nrp, :].^2 / β + crs[Nrp, :].^2 / α) ./ ψ2
  τ3 = (2τscale / hs) * (css[:,   1].^2 / β + crs[:,   1].^2 / α) ./ ψ3
  τ4 = (2τscale / hs) * (css[:, Nsp].^2 / β + crs[:, Nsp].^2 / α) ./ ψ4

  τ1 = sparse(1:Nsp, 1:Nsp, τ1)
  τ2 = sparse(1:Nsp, 1:Nsp, τ2)
  τ3 = sparse(1:Nrp, 1:Nrp, τ3)
  τ4 = sparse(1:Nrp, 1:Nrp, τ4)

  C̃1 =  (Sr0 + Sr0T) + ((csr0 * Qs + QsT * csr0) ⊗ Er0) + ((τ1 * H1) ⊗ Er0)
  C̃2 = -(SrN + SrNT) - ((csrN * Qs + QsT * csrN) ⊗ ErN) + ((τ2 * H2) ⊗ ErN)
  C̃3 =  (Ss0 + Ss0T) + (Es0 ⊗ (crs0 * Qr + QrT * crs0)) + (Es0 ⊗ (τ3 * H3))
  C̃4 = -(SsN + SsNT) - (EsN ⊗ (crsN * Qr + QrT * crsN)) + (EsN ⊗ (τ4 * H4))

  # TODO: Fix minus sign (reverse of the paper)
  G1 = -(Is ⊗ er0T) * Sr0 - ((csr0 * Qs) ⊗ er0T)
  G2 = +(Is ⊗ erNT) * SrN + ((csrN * Qs) ⊗ erNT)
  G3 = -(es0T ⊗ Ir) * Ss0 - (es0T ⊗ (crs0 * Qr))
  G4 = +(esNT ⊗ Ir) * SsN + (esNT ⊗ (crsN * Qr))

  F1 = G1' - ((τ1 * H1) ⊗ er0)
  F2 = G2' - ((τ2 * H2) ⊗ erN)
  F3 = G3' - (es0 ⊗ (τ3 * H3))
  F4 = G4' - (esN ⊗ (τ4 * H4))

  M̃ = Ã + C̃1 + C̃2 + C̃3 + C̃4

  # Modify the operator to handle the boundary conditions
  bctype=(BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE)
  F = (F1, F2, F3, F4)
  τ = (τ1, τ2, τ3, τ4)
  HfI = (H1I, H2I, H3I, H4I)
  # Modify operators for the BC
  for lf = 1:4
    if LFToB[lf] == BC_NEUMANN
      M̃ -= F[lf] * (Diagonal(1 ./ (diag(τ[lf]))) * HfI[lf]) * F[lf]'
    elseif !(LFToB[lf] == BC_DIRICHLET ||
             LFToB[lf] == BC_LOCKED_INTERFACE ||
             LFToB[lf] >= BC_JUMP_INTERFACE)
      error("invalid bc")
    end
  end
  bctype=(LFToB[1], LFToB[2], LFToB[3], LFToB[4])

  # (E, V) = eigen(Matrix(M̃))
  # println((minimum(E), maximum(E)))
  JH = sparse(1:Np, 1:Np, view(J, :)) * (Hs ⊗ Hr)
  (M̃ = M̃,
   F = (F1, F2, F3, F4),
   coord = metrics.coord,
   facecoord = metrics.facecoord,
   JH = JH,
   sJ = metrics.sJ,
   nx = metrics.nx,
   ny = metrics.ny,
   Hf = (H1, H2, H3, H4),
   HfI = (H1I, H2I, H3I, H4I),
   τ = (τ1, τ2, τ3, τ4),
   bctype=bctype)
end
#}}}

#{{{ gloλoperator: Build the trace operators
function threaded_gloλoperator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
  nelems = length(lop)
  nfaces = length(FToB)
  Nλp = zeros(Int64, nfaces)
  FToλstarts = Array{Int64, 1}(undef, nfaces + 1)
  FToλstarts[1] = 1
  IT = Array{Int64,1}(undef,0)
  JT = Array{Int64,1}(undef,0)
  VT = Array{Float64,1}(undef,0)
  VD = Array{Float64,1}(undef,0)
  for f = 1:nfaces
    if FToB[f] == BC_DIRICHLET || FToB[f] == BC_NEUMANN
      FToλstarts[f+1] = FToλstarts[f]
      continue
    end
    (em, ep) = FToE[:, f]
    (fm, fp) = FToLF[:, f]
    Nλp[f] = (fm <= 2 ? Ns[em]+1 : Nr[em]+1)
    @assert Nλp[f] == (fp <= 2 ? Ns[ep]+1 : Nr[ep]+1)
    FToλstarts[f+1] = FToλstarts[f] + Nλp[f]

    @assert EToO[fm, em] && EToS[fm, em] == 1
    Fm = lop[em].F[fm]
    # swap I and J to get transpose
    (Je, Ie, Ve) = findnz(Fm)
    IT = [IT; Ie .+ (FToλstarts[f] - 1)]
    JT = [JT; Je .+ (vstarts[em] - 1)]
    VT = [VT; Ve]

    @assert EToS[fp, ep] == 2
    Fp = lop[ep].F[fp]
    # swap I and J to get transpose
    (Je, Ie, Ve) = findnz(Fp)
    # if element and face orientation do not match, then flip
    if EToO[fp, ep]
      IT = [IT; Ie .+ (FToλstarts[f] - 1)]
      τm = Vector(diag(lop[em].τ[fm]))
      τp = Vector(diag(lop[ep].τ[fp]))
    else
      IT = [IT; FToλstarts[f+1] .- Ie]
      τm = Vector(diag(lop[em].τ[fm]))
      τp = Vector(diag(rot180(lop[ep].τ[fp])))
    end
    JT = [JT; Je .+ (vstarts[ep] - 1)]
    VT = [VT; Ve]

    Hf = Vector(diag(lop[em].Hf[fm]))
    VD = [VD; Hf .* (τm + τp)]
  end
  λNp = FToλstarts[nfaces+1]-1
  VNp = vstarts[nelems+1]-1
  FbarT = sparse(IT, JT, VT, λNp, VNp)
  # Ttranspose = sparse(JT, IT, VT, VNp, λNp)
  (FToλstarts, FbarT, VD)
end
#}}}

#{{{ volbcarray()
function locbcarray!(ge, gδe, lop, LFToB, bc_Dirichlet, bc_Neumann, in_jump,
                     bcargs = ())
  F = lop.F
  (xf, yf) = lop.facecoord
  Hf = lop.Hf
  sJ = lop.sJ
  nx = lop.nx
  ny = lop.ny
  τ = lop.τ
  ge[:] .= 0
  for lf = 1:4
    if LFToB[lf] == BC_DIRICHLET
      vf = bc_Dirichlet(lf, xf[lf], yf[lf], bcargs...)
    elseif LFToB[lf] == BC_NEUMANN
      gN = bc_Neumann(lf, xf[lf], yf[lf], nx[lf], ny[lf], bcargs...)
      vf = sJ[lf] .* gN ./ diag(τ[lf])
    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here
    elseif LFToB[lf] >= BC_JUMP_INTERFACE
      # In this case we need to add in half the jump
      vf = in_jump(lf, xf[lf], yf[lf], bcargs...) / 2
      gδe[lf][:] -= Hf[lf] * τ[lf] * vf
    else
      error("invalid bc")
    end
    ge[:] -= F[lf] * vf
  end
end
#}}}


#{{{ volsourcearray()
function locsourcearray!(ge, source, lop, volargs = ())

  (xloc, yloc) = lop.coord
  JHloc = lop.JH
  ge[:] += JHloc * source(xloc[:], yloc[:], volargs...)

end

#}}}



#{{{
struct SBPLocalOperator1{T<:Real, S<:Factorization}
  offset::Array{Int64,1}
  H::Array{T,1}
  X::Array{T,1}
  Y::Array{T,1}
  E::Array{Int64,1}
  F::Array{S,1}
  SBPLocalOperator1{T,S}(vstarts::Array{Int64,1}, H::Array{T,1}, X::Array{T,1},
                         Y::Array{T,1}, E::Array{Int64,1},
                         F::Array{S,1}) where {T<:Real, S<:Factorization} =
  new(vstarts, H, X, Y, E, F)
end


struct SBPLocalOperator1_forming{T<:Real}
    offset::Array{Int64,1}
    H::Array{T,1}
    X::Array{T,1}
    Y::Array{T,1}
    E::Array{Int64,1}
    SBPLocalOperator1_forming{T}(vstarts::Array{Int64,1}, H::Array{T,1}, X::Array{T,1},
                           Y::Array{T,1}, E::Array{Int64,1}) where {T<:Real} =
    new(vstarts, H, X, Y, E)
end

function SBPLocalOperator1(lop, Nr, Ns, factorization)
  nelems = length(lop)
  vstarts = Array{Int64, 1}(undef, nelems + 1)
  vstarts[1] = 1
  Np = Array{Int64, 1}(undef, nelems)
  VH = Array{Float64,1}(undef,0)
  X = Array{Float64,1}(undef,0)
  Y = Array{Float64,1}(undef,0)
  E = Array{Int64,1}(undef,0)
  FTYPE = typeof(factorization(sparse([1],[1],[1.0])))
  factors = Array{FTYPE, 1}(undef, nelems)
  for e = 1:nelems
    # Fill arrays to build global sparse matrix
    Np[e] = (Nr[e]+1)*(Ns[e]+1)
    vstarts[e+1] = vstarts[e] + Np[e]

    # Global "mass" matrix
    JH = lop[e].JH
    VH = [VH;Vector(diag(JH))]

    # global coordinates and element number array (needed for jump)
    (x,y) = lop[e].coord
    X = [X;x[:]]
    Y = [Y;y[:]]
    E = [E;e * ones(Int64, Np[e])]

    factors[e] = factorization(lop[e].M̃)
  end
  VNp = vstarts[nelems+1]-1 # total number of volume points

  SBPLocalOperator1{Float64, FTYPE}(vstarts, VH, X, Y, E, factors)
end

function SBPLocalOperator1_forming(lop,Nr,Ns)
    nelems = length(lop)
    vstarts = Array{Int64, 1}(undef, nelems + 1)
    vstarts[1] = 1
    Np = Array{Int64, 1}(undef, nelems)
    VH = Array{Float64,1}(undef,0)
    X = Array{Float64,1}(undef,0)
    Y = Array{Float64,1}(undef,0)
    E = Array{Int64,1}(undef,0)
    # FTYPE = typeof(factorization(sparse([1],[1],[1.0])))
    # factors = Array{FTYPE, 1}(undef, nelems)
    for e = 1:nelems
      # Fill arrays to build global sparse matrix
      Np[e] = (Nr[e]+1)*(Ns[e]+1)
      vstarts[e+1] = vstarts[e] + Np[e]

      # Global "mass" matrix
      JH = lop[e].JH
      VH = [VH;Vector(diag(JH))]

      # global coordinates and element number array (needed for jump)
      (x,y) = lop[e].coord
      X = [X;x[:]]
      Y = [Y;y[:]]
      E = [E;e * ones(Int64, Np[e])]

      # factors[e] = factorization(lop[e].M̃)
    end
    VNp = vstarts[nelems+1]-1 # total number of volume points

    SBPLocalOperator1_forming{Float64}(vstarts, VH, X, Y, E)
end


function SBPLocalOperator1_factorization()

end


function threaded_SBPLocalOperator1(lop, Nr, Ns, factorization)
  nelems = length(lop)
  vstarts = Array{Int64, 1}(undef, nelems + 1)
  vstarts[1] = 1
  Np = Array{Int64, 1}(undef, nelems)
  VH = Array{Float64,1}(undef,0)
  X = Array{Float64,1}(undef,0)
  Y = Array{Float64,1}(undef,0)
  E = Array{Int64,1}(undef,0)
  FTYPE = typeof(factorization(sparse([1],[1],[1.0])))
  factors = Array{FTYPE, 1}(undef, nelems)
  @threads for e = 1:nelems
    # Fill arrays to build global sparse matrix
    Np[e] = (Nr[e]+1)*(Ns[e]+1)
    vstarts[e+1] = vstarts[e] + Np[e]

    # Global "mass" matrix
    JH = lop[e].JH
    VH = [VH;Vector(diag(JH))]

    # global coordinates and element number array (needed for jump)
    (x,y) = lop[e].coord
    X = [X;x[:]]
    Y = [Y;y[:]]
    E = [E;e * ones(Int64, Np[e])]

    factors[e] = factorization(lop[e].M̃)
  end
  VNp = vstarts[nelems+1]-1 # total number of volume points

  SBPLocalOperator1{Float64, FTYPE}(vstarts, VH, X, Y, E, factors)
end



#}}}
function LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                              factorization)
  M = SBPLocalOperator1(lop, Nr, Ns, factorization)
  (FToλstarts, FbarT, D) = threaded_gloλoperator(lop, M.offset, FToB, FToE, FToLF, EToO,
                                        EToS, Nr, Ns)
  (M, FbarT, D, M.offset, FToλstarts)
end

function LocalGlobalOperators_forming(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS)
    M = SBPLocalOperator1_forming(lop,Nr,Ns)
    (M)
end

function threaded_LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                              factorization)
  M = SBPLocalOperator1(lop, Nr, Ns, factorization)
  (FToλstarts, FbarT, D) = threaded_gloλoperator(lop, M.offset, FToB, FToE, FToLF, EToO,
                                        EToS, Nr, Ns)
  (M, FbarT, D, M.offset, FToλstarts)
end

function bcstarts(FToB, FToE, FToLF, bctype, Nr, Ns)
  nfaces = length(FToB)
  bcstarts = Array{Int64, 1}(undef, nfaces + 1)
  bcstarts[1] = 1
  for f = 1:nfaces
    if FToB[f] ∈ bctype
      e  = FToE[1,f]
      lf = FToLF[1,f]
      bcstarts[f+1] = bcstarts[f] + (lf ∈ (1,2) ? Ns[e] : Nr[e]) + 1
    else
      bcstarts[f+1] = bcstarts[f]
    end
  end
  bcstarts
end

function LocalToGLobalRHS!(b, g, gδ, u, M, FbarT, vstarts)
  for e = 1:length(M)
    @views u[vstarts[e]:(vstarts[e+1]-1)] = (M[e] \
                                             g[vstarts[e]:(vstarts[e+1]-1)]) # It does back slash here again
  end
  mul!(b, FbarT, u)
  @. b = gδ - b
end

function threaded_LocalToGLobalRHS!(b, g, gδ, u, M, FbarT, vstarts)
  @threads for e = 1:length(M)
    @views u[vstarts[e]:(vstarts[e+1]-1)] = (M[e] \
                                             g[vstarts[e]:(vstarts[e+1]-1)]) # It does back slash here again
  end
  mul!(b, FbarT, u)
  @. b = gδ - b
end


#{{{ assembleλmatrix: Schur complement system
function assembleλmatrix(FToλstarts, vstarts, EToF, FToB, F, D, FbarT)
  nfaces = length(FToλstarts)-1
  nelems = length(vstarts)-1
  λNp = FToλstarts[nfaces+1]-1
  sz = λNp

  for e = 1:nelems
    lλs = Array{Int64, 1}(undef, 4)
    for lf = 1:4
      f = EToF[lf,e]
      lλs[lf] = FToλstarts[f+1] - FToλstarts[f]
    end
    for lf = 1:4
      sz += lλs[lf]*sum(lλs)
    end
  end
  Ie = Array{Int64, 1}(undef, sz)
  Je = Array{Int64, 1}(undef, sz)
  Ve = Array{Float64, 1}(undef, sz)
  Ie[1:λNp] = 1:λNp
  Je[1:λNp] = 1:λNp
  Ve[1:λNp] = D
  offset = λNp
  Fbar = FbarT'
  for e = 1:nelems
    # println((e, nelems))
    vrng = vstarts[e]:(vstarts[e+1]-1)
    for lf = 1:4
      f = EToF[lf,e]
      if FToB[f] == BC_LOCKED_INTERFACE || FToB[f] >= BC_JUMP_INTERFACE
        λrng = FToλstarts[f]:(FToλstarts[f+1]-1)
        B = -(Matrix(F[e]' \ Fbar[vrng, λrng])) # This is where backslash happens
        for lf2 = 1:4
          f2 = EToF[lf2,e]
          if FToB[f2] == BC_LOCKED_INTERFACE || FToB[f2] >= BC_JUMP_INTERFACE
            λrng2 = FToλstarts[f2]:(FToλstarts[f2+1]-1)
            C = -(FbarT[λrng2, vrng] * B)
            λblck = λrng*ones(Int64, 1, length(λrng2))
            λblck2 = ones(Int64, length(λrng), 1) * λrng2'
            last = length(λrng) * length(λrng2)
            Ie[offset.+(1:last)] = λblck[:]
            Je[offset.+(1:last)] = λblck2[:]
            Ve[offset.+(1:last)] = -C'[:]
            offset += last
          end
        end
      end
    end
  end
  @assert offset == sz
  B = sparse(Ie, Je, Ve, λNp, λNp)
  @assert B ≈ B'
  # println((λNp * λNp, nnz(B), nnz(B) / λNp^2))
  B
end


#{{{ assembleλmatrix: Schur complement system
function threaded_assembleλmatrix(FToλstarts, vstarts, EToF, FToB, F, D, FbarT)
  nfaces = length(FToλstarts)-1
  nelems = length(vstarts)-1
  λNp = FToλstarts[nfaces+1]-1
  sz = λNp

  for e = 1:nelems
    lλs = Array{Int64, 1}(undef, 4)
    for lf = 1:4
      f = EToF[lf,e]
      lλs[lf] = FToλstarts[f+1] - FToλstarts[f]
    end
    for lf = 1:4
      sz += lλs[lf]*sum(lλs)
    end
  end
  Ie = Array{Int64, 1}(undef, sz)
  Je = Array{Int64, 1}(undef, sz)
  Ve = Array{Float64, 1}(undef, sz)
  Ie[1:λNp] = 1:λNp
  Je[1:λNp] = 1:λNp
  Ve[1:λNp] = D
  offset = λNp
  Fbar = FbarT'
  for e = 1:nelems
    # println((e, nelems))
    vrng = vstarts[e]:(vstarts[e+1]-1)
    for lf = 1:4
      f = EToF[lf,e]
      if FToB[f] == BC_LOCKED_INTERFACE || FToB[f] >= BC_JUMP_INTERFACE
        λrng = FToλstarts[f]:(FToλstarts[f+1]-1)
        B = -(Matrix(F[e]' \ Fbar[vrng, λrng])) # This is where backslash happens
        for lf2 = 1:4
          f2 = EToF[lf2,e]
          if FToB[f2] == BC_LOCKED_INTERFACE || FToB[f2] >= BC_JUMP_INTERFACE
            λrng2 = FToλstarts[f2]:(FToλstarts[f2+1]-1)
            C = -(FbarT[λrng2, vrng] * B)
            λblck = λrng*ones(Int64, 1, length(λrng2))
            λblck2 = ones(Int64, length(λrng), 1) * λrng2'
            last = length(λrng) * length(λrng2)
            Ie[offset.+(1:last)] = λblck[:]
            Je[offset.+(1:last)] = λblck2[:]
            Ve[offset.+(1:last)] = -C'[:]
            offset += last
          end
        end
      end
    end
  end
  @assert offset == sz
  B = sparse(Ie, Je, Ve, λNp, λNp)
  @assert B ≈ B'
  # println((λNp * λNp, nnz(B), nnz(B) / λNp^2))
  B
end

#}}}

# {{{ Constructor for inp files
function read_inp_2d(T, S, filename::String; bc_map=1:10000)
  # {{{ Read in the file
  f = try
    open(filename)
  catch
    error("InpRead cannot open \"$filename\" ")
  end
  lines = readlines(f)
  close(f)
  # }}}

  # {{{ Read in nodes
  str = "NSET=ALLNODES"
  linenum = SeekToSubstring(lines, str);
  linenum > 0 || error("did not find: $str")
  num_nodes = 0
  for l = linenum+1:length(lines)
    occursin(r"^\s*[0-9]*\s*,.*", lines[l]) ? num_nodes+=1 : break
  end
  Vx = fill(S(NaN), num_nodes)
  Vy = fill(S(NaN), num_nodes)
  Vz = fill(S(NaN), num_nodes)
  init_node_num = -1
  for l = linenum .+ (1:num_nodes)
    node_data = split(lines[l], r"\s|,", keepempty=false)
    (node_num, node_x, node_y, node_z) = try
      (parse(T, node_data[1]),
       parse(S, node_data[2]),
       parse(S, node_data[3]),
       parse(S, node_data[4]))
    catch
      error("cannot parse line $l: \"$(lines[l])\" ")
    end

    if (init_node_num == -1)
        init_node_num = node_num
    end
    node_num = node_num - init_node_num + 1
    Vx[node_num] = node_x
    Vy[node_num] = node_y
    Vz[node_num] = node_z
  end
  # }}}

  # {{{ Read in Elements
  str = "ELEMENT"
  linenum = SeekToSubstring(lines, str);
  num_elm = 0
  while linenum > 0
    for l = linenum .+ (1:length(lines))
      occursin(r"^\s*[0-9]*\s*,.*", lines[l]) ? num_elm+=1 : break
    end
    linenum = SeekToSubstring(lines, str; first=linenum+1)
  end
  num_elm > 0 || error("did not find any element")

  EToV = fill(T(0), 4, num_elm)
  EToBlock = fill(T(0), num_elm)
  linenum = SeekToSubstring(lines, str);
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    B = parse(T, foo[end])
    for l = linenum .+ (1:num_elm)
      elm_data = split(lines[l], r"\s|,", keepempty=false)
      # read into z-order
      (elm_num, elm_v1, elm_v2, elm_v4, elm_v3) = try
        (parse(T, elm_data[1]),
         parse(T, elm_data[2]),
        parse(T, elm_data[3]),
        parse(T, elm_data[4]),
        parse(T, elm_data[5]))
      catch
        break
      end

      EToV[:, elm_num] = [elm_v1, elm_v2, elm_v3, elm_v4] .- init_node_num .+ 1
      EToBlock[elm_num] = B
    end
    linenum = SeekToSubstring(lines, str; first=linenum+1)
  end
  # }}}

  # {{{ Determine connectivity
  EToF = fill(T(0), 4, num_elm)

  VsToF = Dict{Tuple{Int64, Int64}, Int64}()
  numfaces = 0
  for e = 1:num_elm
    for lf = 1:4
      if lf == 1
        Vs = (EToV[1, e], EToV[3, e])
      elseif lf == 2
        Vs = (EToV[2, e], EToV[4, e])
      elseif lf == 3
        Vs = (EToV[1, e], EToV[2, e])
      elseif lf == 4
        Vs = (EToV[3, e], EToV[4, e])
      end
      if Vs[1] > Vs[2]
        Vs = (Vs[2], Vs[1])
      end
      if haskey(VsToF, Vs)
        EToF[lf, e] = VsToF[Vs]
      else
        numfaces = numfaces + 1
        EToF[lf, e] = VsToF[Vs] = numfaces
      end
    end
  end
  #}}}

  # {{{ Read in side set info
  FToB = Array{T, 1}(undef, numfaces)
  fill!(FToB, BC_LOCKED_INTERFACE)
  linenum = SeekToSubstring(lines, "\\*ELSET")
  inp_to_zorder = [3,  2, 4, 1]
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    (bc, face) = try
      (parse(T, foo[1]),
       parse(T, foo[2]))
    catch
      error("cannot parse line $linenum: \"$(lines[linenum])\" ")
    end
    bc = bc_map[bc]
    face = inp_to_zorder[face]
    for l = linenum+1:length(lines)
      if !occursin(r"^\s*[0-9]+", lines[l])
        break
      end
      elms = split(lines[l], r"\s|,", keepempty=false)
      for elm in elms
        elm = try
          parse(T, elm)
        catch
          error("cannot parse line $linenum: \"$(lines[l])\" ")
        end
        if bc == 3
          bc = BC_LOCKED_INTERFACE
        end
        FToB[EToF[face, elm]] = bc
        @assert (bc == BC_DIRICHLET || bc == BC_NEUMANN ||
                 bc == BC_LOCKED_INTERFACE || bc >= BC_JUMP_INTERFACE)
      end
    end
    linenum = SeekToSubstring(lines, "\\*ELSET"; first=linenum+1)
  end
  # }}}

  ([Vx Vy]', EToV, EToF, FToB, EToBlock)
end
read_inp_2d(filename;kw...) = read_inp_2d(Int64, Float64, filename;kw...)

function SeekToSubstring(lines, substring; first=1)
  for l = first:length(lines)
    if occursin(Regex(".*$(substring).*"), lines[l])
      return l
    end
  end
  return -1
end

# }}}

function plot_connectivity(verts, EToV)
  Lx = extrema(verts[1,:])
  Lx = (floor(Int, Lx[1]), ceil(Int, Lx[2]))
  Ly = extrema(verts[2,:])
  Ly = (floor(Int, Ly[1]), ceil(Int, Ly[2]))
  plt = Plot(BrailleCanvas(80, 40,
                           origin_x = Lx[1], origin_y = Ly[1],
                           width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))


  annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt, :br, string(Lx[2]), color = :light_black)
  for e = 1:size(EToV, 2)
    (v1, v2, v3, v4) = EToV[1:4, e]
    x = verts[1, [v1 v2 v4 v3 v1]][:]
    y = verts[2, [v1 v2 v4 v3 v1]][:]
    lineplot!(plt, x, y)
  end
  title!(plt, "connectivity")
  display(plt)
end

function plot_blocks(lop)
  Lx = (floatmax(), -floatmax())
  Ly = (floatmax(), -floatmax())
  for e = 1:length(lop)
    (x, y) = lop[e].coord
    Lxe = extrema(x)
    Lye = extrema(y)
    Lx = (min(Lx[1], Lxe[1]), max(Lx[2], Lxe[2]))
    Ly = (min(Ly[1], Lye[1]), max(Ly[2], Lye[2]))
  end

  Lx = (floor(Int, Lx[1]), ceil(Int, Lx[2]))
  Ly = (floor(Int, Ly[1]), ceil(Int, Ly[2]))

  plt = Plot(BrailleCanvas(80, 40,
                           origin_x = Lx[1], origin_y = Ly[1],
                           width = Lx[2] - Lx[1], height = Lx[2] - Lx[1]))


  annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt, :br, string(Lx[2]), color = :light_black)

  for e = 1:length(lop)
    (xf, yf) = lop[e].facecoord
    bctype = lop[e].bctype
    for lf = 1:length(xf)
      if bctype[lf] == BC_LOCKED_INTERFACE
        lineplot!(plt, xf[lf], yf[lf], color=:blue)
      elseif bctype[lf] == BC_DIRICHLET
        lineplot!(plt, xf[lf], yf[lf], color=:green)
      elseif bctype[lf] == BC_DIRICHLET
        lineplot!(plt, xf[lf], yf[lf], color=:yellow)
      else
        lineplot!(plt, xf[lf], yf[lf], color=:red)
      end
    end
  end
  title!(plt, "mesh")
  display(plt)
end
