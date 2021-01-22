const year_seconds = 31556926

using DifferentialEquations
using Printf

function odefun(dψV,ψδ,p,t)
    RS_FAULT = 7
    VP_FAULT = 8
    reject_step = p.reject_step
    Vp = p.Vp
    lop = p.lop
    
end