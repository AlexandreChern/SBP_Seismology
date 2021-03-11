using DataFrames
using CSV

# cd("BP1")
f = CSV.File("BP1_30.0.dat")
f.shear_stress .= f.shear_stress .- 4.62444 .* 10 .^(f.slip_rate)

CSV.write("BP1_30.0_v2.dat",f,delim=' ')