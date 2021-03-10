using DataFrames
using CSV

cd("BP1")
f = CSV.File("BP1_0.0.dat")
f.shear_stress .= f.shear_stress .- 4.62444 .* 10 .^(f.slip_rate)

CSV.write("BP1_0.0_v2.dat",f,delim=' ')

df = DataFrames.read("BP1_0.0.dat", seperator=' ')