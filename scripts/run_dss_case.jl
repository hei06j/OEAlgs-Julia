using Pkg
Pkg.activate("./")
using OpenDSSDirect
using CSV
using DataFrames
using Plots

const DFs = DataFrames
const DSS = OpenDSSDirect


##
# basepath = joinpath(@__DIR__, "../data/Master_nando_Neutral")
# # basepath = joinpath(@__DIR__, "../data/Master_nando")

# filename = joinpath(basepath, "Master.dss")

# # define time resolution of the data
# time_resolution = 5 # in minutes

# dss("""
#     clear
#     compile "$filename"

#     set VoltageBases=[22.0, 0.400, 0.2309]
#     calcv
    
#     set ControlMode=static
#     reset

#     set Mode=snap number=1 stepsize="$time_resolution"m
#     """)

filename = "/Users/hei06j/Documents/Repositories/Remote/OEAlgs-Julia/data/Master_nando_Neutral/Master.dss"
dss("""
    compile "$filename"
    """)

@show DSS.Circuit.AllBusVMag()

plot(DSS.Circuit.AllBusVMag(), seriestype=:scatter)

# @show DSS.Bus.puVmagAngle()

##
ideal_or_active_customer_kws = []
ideal_or_active_customer_kvars = []
ideal_or_active_customer_voltages = []
for ild in range(1,nrow(names_active_cust))
    DSS.Circuit.SetActiveElement("load."*names_active_cust[!,1][ild]) # select an active customer
    push!(ideal_or_active_customer_kws, real(DSS.CktElement.Powers()[1])) # save the current reactive power of the active customer
    push!(ideal_or_active_customer_kvars, imag(DSS.CktElement.Powers()[1])) # save the current reactive power of the active customer
    push!(ideal_or_active_customer_voltages, abs.(DSS.CktElement.Voltages()[1]))
end
plot(ideal_or_active_customer_voltages, seriestype=:scatter)

plot(ideal_or_active_customer_voltages, ideal_or_active_customer_kvars, seriestype=:scatter)

