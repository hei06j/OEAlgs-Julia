# import Pkg; Pkg.add("OpenDSSDirect")
using CSV
using DataFrames
using Polynomials

const DFs = DataFrames

# Load necessary inputs

tx_rated_capacity = 500

names_active_cust = CSV.File("/Users/abond/OEAlgs-Julia/day_30-15/names_active_cust.csv") |> DFs.DataFrame # active customers in the network
active_cust_der_size = ones(DFs.nrow(names_active_cust)) * 10.0

v_tx_max_exp = CSV.File("/Users/abond/OEAlgs-Julia/day_30-15/v_tx_max_exp.csv") |> DFs.DataFrame # load the max transformer voltage 
delta_v_exp = CSV.File("/Users/abond/OEAlgs-Julia/day_30-15/delta_v_exp.csv") |> DFs.DataFrame # load the delta v of the critical customer
p_agg_tx_exp = CSV.File("/Users/abond/OEAlgs-Julia/day_30-15/p_agg_tx_exp.csv") |> DFs.DataFrame # load the aggregated power of the transformer

forecast_spare_capacity_exp = CSV.File("/Users/abond/OEAlgs-Julia/day_30-15/forecast_spare_capacity_exp.csv") |> DFs.DataFrame # load the spare capacity of the distribution transformer for exports
forecast_agg_net_p_passive_cust = CSV.File("/Users/abond/OEAlgs-Julia/day_30-15/forecast_agg_net_p_passive_cust.csv") |> DFs.DataFrame # load the net power for passive customers

# Create sensitivity curves

pdtx_deltav_sensit_curve_exp = Polynomials.fit(p_agg_tx_exp[!,1], delta_v_exp[!,1], 1) # export
pdtx_vdtx_sensit_curve_exp = Polynomials.fit(p_agg_tx_exp[!,1], v_tx_max_exp[!,1], 1) # export

function ac_deltav_oe_exp_algorithm(
    names_active_cust,
    tx_spare_capacity_exp,
    active_cust_der_size, 
    agg_net_p_passive_cust,
    pdtx_deltav_sensit_curve_exp,
    pdtx_vdtx_sensit_curve_exp)
    """
    Function to calculate the Asset Capacity & Delta Voltage OE proportional allocation for exports in each time step.
    This is used in all time steps of the day.    
    :param names_active_cust: the list of active customers 
    :param tx_spare_capacity_exp: the distribution transformer spare export capacity of the considered time step
    :param active_cust_der_size: the size of active customers' DERs 
    :param agg_net_p_passive_cust: the forecasted aggregated passive customer net demand
    :param pdtx_deltav_sensit_curve_exp: the PDTx-∆V sensitivity curve for exports
    :param pdtx_vdtx_sensit_curve_exp: the PDTx-VDTx sensitivity curve for exports
    :return ac_deltav_oe_exp_values: the calculated AC_∆V OE value for proportional allocation and exports
    """
    
    # First main step: the capacity of the distribution transformer together with its time-varying input data 
    # is used to calculate its spare capacity
    alloc_exp_p_temp = zeros(DFs.nrow(names_active_cust)) # initialisation of the operating envelope value 
    
    # check if there is spare capacity on the distribution transformer, if positive there is spare capacity,
    # if negative or zero the allocated operating envelope is equal to zero (same value as initialisation)
    if tx_spare_capacity_exp > 0
        tx_spare_capacity_temp = tx_spare_capacity_exp # copy the spare capacity to a temporary variable
        
        # run this until there is no more spare capacity or the maximum possible OE value (defined by connection agreement,
        # or fuse of the house) was achieved to all active customers
        while tx_spare_capacity_temp > 0
            for ild in range(1,DFs.nrow(names_active_cust)) # iterate through all active customers to calculate the OE value
                
                # Proportionally allocate spare capacity of the distribution transformer to each active customer
                alloc_exp_p_temp[ild] += tx_spare_capacity_temp * active_cust_der_size[ild] / sum(active_cust_der_size) 
                
                # Check is the allocated spare capacity pass the maximum possible OE value and cap it if needed
                if alloc_exp_p_temp[ild] > 10  # Usually maximum injection allowed by DNSPs per phase (10kW)
                    alloc_exp_p_temp[ild] = 10
                end
            end
            
            # check if the maximum possible OE value (10kW) was achieved to all active customers
            if all(i->(i==10),alloc_exp_p_temp)
                break # the maximum possible OE value (10kW) was achieved to all active customers
            else
                
                # there are still some active customers below the maximum possible OE value 
                # check if there is still some spare capacity to be allocated
                tx_spare_capacity_temp = tx_spare_capacity_exp - sum(alloc_exp_p_temp) 
            end
        end    
        # Second main step: both PDTx-VDTx and PDTx-∆V sensitivity curves are used to estimate the voltage at the 
        # critical customer for a given aggregated active power passing through the DTx
        flag = 0 # flag to indicate when there is a voltage problem (0 = problems; 1 = no problems)
  
        # run the loop until no voltage problem is expected
        while flag == 0
            
            # calculate the expected aggregated active power at the transformer
            # note that 12 is the location of the critical customer in the array
            expected_p_tx = agg_net_p_passive_cust + DFs.nrow(names_active_cust) * alloc_exp_p_temp[12] 
                        
            # use the PDTx-VDTx sensitivity curve to estimate the voltage at the transformer for the expected 
            # aggregated active power at the transformer
            expected_v_tx = pdtx_vdtx_sensit_curve_exp(expected_p_tx)
            
            # use the PDTx-∆V sensitivity curve to estimate the delta voltage between the critical customer and the 
            # distribution transformer for the expected aggregated active power at the transformer
            expected_delta_v = pdtx_deltav_sensit_curve_exp(expected_p_tx)
            
            # use the expected voltage at the distribution transformer and the expected delta voltage to calculate 
            # the expected voltage at the critical customer sensitivity curve to estimate the voltage at the critical customer
            expected_v_crit_cust = expected_v_tx + expected_delta_v
            
            # check if the expected voltage at the critical customer is above 253V  
            if expected_v_crit_cust > 253
                alloc_exp_p_temp[12] -= 0.5 # if above the limit, reduce the allocated OE by 0.5kW (pre-defined reduction step)
                
                # check if the OE has achieved the minimum possible OE value (zero)
                if alloc_exp_p_temp[12] < 0
                    alloc_exp_p_temp[12] = 0 # if below zero, allocate zero as OE value
                    
                    # voltage problems are still expected even after OE value is equal to zero. OE calculation finished.
                    break
                end
            else
                flag = 1 # no voltage problems are expected. OE calculation finished.
            end
        end
    end
    ac_deltav_oe_exp_values = alloc_exp_p_temp[12] # rename the allocated OE value
    
    return ac_deltav_oe_exp_values # return the calculated OE value for the export for the considered time step
end

# define number of time steps in the day (24h) for the corresponding time resolution
num_of_time_step = 288

# initialise arrays to save the OE values for each time step of the day
ac_deltav_oe_day_exp_values = zeros(DFs.nrow(names_active_cust), num_of_time_step) # for exports

# Calculate the OE values for each time step of the day in a for loop
for itime in range(1,num_of_time_step)
    # println(itime,forecast_spare_capacity_exp[!,1][1])
    tx_spare_capacity_exp = forecast_spare_capacity_exp[!,1][itime] # separate the transformer export capacity for the current time step
    agg_net_p_passive_cust = forecast_agg_net_p_passive_cust[!,1][itime] # separate the transformer import capacity for the current time step
    # call the AC CrV OE functions to calculate the OE value to each time step
    ac_deltav_oe_day_exp_values[:,itime] .= ac_deltav_oe_exp_algorithm(names_active_cust, tx_spare_capacity_exp, active_cust_der_size, agg_net_p_passive_cust, pdtx_deltav_sensit_curve_exp, pdtx_vdtx_sensit_curve_exp) # export
end

active_exp_val_path = "/Users/abond/OEAlgs-Julia/data/ac_deltav/active_exp_values.csv"
active = DFs.DataFrame(ac_deltav_oe_day_exp_values,Symbol.(Vector(range(1,num_of_time_step))))

CSV.write(active_exp_val_path, active)

using OpenDSSDirect

const DSS = OpenDSSDirect

# define time resolution of the data
time_resolution = 5 # in minutes

# define number of time steps in the day (24h) for the corresponding time resolution
num_of_time_step = 288

# path to master file
filename = "/Users/abond/OEAlgs-Julia/Master.txt"

dss("""
    Clear
    Compile "$filename"
    """)
DSS.Text.Command("Set VoltageBases=[22.0, 0.400, 0.2309]")
DSS.Text.Command("calcv")
DSS.Text.Command("Set ControlMode=static")
DSS.Text.Command("Reset")
DSS.Text.Command("Set Mode=daily number=1 stepsize=5m")

# load LV network data from the OpenDSS model
load_list = DSS.Loads.AllNames() # list of loads

# neutral conductor can go up to 2-5V 10% difference 
# phase - neutral conenction vs phase - ground. opendss is neutral by default
# algorithm should know what inverter controls are being used otherwise need to figure out using smart meter data

valid_exp_voltage_lv_cust = zeros(length(load_list), num_of_time_step)

basepath = "/Users/abond/OEAlgs-Julia/"
active_cust = joinpath(basepath, "names_active_cust.csv")
tx_angles = joinpath(basepath, "tx_p_angles.csv")
tx_voltages = joinpath(basepath, "tx_p_voltages.csv")

# read in data for inputs
names_active_cust = CSV.File(active_cust) |> DFs.DataFrame # active customers in the network
active_cust_der_size = ones(DFs.nrow(names_active_cust)) * 10.0

# load voltage magnitudes and angles at the primary side of the distribution transformer which are affected by the 
# interactions with the upstream HV network (collected in another platform where the HV-LV network is fully modelled)
tx_pri_voltages_day = CSV.File(tx_voltages) |> DFs.DataFrame # voltage magnitudes
tx_pri_angles_day = CSV.File(tx_angles) |> DFs.DataFrame # voltage angles

# use OpenDSS to calculate customer voltages
for itime in range(1,num_of_time_step)
    # Set vsource with the distribution transformer voltage magnitudes and angles for the time step at the primary side
    temp1 = tx_pri_voltages_day[!,1][itime] / (22000 / sqrt(3))
    temp2 = tx_pri_angles_day[!,1][itime]
    temp3 = tx_pri_voltages_day[!,2][itime] / (22000 / sqrt(3))
    temp4 = tx_pri_angles_day[!,2][itime]
    temp5 = tx_pri_voltages_day[!,3][itime] / (22000 / sqrt(3))
    temp6 = tx_pri_angles_day[!,3][itime]

    DSS.Text.Command("""edit vsource.source bus1=sourcebus.1 basekv=12.701706 pu="$temp1" angle="$temp2" phases=1""")
    DSS.Text.Command("""edit vsource.source2 bus1=sourcebus.2 basekv=12.701706 pu="$temp3" angle="$temp4" phases=1""")
    DSS.Text.Command("""edit vsource.source3 bus1=sourcebus.3 basekv=12.701706 pu="$temp5" angle="$temp6" phases=1""")

    DSS.Solution.Solve() # solve the power flow in OpenDSS

    # implement the calculated OE value to each active customer 
    for ild in range(1,DFs.nrow(names_active_cust))
        DSS.Circuit.SetActiveElement("load."*string(names_active_cust[!,1][ild]))
        temp = imag(DSS.CktElement.Powers()[1]) # save the current reactive power of the active customer
        DSS.Properties.Value("status","fixed") # fix the load status to run with the set value instead of profile
        DSS.Properties.Value("kw",string(-1*ac_deltav_oe_day_exp_values[ild, itime])) # set the kW to the calculated OE value
        DSS.Properties.Value("kvar",string(temp)) # set the kvar to be the same as before
    end
        
    DSS.Solution.SolveSnap() # solve the power flow for the time step with OE values in place

    # verify voltages in all customers
    all_volt_temp = [] # initialise array to collect the voltages
    reactive_p_temp = []
    for ild in range(1,length(load_list))
        DSS.Circuit.SetActiveElement("load."*string(load_list[ild])) # select the customer
        push!(reactive_p_temp, imag(DSS.CktElement.Powers()[1]))
        push!(all_volt_temp, DSS.CktElement.VoltagesMagAng()[1]) # extract its voltage magnitude
    end
     
    valid_exp_voltage_lv_cust[:, itime] = reactive_p_temp # save voltages from all customers for the current time step 

    # Reset active customers to profile values for the next time step
    default_kw_kvar = 1
    for ild in range(1,DFs.nrow(names_active_cust))
        DSS.Circuit.SetActiveElement("load."*string(names_active_cust[!,1][ild]))
        DSS.Properties.Value("status","variable")
        DSS.Properties.Value("kW",string(default_kw_kvar))
        DSS.Properties.Value("kvar",string(default_kw_kvar))
    end
end

#save customer voltage to csv for plotting
cust_volt_path = "/Users/abond/OEAlgs-Julia/data/ac_deltav/all_cust_reactive_values.csv"
cust_volts = DFs.DataFrame(valid_exp_voltage_lv_cust,Symbol.(Vector(range(1,num_of_time_step))))

CSV.write(cust_volt_path, cust_volts)

#visualise reactive power

# Check network-wide voltage compliance
valid_exp_voltage_lv_cust_10min = zeros(length(load_list), 144)
valid_exp_voltage_lv_cust_10min_sorted = zeros(length(load_list), 144)
global cont_exp = 0
for ild in range(1,length(load_list))
    for i in range(1,144)
        valid_exp_voltage_lv_cust_10min[ild,i] = (valid_exp_voltage_lv_cust[ild,2*i-1] + valid_exp_voltage_lv_cust[ild,2*i]) / 2
    end
    valid_exp_voltage_lv_cust_10min_sorted[ild,:] = sort(valid_exp_voltage_lv_cust_10min[ild,:])
    if (valid_exp_voltage_lv_cust_10min_sorted[ild,143] > 253) || (valid_exp_voltage_lv_cust_10min_sorted[ild,2] < 216)
        global cont_exp += 1
    end
end

exp_voltage_compliance = (1 - ((cont_exp) / length(load_list))) * 100

println("The voltage compliance is: "*exp_voltage_compliance)

println("The maximum export voltage was: "*maximum(valid_exp_voltage_lv_cust))

# graph all volt_cust in pandas easier

# delta v + asset capacity is more accurate because the sensitivity curves have slightly different slopes, so you get more accurate picture from data

# FINISH SCRIPT, DO PIPELINE WORK (LUKE's STUFF), SORT OUT DAILY DATA UPLOAD TO MARKUS' SERVER, NEED quarterly update for ee
