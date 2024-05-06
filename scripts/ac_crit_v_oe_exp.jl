using CSV
using DataFrames
using Polynomials

const DFs = DataFrames

# Load necessary inputs

tx_rated_capacity = 500

names_active_cust = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/names_active_cust.csv") |> DFs.DataFrame # active customers in the network
active_cust_der_size = ones(DFs.nrow(names_active_cust)) * 10.0

hist_voltage_crit_cust = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/hist_voltage_crit_cust.csv") |> DFs.DataFrame # load the historical voltage magnitudes of the critical customer
hist_net_demand_crit_cust = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/hist_net_demand_crit_cust.csv") |> DFs.DataFrame # load the historical net power demand of the critical customer

forecast_spare_capacity_exp = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/forecast_spare_capacity.csv") |> DFs.DataFrame # load the spare capacity of the distribution transformer for exports

# tx_pri_angles_day = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/tx_pri_angles.csv") |> DFs.DataFrame # load voltage angles at the primary side of the distribution transformer
# tx_pri_voltages_day = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/tx_pri_voltages.csv") |> DFs.DataFrame # load voltage magnitudes at the primary side of the distribution transformer

# Create sensitivity curve from historical data

p_v_sensit_curve = Polynomials.fit(hist_net_demand_crit_cust[!,1], hist_voltage_crit_cust[!,1], 1)

# Define AC OE Export Algorithm
function ac_crv_oe_exp_algorithm(
    names_active_cust,
    tx_spare_capacity_exp,
    active_cust_der_size,
    p_v_sensit_curve)
    """
    Function to calculate the Asset Capacity & Critical Voltage OE proportional allocation for exports in each time step.
    This is used in all time steps of the day.    
    :param names_active_cust: the list of active customers 
    :param tx_spare_capacity_exp: the distribution transformer spare export capacity of the considered time step
    :param active_cust_der_size: the size of active customers' DERs 
    :param p_v_sensit_curve: the P-V sensitivity curve of the critical customer
    :return ac_crv_oe_exp_values: the calculated AC_CrV OE value for proportional allocation and exports
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
            
            # check if the maximum possible OE value (10kW) was achieved to all acticve customers
            if all(i->(i==10),alloc_exp_p_temp)
                break # the maximum possible OE value (10kW) was achieved to all acticve customers
            else
                
                # there are still some active customers below the maximum possible OE value 
                # check if there is still some spare capacity to be allocated
                tx_spare_capacity_temp = tx_spare_capacity_exp - sum(alloc_exp_p_temp) 
            end
        end    
        # Second main step: the P-V sensitivity curve for the critical customer is used to estimate its voltage
        flag = 0 # flag to indicate when there is a voltage problem (0 = problems; 1 = no problems)
        
        # select the expected active power export for the critical customer, which is to use the full 
        # allocation of spare capacity made above 
        expected_pnet_crit_cust = alloc_exp_p_temp[12] # 12 is the location of the critical customer in the array
        
        # run the loop until no voltage problem is expected
        while flag == 0
            
            # use the P-V sensitivity curve to estimate the voltage at the critical customer for the expected 
            # active power export
            expected_v_crit_cust = p_v_sensit_curve(expected_pnet_crit_cust)
            
            # check if the expected voltage at the critical customer is above 253V  
            if expected_v_crit_cust > 253
                expected_pnet_crit_cust -= 0.5 # if above the limit, reduce the allocated OE by 0.5kW (pre-defined reduction step)
                
                # check if the OE has achieved the minimum possible OE value (zero)
                if expected_pnet_crit_cust < 0
                    expected_pnet_crit_cust = 0 # if below zero, allocate zero as OE value
                    
                    # voltage problems are still expected even after OE value is equal to zero. OE calculation finished.
                    break
                end
            else
                flag = 1 # no voltage problems are expected. OE calculation finished.
            end
        end
    end
    ac_crv_oe_exp_values = expected_pnet_crit_cust # rename the allocated OE value
    
    return ac_crv_oe_exp_values # return the calculated OE value for the export for the considered time step
end

# Run the AC & CV OE Algs for the Day

# define number of time steps in the day (24h) for the corresponding time resolution
num_of_time_step = 288

# initialise arrays to save the OE values for each time step of the day
ac_crv_oe_day_exp_values = zeros(DFs.nrow(names_active_cust), num_of_time_step) # for exports

# Calculate the OE values for each time step of the day in a for loop
for iTime in range(1,num_of_time_step)
    tx_spare_capacity_exp = forecast_spare_capacity_exp[!,1][iTime] # separate the transformer capacity export for the current time step
    # call the Asset Capacity OE functions to calculate the OE value to each time step
    ac_crv_oe_day_exp_values[:,iTime] .= ac_crv_oe_exp_algorithm(names_active_cust, tx_spare_capacity_exp, active_cust_der_size, p_v_sensit_curve) # export
end

active_exp_val_path = "/Users/abond/OEAlgs-Julia/data/ac_crit_v/active_exp_values.csv"
active = DFs.DataFrame(ac_crv_oe_day_exp_values,Symbol.(Vector(range(1,num_of_time_step))))

CSV.write(active_exp_val_path, active)

# Performance Check (do later)