using CSV
using DataFrames

const DFs = DataFrames

# Load necessary inputs

tx_rated_capacity = 500

names_active_cust = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/names_active_cust.csv") |> DFs.DataFrame # active customers in the network
active_cust_der_size = ones(DFs.nrow(names_active_cust)) * 10.0

# forecast_agg_net_p_act_cust = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/forecast_agg_net_p_act_cust.csv") |> DFs.DataFrame # load the forecasted aggregated active customer net demand
# forecast_agg_net_p_passive_cust = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/forecast_agg_net_p_passive_cust.csv") |> DFs.DataFrame # load the forecasted aggregated passive customer net demand

# forecast_p_tx = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/forecast_p_tx.csv") |> DFs.DataFrame # load the transformer active (P) power demand
# forecast_q_tx = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/forecast_q_tx.csv") |> DFs.DataFrame # load the transformer reactive (Q) power demand
# forecast_s_tx = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/forecast_s_tx.csv") |> DFs.DataFrame # load the transformer apparent (S) power demand

forecast_spare_capacity_exp = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/forecast_spare_capacity.csv") |> DFs.DataFrame # load the spare capacity of the distribution transformer for exports

tx_pri_angles_day = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/tx_pri_angles.csv") |> DFs.DataFrame # load voltage angles at the primary side of the distribution transformer
tx_pri_voltages_day = CSV.File("/Users/abond/OEAlgs-Julia/data/day_30-15/tx_pri_voltages.csv") |> DFs.DataFrame # load voltage magnitudes at the primary side of the distribution transformer

# Define AC OE Export Algorithm

function ac_oe_exp_algorithm(
    names_active_cust,
    tx_spare_capacity_exp,
    active_cust_der_size)
    """
    Function to calculate the Asset Capacity OE proportional allocation for exports in each time step.
    This is used in all time steps of the day.    
    :param names_active_cust: the list of active customers 
    :param tx_spare_capacity_exp: the distribution transformer spare export capacity of the considered time step
    :param active_cust_der_size: the size of active customers' DERs 
    :return ac_oe_exp_values: the calculated Asset Capacity OE value for proportional allocation and exports
    """
    DFs = DataFrames
    # initialisation of the operating envelope values 
    alloc_exp_p_temp = zeros(DFs.nrow(names_active_cust)) 
    
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
                if alloc_exp_p_temp[ild] > 10 # Usually maximum injection allowed by DNSPs per phase (10kW)
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
    end

    ac_oe_exp_values = alloc_exp_p_temp # rename the allocated OE value
    
    return ac_oe_exp_values # return the calculated OE value for the export for the considered time step
end

# Run the AC OE Algs for the Day

# define number of time steps in the day (24h) for the corresponding time resolution
num_of_time_step = 288

# initialise arrays to save the OE values for each time step of the day
ac_oe_day_exp_values = zeros(DFs.nrow(names_active_cust), num_of_time_step) # for exports

# Calculate the OE values for each time step of the day in a for loop
for iTime in range(1,num_of_time_step)
    tx_spare_capacity_exp = forecast_spare_capacity_exp[!,1][iTime] # separate the transformer capacity export for the current time step
    # println(tx_spare_capacity_exp)
    # call the Asset Capacity OE functions to calculate the OE value to each time step
    ac_oe_day_exp_values[:,iTime] = ac_oe_exp_algorithm(names_active_cust, tx_spare_capacity_exp, active_cust_der_size) # export
end

# save the OE value for later    
if tx_rated_capacity == 500
    ac_oe_day_exp_values_original = ac_oe_day_exp_values
end

active_exp_val_path = "/Users/abond/OEAlgs-Julia/data/asset_cap/active_exp_values.csv"
active = DFs.DataFrame(ac_oe_day_exp_values,Symbol.(Vector(range(1,num_of_time_step))))

CSV.write(active_exp_val_path, active)

# Performance Check (do later)