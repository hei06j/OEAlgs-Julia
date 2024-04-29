#Proportional
function ideal_oe_prop_exp_algorithm(
    names_active_cust, 
    active_cust_der_size,
    dss_engine)
    """
    Function to calculate the Ideal OE proportional allocation for exports in each time step.
    This is used in all time steps of the day.    
    :param names_active_cust: the list of active customers (Julia DataFrame)
    :param active_cust_der_size: the size of active customers' DERs 
    :param dss_engine: the OpenDSS engine
    :return ideal_oe_prop_exp_values: the calculated ideal OE value for proportional allocation and exports
    :return lv_tx_util: the transformer utilisation when using the calculated OE value
    :return lv_hof_util_max: the maximum utilisation of the LV head of feeder when using the calculated OE value
    :return volt_all_cust_temp: voltages on all customers when using the calculated OE value 
    """
    DSS = dss_engine

    # load LV network data from the OpenDSS model
    load_list = DSS.Loads.AllNames() # list of loads
    line_list = DSS.Lines.AllNames() # list of lines

    ### Check whether voltage limits have been breached
    
    max_voltage = 253
    oe_step_size = 0.5

    # must initialise variables outside of loop for variables to be stored
    flag_thermal_line = 0
    flag_voltage = 0
    flag_thermal_tx = 0
    lv_tx_util = 0 
    lv_hof_util_max = 0 
    volt_all_cust_temp = [] 

    # lv_hof_util_f0 = 0
    # lv_hof_util_f1 = 0
    # lv_hof_util_f2 = 0

    alloc_exp_p_temp = active_cust_der_size[1] # allocated export power (temporarily set to a static 10 kw / can change according to different DER sizes)

    # iteratively run power flow calculations with reducing OE values until the point that no network limit (voltages or 
    # thermal) in any part of the LV network is breached.
    flag = 0 # flag to indicate when there is a technical problem (0 = problems; 1 = no problems)

    while flag==0

        for ild in range(1,nrow(names_active_cust))
            DSS.Circuit.SetActiveElement("load."*names_active_cust[!,1][ild]) # select an active customer
            temp = imag(DSS.CktElement.Powers()[1]) # save the current reactive power of the active customer
            DSS.Properties.Value("status","fixed") # fix the load status to run with the set value instead of profile
            DSS.Properties.Value("kW",string(-1*alloc_exp_p_temp)) # set the kW to the allocated OE value
            DSS.Properties.Value("kvar",string(temp)) # set the kvar to be the same as before
        end
        
        DSS.Solution.SolveSnap() # solve the power flow in OpenDSS without changing the time of the day

        # verify voltages on all customers
        
        # initialisation of variable to save voltage on all customers for the current OE value
        volt_all_cust_temp = [] 
        
        # collect the voltage from all customers
        for ild in range(1,length(load_list))
            DSS.Circuit.SetActiveElement("load."*string(load_list[ild])) # select a customer
            println(DSS.CktElement.VoltagesMagAng()[1])
            push!(volt_all_cust_temp, DSS.CktElement.VoltagesMagAng()[1]) # extract its voltage magnitude
        end
        
        # verify if the voltage of any of the customers is above the maximum statutory limit
        if maximum(volt_all_cust_temp) > max_voltage
            flag_voltage = 0 # flag to indicate when there is a voltage problem (0 = problems; 1 = no problems)
        else
            flag_voltage = 1 # flag to indicate when there is a voltage problem (0 = problems; 1 = no problems)
        end
        
        # initialise variables
        lv_hof_util_f0 = 0
        lv_hof_util_f1 = 0
        lv_hof_util_f2 = 0

        # verify each line for thermal issues
        for iline in range(1,length(line_list))
            DSS.Circuit.SetActiveElement("line."*string(line_list[iline])) # select a line

            currents = DSS.CktElement.CurrentsMagAng()
            
            # extract the current passing through the lines in phases A, B, and C
            Ia_line = maximum((currents[1],currents[7])) 
            Ib_line = maximum((currents[3],currents[9]))
            Ic_line = maximum((currents[5],currents[11]))
            
            # select the corresponding line
            DSS.Lines.Name(string(line_list[iline]))

            # extract the rated current of the corresponding line
            I_rated = DSS.Lines.NormAmps()
            
            # find the maximum per phase current passing through the line
            I_max_line_temp = maximum((Ia_line, Ib_line, Ic_line)) 
            
            # check if the line current is above the rated current of the corresponding line
            if I_max_line_temp > I_rated
                flag_thermal_line = 0 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
                break # in case any line has a thermal problem, the for loop stops (no need to check other lines) 
            else
                flag_thermal_line = 1 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
                # calculate the utilisation of the head of each feeder (there are 3 feeders in this LV network)
                if line_list[iline] == "hv_f0_lv28_f0_l0"
                    lv_hof_util_f0 = 100 * I_max_line_temp / I_rated
                elseif line_list[iline] == "hv_f0_lv28_f1_l0"
                    lv_hof_util_f1 = 100 * I_max_line_temp / I_rated
                elseif line_list[iline] == "hv_f0_lv28_f2_l0"
                    lv_hof_util_f2 = 100 * I_max_line_temp / I_rated
                end
            end
        end
        
        # keep the maximum utilisation on the head of feeder
        lv_hof_util_max = maximum((lv_hof_util_f0, lv_hof_util_f1, lv_hof_util_f2))

        # verify thermal limits of transformer
        
        # select the distribution transformer and extract rated voltage at the secondary
        DSS.Circuit.SetActiveElement("transformer.hv_f0_lv28_tx") 
        DSS.Transformers.Wdg(2.0) # select secondary winding (LV side)
        tx_rated_volt_sec = DSS.Transformers.kV() / sqrt(3) # Extract the transformer rated voltage at the secondary side
        
        # extract the distribution transformer rated power capacity
        tx_rated_power = parse(Float64, split(strip(DSS.Properties.Value("kVAs"), '['),", ")[1]) # why do we need this? never used

        # calculated distribution transformer current capacity (Amperes) per phase
        tx_amp_capacity_phase = (parse(Float64, split(strip(DSS.Properties.Value("kVAs"), '['),", ")[1]) 
                                / 3) / tx_rated_volt_sec 

        # extract the current passing through the distribution transformer for phases A, B, and C
        currents = DSS.CktElement.CurrentsMagAng()
        Ia_tx = currents[9]
        Ib_tx = currents[11]
        Ic_tx = currents[13]

        # check if the current passing through the distribution transformer is above the calculated current capacity
        if (Ia_tx + Ib_tx + Ic_tx) > (3 * tx_amp_capacity_phase)
            flag_thermal_tx = 0 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
        else
            flag_thermal_tx = 1 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
            lv_tx_util = 100 * (Ia_tx + Ib_tx + Ic_tx) / (3 * tx_amp_capacity_phase) # calculate the utilisation on the transformer
            # println(lv_tx_util)
        end

        # check if any technical isssue was found in the network           
        flag = flag_thermal_line * flag_voltage * flag_thermal_tx # flag to indicate when there is a technical problem (0 = problems; 1 = no problems)
        if flag == 0
            alloc_exp_p_temp -= oe_step_size # if above any limit, reduce the allocated OE by 0.5kW (pre-defined reduction step)
            
            # check if the OE has achieved the minimum possible OE value (zero)
            if alloc_exp_p_temp < 0
                alloc_exp_p_temp = 0 # if below zero, allocate zero as OE value
                break # voltage problems are still expected even after OE value is equal to zero. OE calculation finished.
            end
        end
    end

    ideal_oe_prop_exp_values = alloc_exp_p_temp # rename the allocated OE value

    return ideal_oe_prop_exp_values, lv_tx_util, lv_hof_util_max, volt_all_cust_temp # return the calculated OE value for the export for the considered time step
end

#Max Allocation
function ideal_oe_max_exp_algorithm(
    names_active_cust,
    active_cust_der_size,
    dss_engine = OpenDSSDirect,
    data_frames = DataFrames)
    """
    Function to calculate the Ideal OE maximum allocation for exports in each time step.
    This is used in all time steps of the day.    
    :param names_active_cust: the list of active customers 
    :param active_cust_der_size: the size of active customers' DERs 
    :param dss_engine: the OpenDSS engine
    :return ideal_oe_max_exp_values: the calculated ideal OE value for maximum allocation and exports
    :return lv_tx_util: the transformer utilisation when using the calculated OE value
    :return lv_hof_util_max: the maximum utilisation of the LV head of feeder when using the calculated OE value
    :return volt_all_cust_temp: voltages on all customers when using the calculated OE value
    """
    
    DSS = dss_engine
    DFs = data_frames

    # load LV network data from the OpenDSS model
    load_list = DSS.Loads.AllNames() # list of loads
    line_list = DSS.Lines.AllNames() # list of lines

    ### Check whether voltage limits have been breached
    
    max_voltage = 253
    oe_step_size = 0.5

    # must initialise variables outside of loop for variables to be stored
    flag_thermal_line = 0
    flag_voltage = 0
    flag_thermal_tx = 0
    lv_tx_util = -999
    lv_hof_util_max = -999 
    volt_all_cust_temp = [] 

    # changing this to a list will make this a true proportional allocation algorithm
    alloc_exp_p_temp = ones(DFs.nrow(names_active_cust))*10.0 # set to 10 kW for now

    #### INVESTIGATE AND FIND OUT WHY active_cust_der_size KEPT CHANGING!!!! WEIRD

    ### PLOT AND SEE IF SAME AS PYTHON NOTEBOOK
    # alloc_exp_p_temp = active_cust_der_size

    # println(alloc_exp_p_temp)

    # iteratively run power flow calculations with reducing OE values until the point that no network limit (voltages or 
    # thermal) in any part of the LV network is breached.    
    flag = 0 # flag to indicate when there is a technical problem (0 = problems; 1 = no problems)
    while flag == 0
        
        for ild in range(1,DFs.nrow(names_active_cust))
            DSS.Circuit.SetActiveElement("load."*names_active_cust[!,1][ild]) # select an active customer
            temp = imag(DSS.CktElement.Powers()[1]) # save the current reactive power of the active customer
            DSS.Properties.Value("status","fixed") # fix the load status to run with the set value instead of profile
            DSS.Properties.Value("kW",string(-1*alloc_exp_p_temp[ild])) # set the kW to the allocated OE value
            DSS.Properties.Value("kvar",string(temp)) # set the kvar to be the same as before
            # println(temp)
        end  # solve the power flow without changing the time of the day

        DSS.Solution.SolveSnap()

        # verify thermal limits of lines

        # initialise variables
        lv_hof_util_f0 = 0
        lv_hof_util_f1 = 0
        lv_hof_util_f2 = 0        
        
        # verify each line for thermal issues
        for iline in range(1,length(line_list))
            DSS.Circuit.SetActiveElement("line."*string(line_list[iline])) # select a line
            
            currents = DSS.CktElement.CurrentsMagAng()
            # extract the current passing through the lines in phases A, B, and C
            Ia_line = currents[1]
            Ib_line = currents[3]
            Ic_line = currents[5]

            DSS.Lines.Name(string(line_list[iline]))
            
            # extract the rated current of the corresponding line
            I_rated = DSS.Lines.NormAmps()
            
            # find the maximum per phase current passing through the line
            I_max_line_temp = maximum((Ia_line, Ib_line, Ic_line))
            
            # check if the line current is above the rated current of the corresponding line
            if I_max_line_temp > I_rated
                flag_thermal_line = 0 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
                break # in case any line has a thermal problem, the for loop stops (no need to check other lines) 
            else
                flag_thermal_line = 1 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
                
                # calculate the utilisation of the head of each feeder (there are 3 feeders in this LV network)
                if line_list[iline] == "hv_f0_lv28_f0_l0"
                    lv_hof_util_f0 = 100 * I_max_line_temp / I_rated
                elseif line_list[iline] == "hv_f0_lv28_f1_l0"
                    lv_hof_util_f1 = 100 * I_max_line_temp / I_rated
                elseif line_list[iline] == "hv_f0_lv28_f2_l0"
                    lv_hof_util_f2 = 100 * I_max_line_temp / I_rated
                end
            end
        end
        
        # keep the maximum utilisation on the head of feeder
        lv_hof_util_max = maximum((lv_hof_util_f0, lv_hof_util_f1, lv_hof_util_f2))
            
        # check if a thermal problem was found and reduce the OE value if needed
        if flag_thermal_line == 0 # if True, there are thermal problems
            alloc_exp_p_temp .-= oe_step_size # reduce the OE value of all active customers by the pre-defined step (0.5kW)
            
            # check if the OE has achieved the minimum possible OE value (zero) for all active customers
            if all(i->(i==0),alloc_exp_p_temp)
                break # thermal problems are still expected even after OE value is equal to zero. OE calculation finished.
            # check if the OE is below the minimum possible OE value (zero)
            elseif any(i->(i<0),alloc_exp_p_temp)
                # if below zero, allocate zero as OE value
                ifelse.(alloc_exp_p_temp .< 0, 0, alloc_exp_p_temp)
                continue
            end
            continue
        end
 
        # verify thermal limits of transformer
 
        # select the distribution transformer and extract rated voltage at the secondary
        DSS.Circuit.SetActiveElement("transformer.hv_f0_lv28_tx") 
        DSS.Transformers.Wdg(2.0) # select secondary winding (LV side)
        tx_rated_volt_sec = DSS.Transformers.kV() / sqrt(3)

        # extract the distribution transformer rated power capacity
        # tx_rated_power = parse(Float64, split(strip(DSS.Properties.Value("kVAs"), '['),", ")[1]) # why do we need this? never used

        # calculated distribution transformer current capacity (Amperes) per phase
        tx_amp_capacity_phase = (parse(Float64, split(strip(DSS.Properties.Value("kVAs"), '['),", ")[1]) 
                                / 3) / tx_rated_volt_sec 

        # extract the current passing through the distribution transformer for phases A, B, and C
        currents = DSS.CktElement.CurrentsMagAng()
        Ia_tx = currents[9]
        Ib_tx = currents[11]
        Ic_tx = currents[13]
               
        # check if the current passing through the distribution transformer is above the calculated current capacity
        if (Ia_tx + Ib_tx + Ic_tx) > (3 * tx_amp_capacity_phase)
            flag_thermal_tx = 0 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
        else
            flag_thermal_tx = 1 # flag to indicate when there is a thermal problem (0 = problems; 1 = no problems)
            lv_tx_util = 100 * (Ia_tx + Ib_tx + Ic_tx) / (3 * tx_amp_capacity_phase) # calculate the utilisation on the transformer
        end
        # check if a thermal problem was found and reduce the OE value if needed
        if flag_thermal_tx == 0 # if True, there are thermal problems
            alloc_exp_p_temp .-= oe_step_size # reduce the OE value of all active customers by the pre-defined step (0.5kW)
            
            # check if the OE has achieved the minimum possible OE value (zero) for all active customers
            if all(i->(i==0),alloc_exp_p_temp)
                break # thermal problems are still expected even after OE value is equal to zero. OE calculation finished. 
            # check if the OE is below the minimum possible OE value (zero)
            elseif any(i->(i<0),alloc_exp_p_temp)
                # if below zero, allocate zero as OE value
                ifelse.(alloc_exp_p_temp .< 0, 0, alloc_exp_p_temp)
                continue
            end
            continue
        end
            
        # verify voltages on all customers
        
        # initialisation of variable to save voltage on all customers for the current OE value
        volt_all_cust_temp = [] 
                
        # collect the voltage from all customers
        for ild in range(1,length(load_list))
            DSS.Circuit.SetActiveElement("load."*string(load_list[ild])) # select a customer
            push!(volt_all_cust_temp, DSS.CktElement.VoltagesMagAng()[1]) # extract its voltage magnitude
        end

        # check if a there is any voltage problem and reduce the OE value if needed
        if maximum(volt_all_cust_temp) > max_voltage
            
            # check if the OE has achieved the minimum possible OE value (zero) for all active customers
            if all(i->(i==0),alloc_exp_p_temp)
                break # voltage problems are still expected even after OE value is equal to zero. OE calculation finished.
            else
                alloc_sort_temp = sort(volt_all_cust_temp) # sort voltages
                loc_sort_temp = 0
                cust_size = length(alloc_sort_temp)
                flag_update = 0 # flag to indicate when there is a technical problem (0 = problems; 1 = no problems)
                while flag_update == 0
                    max_voltage_temp = alloc_sort_temp[cust_size-loc_sort_temp]
                    # indices ensures that all of the customers are caught (previous method might miss an active cust)
                    indices = [i for (i, v) in enumerate(volt_all_cust_temp) if v == max_voltage_temp]
                    
                    for loc_max_temp in indices
                        # find index of maximum voltage
                        cust_temp = load_list[loc_max_temp]
                        # check if the maximum voltage is happening is a active customer
                        if cust_temp in names_active_cust[:,1]
                            loc_active_temp = findfirst(==(cust_temp),names_active_cust[:,1])
                            
                            # check if the OE value is already the minimum possible value (zero)
                            if alloc_exp_p_temp[loc_active_temp] != 0
                                alloc_exp_p_temp[loc_active_temp] -= oe_step_size # reduce the OE value of for the active customers with voltage probelms by the pre-defined step (0.5kW)
                                flag_voltage = 0 # update flag 
                                flag_update = 1 # update flag
                                
                                # check if the OE is below the minimum possible OE value (zero)
                                if any(i->(i<0),alloc_exp_p_temp)
                                    # if below zero, allocate zero as OE value
                                    ifelse.(alloc_exp_p_temp .< 0, 0, alloc_exp_p_temp)
                                end
                            # check if the OE has achieved the minimum possible OE value (zero) for all active customers
                            elseif all(i->(i==0),alloc_exp_p_temp)
                                break # voltage problems are still expected even after OE value is equal to zero. OE calculation finished.
                            end
                        end
                    end
                    loc_sort_temp += 1
                end
            end
        else
            flag_voltage = 1 # update flag
        end

        flag = flag_thermal_line * flag_thermal_tx * flag_voltage # update general flag to exit the while
    end

    ideal_oe_max_exp_values = alloc_exp_p_temp # rename the allocated OE value
    
    return ideal_oe_max_exp_values, lv_tx_util, lv_hof_util_max, volt_all_cust_temp # return the calculated OE value for the export for the considered time step
end

#Asset Capacity
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

#Asset Capacity Critical Voltage
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

#Asset Capacity and Delta V
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

#Produce Output
function run_oe_algs(
    ideal=True,
    asset_cap=True,
    ac_crv=True,
    ac_deltav=True)
    """
    Function to run all oe algs
    """
    #initialise inputs

    #run all
    if ideal
        ideal_oe_prop_exp_algorithm()
        ideal_oe_max_exp_algorithm()
    end
    if asset_cap
        ac_oe_exp_algorithm()
    end
    if ac_crv
        ac_crv_oe_exp_algorithm()
    end
    if ac_deltav
        ac_deltav_oe_exp_algorithm()
    end
end

#change inverter ctrl and see differences -- what are good metrics to compare?
