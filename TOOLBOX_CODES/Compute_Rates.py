import pandas as pd
import numpy as np
import os

"""
H2O-H2 State-to-State Thermal Rate Coefficients Calculator

This script calculates effective and thermal rate coefficients for H2O-H2 collisions 
using state-to-state rates.

References:
1. Joy, C., Bostan, D., Mandal, B., & Babikov, D. (2024). 
Rate coefficients for rotational state-to-state transitions in H2O+ H2 collisions as predicted by mixed quantum–classical theory. 
Astronomy & Astrophysics, 692, A229.

Author: Carolin
Date: January 2025
Version: 1.0
"""

def read_state_to_state_data(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    start_idx = 0
    for i, line in enumerate(lines):
        if 'T_10.0' in line:
            start_idx = i + 1
            break

    data = []
    for line in lines[start_idx:]:
        values = line.split()
        if len(values) == 29:
            try:
                row = [float(x.replace('E', 'e')) for x in values]
                data.append(row)
            except ValueError:
                continue

    data = np.array(data)
    transition_water = data[:, [0,1,2,4,5,6]]  # J1,KA1,KC1,J1p,KA1p,KC1p
    transition_h2 = data[:, [3,7]]  # J2,J2p
    
    temperature = [10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 100, 150, 200, 
                  300, 400, 500, 800, 1000, 1200, 1500, 2000]
    rates_dict = {temp: data[:, i+8] for i, temp in enumerate(temperature)}
    
    return temperature, transition_water, transition_h2, rates_dict

def compute_effective_rates(transition_water, transition_h2, rates_dict):
    effective_rates = {}
    
    for temp in rates_dict.keys():
        effective_rates[temp] = {}
        
        # Get unique transitions based on input format
        unique_transitions = set()
        for i in range(len(transition_water)):
            trans = (transition_water[i,0], transition_water[i,1], transition_water[i,2], 
                    transition_h2[i,0], transition_water[i,3], transition_water[i,4], 
                    transition_water[i,5])
            unique_transitions.add(trans)
        
        for trans in unique_transitions:
            j1, ka1, kc1, j2, j1p, ka1p, kc1p = trans
            # Sum over J2p for each unique transition
            mask = ((transition_water[:, 0] == j1) & 
                   (transition_water[:, 1] == ka1) & 
                   (transition_water[:, 2] == kc1) &
                   (transition_h2[:, 0] == j2) &
                   (transition_water[:, 3] == j1p) &
                   (transition_water[:, 4] == ka1p) &
                   (transition_water[:, 5] == kc1p))
            
            effective_rates[temp][trans] = np.sum(rates_dict[temp][mask])
    
    return effective_rates

def compute_thermal_rates(effective_rates, temperature):
    kb = 0.6950386692  # Boltzmann constant in cm^-1/K
    
    h2_levels = {
        0: 0.0,
        2: 354.35,
        4: 1168.78,
        6: 2414.76,
        8: 4051.73
    }
    
    thermal_rates = {}
    for temp in temperature:
        thermal_rates[temp] = {}
        kbt = kb * temp
        
        # Get unique water transitions (J1,KA1,KC1,J1p,KA1p,KC1p)
        water_transitions = set()
        for k in effective_rates[temp].keys():
            water_trans = (k[0], k[1], k[2], k[4], k[5], k[6])
            water_transitions.add(water_trans)
        
        # Calculate partition function
        zpart = sum((2 * j2 + 1) * np.exp(-energy / kbt) 
                for j2, energy in h2_levels.items())
        
        for h2o_trans in water_transitions:
            j1, ka1, kc1, j1p, ka1p, kc1p = h2o_trans
            sum_rates = 0
            
            last_available_rate = None
            
            # Average over initial J2 states
            for j2, energy in h2_levels.items():
                weight = (2 * j2 + 1) * np.exp(-energy / kbt)
                trans = (j1, ka1, kc1, j2, j1p, ka1p, kc1p)
                
                rate = effective_rates[temp].get(trans, 0.0)
                
                if rate > 0:
                    last_available_rate = rate
                elif last_available_rate is not None:
                    rate = last_available_rate
                
                sum_rates += weight * rate

            thermal_rate = sum_rates / zpart
            thermal_rates[temp][h2o_trans] = thermal_rate
    
    return thermal_rates

def write_effective_rates(effective_rates, temperature, filename):
    """Write effective rates to file."""
    with open(filename, 'w') as f:
        f.write("! EFFECTIVE RATE COEFFICIENTS (cm^3 s^-1)\n\n")
        f.write("Temperature (K):" + "".join(f"{temp:15.1f}" for temp in temperature) + "\n\n")
        f.write("J1  KA1 KC1 J1p KA1p KC1p J2\n")
        
        # Sort transitions for consistent output
        all_transitions = sorted(effective_rates[temperature[0]].keys())
        
        for trans in all_transitions:
            j1, ka1, kc1, j1p, ka1p, kc1p, j2 = trans
            rates = [effective_rates[temp][trans] for temp in temperature]
            
            f.write(f"{int(j1):3d}{int(ka1):4d}{int(kc1):4d}"
                   f"{int(j1p):4d}{int(ka1p):5d}{int(kc1p):4d}{int(j2):4d}")
            f.write("".join(f"{rate:15.6E}" for rate in rates) + "\n")


def write_thermal_rates(thermal_rates, temperature, filename):
    """Write thermal rates to file."""
    with open(filename, 'w') as f:
        f.write("! THERMAL RATE COEFFICIENTS (cm^3 s^-1)\n\n")
        f.write("Temperature (K):" + "".join(f"{temp:15.1f}" for temp in temperature) + "\n\n")
        f.write("J1  KA1 KC1 J1p KA1p KC1p\n")
        
        # Sort transitions for consistent output
        all_transitions = sorted(thermal_rates[temperature[0]].keys())
        
        for trans in all_transitions:
            j1, ka1, kc1, j1p, ka1p, kc1p = trans
            rates = [thermal_rates[temp][trans] for temp in temperature]
            
            f.write(f"{int(j1):3d}{int(ka1):4d}{int(kc1):4d}"
                   f"{int(j1p):4d}{int(ka1p):5d}{int(kc1p):4d}")
            f.write("".join(f"{rate:15.6E}" for rate in rates) + "\n")

def main():
    # Input file
    input_file = "MQCT_DATA/MQCT_data_paraH2O_paraH2.txt" 

    # Create output directory 
    os.makedirs("./Outputs", exist_ok=True)    

    # Output files
    output_effective = "Outputs/output_effective_rates.dat"
    output_thermal = "Outputs/output_thermal_rates.dat"
    
    # Read and process data
    temperature, transition_h2o, transition_h2, rates_dict = read_state_to_state_data(input_file)
    
    # Compute effective rates
    effective_rates = compute_effective_rates(transition_h2o, transition_h2, rates_dict)
    
    # Compute thermal rates
    thermal_rates = compute_thermal_rates(effective_rates, temperature)
    
    # Write results to files
    write_effective_rates(effective_rates, temperature, output_effective)
    write_thermal_rates(thermal_rates, temperature, output_thermal)

if __name__ == "__main__":
    main() 