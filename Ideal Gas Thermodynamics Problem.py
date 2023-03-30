import math

# Constants
R = 8.314 # Gas constant, J/(mol*K)
gamma = 1.4 # Heat capacity ratio
M_A = 32 # Molecular mass of gas A, g/mol
M_B = 28 # Molecular mass of gas B, g/mol
M_C = 28 # Molecular mass of gas C, g/mol

# Initial conditions
V = 1 # Volume, m^3
T_A = 400 # Temperature of gas A, K
T_B = T_C = 300 # Temperature of gases B and C, K
P_A = 5e6 # Pressure of gas A, Pa
P_B = P_C = 1e6 # Pressure of gases B and C, Pa
V_A = V/2 # Volume of gas A, m^3
V_B = V_C = V/4 # Volume of gases B and C, m^3

# Final conditions
T_f = 350 # Final temperature, K

# Step 1
n_A = P_A*V_A/(R*T_A) # Number of moles of gas A
V_Af = n_A*R*T_A/P_A # Final volume of gas A
W1 = -n_A*R*T_A*math.log(V_Af/V_A) # Work done by gas A
Q1 = 0 # Heat transferred

# Step 2
T_Af = T_A*(V_A/V_Af)**(gamma-1) # Final temperature of gas A
W2 = -n_A*R*(T_Af-T_A)/(gamma-1) # Work done by gas A
Q2 = 0 # Heat transferred

# Step 3
T_Aff = T_Af*(P_A/P_A)**((gamma-1)/gamma) # Final temperature of gas A
W3 = -n_A*R*(T_Aff-T_Af)/(gamma-1) # Work done by gas A
Q3 = 0 # Heat transferred

# Step 4
n_B = P_B*V_B/(R*T_B) # Number of moles of gas B
n_C = P_C*V_C/(R*T_C) # Number of moles of gas C
n_BC = n_B + n_C # Total number of moles of gases B and C
T_BCf = (n_B*T_B + n_C*T_C)/n_BC # Final temperature of gases B and C
W4 = 0 # Work done
Q4 = 0 # Heat transferred

# Step 5
V_BCf = n_BC*R*T_BCf/(2*P_B) # Final volume of gases B and C
W5 = -n_BC*R*T_BCf/(gamma-1)*(V_BCf - n_BC*R*T_BCf/P_B) # Work done on gases B and C
Q5 = 0 # Heat transferred

# Final step
n_Bf = n_B*n_BC/V_BCf # Final number of moles of gas B
n_Cf = n_C*n_BC/V_BCf # Final number of moles of gas C
n_Af = n_A # Final number of moles of gas A
V_Bf = V_Cf = V_BCf/2 # Final volume of gases B and C
P_Bf = P_Cf = n_BC*R*T_BCf/V_BCf # Final pressure of gases B and C
T_Af = T_Aff = T_Bf = T_Cf = T_f # Final temperatures of all gases

# Entropy changes
S_A = n_A*R*math.log(V_Af/V_A)
S_B = n_B*R*math.log(V_Bf/V_B) + n_B*R*math.log(P_B/P_Bf)
S_C = n_C*R*math.log(V_Cf/V_C) + n_C*R*math.log(P_C/P_Cf)

# Printing results
print("Final conditions:")
print(f"Number of moles of gas A: {n_Af:.4f}")
print(f"Number of moles of gas B: {n_Bf:.4f}")
print(f"Number of moles of gas C: {n_Cf:.4f}")
print(f"Volume of gases B and C: {V_BCf:.4f} m^3")
print(f"Volume of gas A: {V_Af:.4f} m^3")
print(f"Pressure of gases B and C: {P_Bf:.4f} Pa")
print(f"Temperature of all gases: {T_f} K")
print()
print("Entropy changes:")
print(f"Change in entropy of gas A: {S_A:.4f} J/K")
print(f"Change in entropy of gas B: {S_B:.4f} J/K")
print(f"Change in entropy of gas C: {S_C:.4f} J/K")
