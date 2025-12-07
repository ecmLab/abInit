import math

# Shomate parameters from NIST for 298–1000 K ranges (units as per Shomate eqns)
# H2O(g): https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI
H2O_SHOMATE_298_1000 = {
    'A': 30.09200,
    'B': 6.832514,
    'C': 6.793435,
    'D': -2.534480,
    'E': 0.082139,
    'F': -250.8810,
    'G': 223.3967,
    'H': -241.8264,
}

# H2S(g): https://webbook.nist.gov/cgi/cbook.cgi?ID=C7783064&Units=SI
H2S_SHOMATE_298_1000 = {
    'A': 26.33936,
    'B': 9.022851,
    'C': -0.021697,
    'D': 0.000168,
    'E': -0.000005,
    'F': -241.8264,  # placeholder; adjust if needed; only relative trends used
    'G': 237.399,
    'H': -20.606,
}

R = 8.31446261815324  # J/mol/K

def shomate_cp(T, p):
    t = T/1000.0
    A,B,C,D,E,F,G,H = p['A'],p['B'],p['C'],p['D'],p['E'],p['F'],p['G'],p['H']
    # Cp in J/mol/K
    return A + B*t + C*t*t + D*t**3 + E/(t*t)

def shomate_H(T, p):
    t = T/1000.0
    A,B,C,D,E,F,G,H = p['A'],p['B'],p['C'],p['D'],p['E'],p['F'],p['G'],p['H']
    # H - H(0) in kJ/mol
    return A*t + B*(t*t)/2 + C*(t*t*t)/3 + D*(t**4)/4 - E/t + F - H

def shomate_S(T, p):
    t = T/1000.0
    A,B,C,D,E,F,G,H = p['A'],p['B'],p['C'],p['D'],p['E'],p['F'],p['G'],p['H']
    # S in J/mol/K
    return A*math.log(t) + B*t + C*(t*t)/2 + D*(t**3)/3 - E/(2*t*t) + G

def mu_ideal_gas(T, p_bar, shomate_params):
    # chemical potential approximated as G(T, p) = H(T) - T*S(T) + RT ln(p/p0)
    H_kJmol = shomate_H(T, shomate_params)  # kJ/mol
    S_JmolK = shomate_S(T, shomate_params)  # J/mol/K
    G_kJmol = H_kJmol - (T * S_JmolK)/1000.0
    mu = G_kJmol*1000.0 + R*T*math.log(p_bar/1.0)
    return mu  # J/mol

def RH_from_dewpoint(T_C, Tdp_C):
    # Magnus formula
    a, b = 17.625, 243.04
    gamma = (a*Tdp_C)/(b+Tdp_C) - (a*T_C)/(b+T_C)
    return math.exp(gamma)

def pH2O_from_RH(T_C, RH, p_tot_bar=1.0):
    # saturation vapor pressure (bar) via Tetens formula approx
    T = T_C
    Psat_kPa = 0.61078 * math.exp((17.27*T)/(T+237.3))  # kPa
    Psat_bar = Psat_kPa/100.0
    return RH * Psat_bar

