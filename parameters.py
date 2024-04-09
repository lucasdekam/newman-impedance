from numpy import sqrt

NMAX = 10

# Estimate from Wikipedia https://en.wikipedia.org/wiki/Conductivity_(electrolytic):
CONDUCTIVITY_THEOR = (34.982 + 6.80) * 1e-3 * 1e-4 * 1e3  # S/m
CONDUCTIVITY_MEAS = 2.8e-3  # measured (S/m)

# measured capacitance at PZC (uF/cm^2) * conversion uF -> F * conversion /cm^2 -> /m^2
CAPACITANCE = 35 * 1e-6 * 1e4  # F/m^2

RADIUS = 1.5e-3  # disk radius in m

# Frequency range
FREQ_MIN_HZ = 0.5  # Hz
FREQ_MAX_HZ = 1000  # Hz

# Gouy-Chapman capacitance
E_0 = 1.602e-19  # C
N_0 = 1e-4 * 1e3 * 6.022e23  # ion number density in the bulk, /m^3
K_B = 1.38e-23  # J/K
T = 298  # K
EPSILON = 8.854e-12 * 80  # water permittivity, F/m
LAMBDA_DEBYE = sqrt(K_B * T * EPSILON / (2 * E_0**2 * N_0))
GC_MIN_CAPACITANCE = EPSILON / LAMBDA_DEBYE
