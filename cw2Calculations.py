import numpy as np

airRI = 1
SiRI = 3.46
SiO2RI = 1.45
GaAsRI = 3.5
A1AsRI = 3.0
Na3A1F6RI = 1.34
ZnSeRI = 2.5

a = 50.0
b = 173.0

def calc_r (n_a, n_b):
    return np.abs((n_a - n_b)/(n_a + n_b))

def calc_g (n_a, n_b, a, b):
    numerator = (n_a * a) - (n_b * b)
    denominator = (n_a * a) + (n_b * b)
    return np.abs(numerator/denominator)

def calc_omega_0 (c, n_a, n_b, a, b):
    numerator = np.pi * c
    denominator = (n_a * a) + (n_b * b)
    return numerator / denominator

# for when g = 0
def calc_approx_width (r):
    return (4 * r) / np.pi

# for non zero g
def approx_width_nonzero (r, g):
    numerator = 4 * r * np.sqrt(1 - g**2 * ((np.pi ** 2 / 4) - r**2))
    denominator = np.pi * (1 + r**2 * g**2)
    return numerator / denominator

# wrong
def calc_exact_width (r):
    return np.arccos(((2 * r **2) - 1)/np.pi)

r = calc_r(SiRI, airRI)
g = calc_g(SiRI, airRI, a, b)
omega_0 = calc_omega_0(3E8, SiRI, airRI, a, b)
approx_width = calc_approx_width(r)
omega_0_wavelength = 3E8 / omega_0
approx_width_wavelength = approx_width * omega_0_wavelength
exact_width = calc_exact_width(r)

print('r is', r)
print('g is', g)
print('omega_0 is', omega_0)
print('approx width is', approx_width)
print('approx width is', approx_width_wavelength)
print('exact width is', exact_width)