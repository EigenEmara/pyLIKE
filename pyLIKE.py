import warnings

import numpy as np
from CoolProp.CoolProp import PropsSI
from tabulate import tabulate


# Enum definitions for geometry type, for integral length scale estimation.
class GeometryType:
    INTERNAL_FLOW = 0   # Internal flow (Hydraulic diameter based).
    EXTERNAL_FLOW = 1   # External flow (Boundary layer thickness based).
    WIND_TUNNEL = 2     # Wind tunnel (Grid spacing based).


def reynolds(rho, v, D, mu):
    return (rho * v * D) / mu


# Estimates the size of integral eddies based on system and geometry type:
# internal flow, external flow or wind tunnel.
def integral_length(geometry_type, x_char):
    if geometry_type == GeometryType.INTERNAL_FLOW:
        dh = x_char
        return 0.07 * dh

    if geometry_type == GeometryType.EXTERNAL_FLOW:
        delta = x_char
        return 0.45 * delta

    if geometry_type == GeometryType.WIND_TUNNEL:
        s = x_char
        return 0.2 * s


# Estimates turbulence intensity based on Reynolds number
def turbulence_intensity(Re):
    # | System or general condition                                             | I           |
    # |-------------------------------------------------------------------------+-------------|
    # | Aerodynamic airfoils                                                    | 0.003       |
    # | Core flows                                                              | 0.02 - 0.05 |
    # | Pipes, ducts, simple internal flows at intermediate to high Re          | 0.02 - 0.12 |
    # | High-velocity in complex systems (turbomachinery, heat exchangers,
    #   baffles, high swirl devices, near the wall                              | 0.05 - 0.2  |
    # | Atmospheric boundary layer flows, gusting winds, tornadoes, hurricanes,
    #   oscillating boundaries                                                  | 0.3         |

    # Sandborn (1955), Measuerd on pipe center for fully developed turbulent flow.
    I_sandborn = 0.144 * Re ** (-0.146)

    # One-Eighth power law
    I_oeigth = 0.16 * Re ** (-1/8)

    # Russo & Basse (2016) and Basse (2017) , Experimental and CFD estimations.
    # Measured, pipe axis (recommended for smooth pipes).
    I_russo_meas_pipe_axis = 0.055 * Re ** (-0.0407)

    # Measured, pipe area (higher intensity near the wall, except for viscous BL).
    I_basse_meas_pipe_area = 0.227 * Re ** (-0.1)

    # Incompressible CFD, pipe area
    I_cfd_pipe_area = 0.14 * Re ** (-0.079)

    return [I_sandborn, I_oeigth, I_russo_meas_pipe_axis,
            I_basse_meas_pipe_area, I_cfd_pipe_area]


# Turbulence kinetic energy K
def tke(intensity, v):
    return (3.0 / 2.0) * (intensity * v) ** 2


# Turbulence dissipation rate eps
def eps(k, l):
    return 0.09 * ((k ** (3.0/2.0)) / l)


# Integral eddy time and velocity scales.
def integral_scales(l, k, e):
    t = 0.09 * (k / e)      # Integral eddy time
    u = l / t               # Integral eddy velocity
    return t, u


# Taylor eddy time and velocity scales.
def taylor_scales(k, e, nu):
    l = np.sqrt(10 * k * nu / e)    # Taylor eddy length
    t = np.sqrt(15 * nu / e)        # Taylor eddy time
    u = l / t                       # Taylor eddy velocity
    return l, t, u


# Kolmogorov eddy time and velocity scales.
def kolmogorov_scales(e, nu):
    l = np.power((nu**3 / e), 0.25)   # Kolmogorov eddy length
    t = np.sqrt(nu / e)               # Kolmogorov eddy time
    u = np.power(nu*e, 0.25)          # Kolmogorov eddy velocity
    return l, t, u


# Specific dissipation rate
def omega(k, e):
    return e / (0.09 * k)


# Turbulent kinematic viscosity
def nu_t(k, e):
    return 0.09 * k * k / e


# Peak turbulent velocity
def u_turb_max(v):
    return 5.0/4.0 * v


# Calculates skin friction Cf, using Schlichting formula, external smooth plate flow.
def cf_external(Re):
    if Re > 1e9:
        raise ValueError('Schlichting formula is valid for Re <= 1e9')
    return np.power((2 * np.log10(Re) - 0.65), -2.3)


# Calculates skin firction Cf, using Moody formula, internal smooth and rough flow.
# ks: surface roughness and D is the hydraulic diameter
def cf_internal(Re, D, ks=0.0):
    b = np.power((2e4 * ks / D) + (1e6 / Re), 1.0/3.0)
    return 0.0055 * (1 + b)


def wall_shear_stress(cf, rho, u_char):
    return 0.5 * cf * rho * (u_char ** 2)


# Calculates shear velocity U-star dimensionless number
def u_star(tw, rho):
    return np.sqrt(tw / rho)


def y(yplus, ustar, nu):
    return yplus * nu / ustar


def y_plus(y, ustar, nu):
    return y * ustar / nu


class TurbulentCase:
    def __init__(self, x_char, u_char, substance_name, T, P, geo_type=GeometryType.INTERNAL_FLOW):
        self.x = x_char
        self.u = u_char
        self.T = T
        self.P = P
        self.rho = PropsSI('D', 'T', T, 'P', P, substance_name)
        self.mu = PropsSI('V', 'T', T, 'P', P, substance_name)
        self.nu = self.mu / self.rho
        self.Re = reynolds(self.rho, self.u, self.x, self.mu)

        if self.Re < 2300:
            warnings.warn('Reynolds number less than 2300.')

        self.geometry_type = geo_type

        # LIKE Algorithm
        self.l = integral_length(geo_type, self.x)
        self.I = turbulence_intensity(self.Re)[1]  # Sandborn
        self.k = tke(self.I, self.u)
        self.e = eps(self.k, self.l)

        # Eddy scales
        self.it, self.iu = integral_scales(self.l, self.k, self.e)
        self.tl, self.tt, self.tu = taylor_scales(self.k, self.e, self.nu)
        self.kl, self.kt, self.ku = kolmogorov_scales(self.e, self.nu)

    def report(self):
        data = [
            ['Reynolds Number (Re)', self.Re],
            ['Turbulence intensity (I)', self.I],
            ['Turbulence kinetic energy (k) [m2/s2]', self.k],
            ['Dissipation rate (Epsilon) [m2/s3]', self.e],
            ['Integral eddy length scale [m]', self.l],
            ['Integral eddy time scale [s]', self.it],
            ['Integral eddy velocity scale [m/s]', self.iu],
            ['Taylor eddy length scale [m]', self.tl],
            ['Taylor eddy time scale [s]', self.tt],
            ['Taylor eddy velocity scale [m/s]', self.tu],
            ['Kolmogorov eddy length scale [m]', self.kl],
            ['Kolmogorov eddy time scale [s]', self.kt],
            ['Kolmogorov eddy velocity scale [m/s]', self.ku],
            ['Turbulence kinematic viscosity [m2/s]', nu_t(self.k, self.e)],
            ['Specific dissipation rate (Omega)', omega(self.k, self.e)],
            ['Peak turbulent velocity [m/s]', u_turb_max(self.u)],
        ]
        print(tabulate(data, headers=['Parameter', 'Value']))

    def y(self, yplus, ks=0.0):
        cf = 0.0
        if self.geometry_type == GeometryType.INTERNAL_FLOW:
            cf = cf_internal(self.Re, self.x, ks)
        elif self.geometry_type == GeometryType.EXTERNAL_FLOW:
            cf = cf_external(self.Re)
        else:
            raise NotImplementedError('Calculating skin friction coefficient for unsupported geometry type.')

        tw = wall_shear_stress(cf, self.rho, self.u)
        ustar = u_star(tw, self.rho)
        return (yplus * self.nu) / ustar


if __name__ == '__main__':
    d = 1.0
    v = 0.012
    T = 300
    P = 1.013e5

    case = TurbulentCase(d, v, 'Air', T, P)
    case.report()
    print(case.y(30))

