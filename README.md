# pyLIKE
 Implementation of LIKE algorithm by Dr. Sal Rodriguez to estimate CFD turbulence case initial conditions

## Example Case
Air flow through a pipe with diamter = 0.5 m, velocity = 2.0 m/s, T = 300 K at atmoshperic pressure.


    if __name__ == '__main__':
        d = 0.5
        v = 2.0
        T = 300
        P = 1.013e5

        case = TurbulentCase(d, v, 'Air', T, P)
        case.report()

## Results

    Parameter                                       Value
    -------------------------------------  --------------
    Reynolds Number (Re)                       63477.6
    Turbulence intensity (I)                   0.0401599
    Turbulence kinetic energy (k) [m2/s2]      0.0096769
    Dissipation rate (Epsilon) [m2/s3]         0.00244782
    Integral eddy length scale [m]             0.035
    Integral eddy time scale [s]               0.355795
    Integral eddy velocity scale [m/s]         0.0983712
    Taylor eddy length scale [m]               0.0249556
    Taylor eddy time scale [s]                 0.310703
    Taylor eddy velocity scale [m/s]           0.0803198
    Kolmogorov eddy length scale [m]           0.00112419
    Kolmogorov eddy time scale [s]             0.0802233
    Kolmogorov eddy velocity scale [m/s]       0.0140133
    Turbulence kinematic viscosity [m2/s]      0.00344299
    Specific dissipation rate (Omega)          2.81061
    Peak turbulent velocity [m/s]              2.5
