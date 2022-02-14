import numpy as np

def M_conservative_upsample(M, dp):
    """
    Upsample a distribution with total mass and number concentrations
    conservation while keeping the mean mass diameter equal to the geometric
    mean diameter, assuming logarithmic spacing of section boundaries.

    M = N * pi/6 * rho * (dmin*dmax)**(3/2)
    and
    M_i = N_i * pi/6 * rho * (d_{i-1}*d_i)**(3/2)

    > Input :
    M : float, total mass concentration
    dp : float array, section boundaries

    > Output :
    [M_i] float array, mass concentration in successive sections
    [N_i] float array, number concentration in successive sections
    """

    Mi = M * (dp[1:]**1.5 - dp[:-1]**1.5) / (dp[-1]**1.5 - dp[0]**1.5)

    return Mi

def M_N_conservative_upsample(M, rho, dp):
    """
    Upsample a distribution with total mass and number concentrations
    conservation while keeping the mean mass diameter equal to the geometric
    mean diameter, assuming logarithmic spacing of section boundaries.

    M = N * pi/6 * rho * (dmin*dmax)**(3/2)
    and
    M_i = N_i * pi/6 * rho * (d_{i-1}*d_i)**(3/2)

    > Input :
    M : float, total mass concentration
    rho : float, density
    dp : float array, section boundaries

    > Output :
    [Mi] array of mass concentration in successive sections
    [Ni] array of number concentration in successive sections
    """

    Mi = M_conservative_upsample(M, dp)

    Ni = Mi / (rho*np.pi/6) * (dp[1:] * dp[:-1])**(-1.5)

    return Mi, Ni

def M_conservative_upsample_logspace(M, dmin, dmax, ns):
    """
    Upsample a distribution with total mass and number concentrations
    conservation while keeping the mean mass diameter equal to the geometric
    mean diameter, assuming logarithmic spacing of section boundaries.

    M = N * pi/6 * rho * (dmin*dmax)**(3/2)
    and
    M_i = N_i * pi/6 * rho * (d_{i-1}*d_i)**(3/2)

    > Input :
    M : float, total mass concentration
    dmin : float, lowest boundary of diameter range
    dmax : float, higher boundary of diameter range
    ns : int, number of sections

    > Output :
    [Mi] array of mass concentration in successive sections
    [Ni] array of number concentration in successive sections
    """

    dp = np.geomspace(dmin, dmax, num=ns+1)

    return M_conservative_upsample(M, dp)

def M_N_conservative_upsample_logspace(M, rho, dmin, dmax, ns):
    """
    Upsample a distribution with total mass and number concentrations
    conservation while keeping the mean mass diameter equal to the geometric
    mean diameter, assuming logarithmic spacing of section boundaries.

    M = N * pi/6 * rho * (dmin*dmax)**(3/2)
    and
    M_i = N_i * pi/6 * rho * (d_{i-1}*d_i)**(3/2)

    > Input :
    M : float, total mass concentration
    N : float, total number concentration
    rho : float, density
    dmin : float, lowest boundary of diameter range
    dmax : float, higher boundary of diameter range
    ns : int, number of sections

    > Output :
    [Mi] array of mass concentration in successive sections
    [Ni] array of number concentration in successive sections
    """

    dp = np.geomspace(dmin, dmax, num=ns+1)

    return M_N_conservative_upsample(M, rho, dp)
