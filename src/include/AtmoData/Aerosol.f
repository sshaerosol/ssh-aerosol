C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
C-----------------------------------------------------------------------

c     File: Aerosol.f

c     Function: compute_gas_diffusivity
c
c     Computes gas-phase diffusivity in air.
c
c     Parameters:
c     temperature - Temperature (K).
c     pressure - Atmospheric pressure (Pa).
c     diameter - Molecular diameter (Angstrom).*
c     weight - Molecular weight (g/mol).
c     collision - Collision factor.
c
c     Returns:
c     diffusivity - Gas-phase diffusivity (m^2/s).
      subroutine ssh_compute_gas_diffusivity(temperature, pressure,
     $     diameter, weight,collision, diffusivity)

      double precision temperature
      double precision pressure
      double precision diameter
      double precision weight
      double precision collision

      double precision diffusivity

      double precision x
      double precision collision_integral
      double precision sigma

      x = temperature / dsqrt(collision * 9.7d1)
      call ssh_compute_collision_integral(x, collision_integral)

      sigma = .5d0 * (diameter + 3.617d0)

      diffusivity =
     $     1.8829225d-2 / pressure * temperature * dsqrt(temperature)
     $     * dsqrt((weight + 2.897d1) / (weight * 2.897d1))
     $     / (sigma * sigma * collision_integral)

      end


c     Function: gerber_wet_diameter
c
c     Computes wet diameter based on Gerber's formula.
c
c     Parameters:
c     relative_humidity - Relative humidity.
c     temperature - Temperature (K).
c     dry_diameter - Aerosol dry diameter (micrometer).
c
c     Returns:
c     wet_diameter - Aerosol wet diameter (micrometer).
      subroutine ssh_gerber_wet_diameter(relative_humidity, temperature,
     $     dry_diameter, wet_diameter)

      double precision relative_humidity
      double precision temperature
      double precision dry_diameter

      double precision wet_diameter

      double precision dry_radius
      double precision correction

c     Computes the dry radius in centimeters.
      dry_radius = dry_diameter * 0.5d-4
c     Temperature correction.
      correction = 0.5372215625062934D-12
     $     * (1.d0 + 0.3942463621284677d-02 * (298.d0 - temperature))

c     The last factor (2.d-2) converts back to a diameter in
c     micrometers.
      wet_diameter =
     $     (0.4989352162271429d0 * dry_radius**0.3026183900844475d1
     $     / (correction * dry_radius**(-0.1371059101078550d1)
     $     - dlog(relative_humidity)) + dry_radius * dry_radius
     $     * dry_radius) ** 0.3333333333333333333d0 * 2.d4

      end


c     Function: compute_collision_integral
c
c     Computes the "collision integral" for gas-phase diffusion, based
c     on Lennard-Jones potential. Computation is based on tabulations
c     from: Hirschfelder, J.O., Curtiss, C.O., and Bird, R.B, Molecular
c     Theory of Gases and Liquids, Wiley, New York, 1954. See also
c     Poling, B.E., Prausnitz, J.M., and O'Connell, J.P, The Properties
c     of Gases and Liquids, McGraw-Hill, 2000.
c
c     Parameter:
c     x - x is equal to "k T / e" where k is the Boltzmann constant, T
c     is the temperature and e is the Lennard-Jones molecular
c     interaction parameter.
c
c     Returns:
c     omega - The collision integral.
      subroutine ssh_compute_collision_integral(x, omega)

      double precision x

      double precision omega

      integer i

c     Tabulated quantities.
      double precision x_tab(81), omega_tab(81)

      data x_tab / 0.30d0, 0.35d0, 0.40d0, 0.45d0, 0.50d0, 0.55d0,
     $     0.60d0, 0.65d0, 0.70d0, 0.75d0, 0.80d0, 0.85d0, 0.90d0,
     $     0.95d0, 1.00d0, 1.05d0, 1.10d0, 1.15d0, 1.20d0, 1.25d0,
     $     1.30d0, 1.35d0, 1.40d0, 1.45d0, 1.50d0, 1.55d0, 1.60d0,
     $     1.65d0, 1.70d0, 1.75d0, 1.80d0, 1.85d0, 1.90d0, 1.95d0,
     $     2.00d0, 2.10d0, 2.20d0, 2.30d0, 2.40d0, 2.50d0, 2.60d0,
     $     2.70d0, 2.80d0, 2.90d0, 3.00d0, 3.10d0, 3.20d0, 3.30d0,
     $     3.40d0, 3.50d0, 3.60d0, 3.70d0, 3.80d0, 3.90d0, 4.00d0,
     $     4.10d0, 4.20d0, 4.30d0, 4.40d0, 4.50d0, 4.60d0, 4.70d0,
     $     4.80d0, 4.90d0, 5.00d0, 6.00d0, 7.00d0, 8.00d0, 9.00d0,
     $     1.0d01, 2.0d01, 3.0d01, 4.0d01, 5.0d01, 6.0d01, 7.0d01,
     $     8.0d01, 9.0d01, 1.0d02, 2.0d02, 4.0d02 /

      data omega_tab / 2.662d0, 2.476d0, 2.318d0, 2.184d0, 2.066d0,
     $     1.966d0, 1.877d0, 1.798d0, 1.729d0, 1.667d0, 1.612d0,
     $     1.562d0, 1.517d0, 1.476d0, 1.439d0, 1.406d0, 1.375d0,
     $     1.346d0, 1.320d0, 1.296d0, 1.273d0, 1.253d0, 1.233d0,
     $     1.215d0, 1.198d0, 1.182d0, 1.167d0, 1.153d0, 1.140d0,
     $     1.128d0, 1.116d0, 1.105d0, 1.094d0, 1.084d0, 1.075d0,
     $     1.057d0, 1.041d0, 1.026d0, 1.012d0, 9.996d-01, 9.878d-01,
     $     9.770d-01, 9.672d-01, 9.576d-01, 9.490d-01, 9.406d-01,
     $     9.328d-01, 9.256d-01, 9.186d-01, 9.120d-01, 9.058d-01,
     $     8.998d-01, 8.942d-01, 8.888d-01, 8.836d-01, 8.788d-01,
     $     8.740d-01, 8.694d-01, 8.652d-01, 8.610d-01, 8.568d-01,
     $     8.530d-01, 8.492d-01, 8.456d-01, 8.422d-01, 8.124d-01,
     $     7.896d-01, 7.712d-01, 7.556d-01, 7.424d-01, 6.640d-01,
     $     6.232d-01, 5.960d-01, 5.756d-01, 5.596d-01, 5.464d-01,
     $     5.352d-01, 5.256d-01, 5.130d-01, 4.644d-01, 4.170d-01 /

      if (x.le.x_tab(1)) then
         omega = omega_tab(1)
      else if (x.ge.x_tab(81)) then
         omega = omega_tab(81)
      else
         i = 2
c     Searches for i so that x_tab(i-1) < x <= x_tab(i).
         do while (x.gt.x_tab(i))
            i = i + 1
         end do
         omega = omega_tab(i-1) + (omega_tab(i) - omega_tab(i-1))
     $        / (x_tab(i) - x_tab(i-1)) * (x - x_tab(i-1))
      endif

      end


c     Function: compute_condensation_transfer_rate
c
c     Computes condensation transfer rate for a given species.
c
c     Parameters:
c     diffusivity - Gas-phase diffusivity (m^2/s).
c     velocity - Quadratic mean velocity (m/s).
c     accomodation - Accomodation coefficient.
c     wet_diameter - Aerosol wet diameter (micrometer).
c
c     Returns:
c     rate - Condensation/evaporation transfer rate (m^3/s). The
c     condensation growth rate is then given by "rate * (c^g -- c^s)"
c     where c^g is the gaseous concentration and c^s is the
c     concentration at the aerosol surface.
      subroutine ssh_compute_condensation_transfer_rate(diffusivity,
     $     velocity, accomodation, wet_diam, rate)

      double precision diffusivity
      double precision velocity
      double precision accomodation
      double precision wet_diam

      double precision rate

c     Wet diameter in meter.
      double precision diameter_m
c     Knudsen number.
      double precision knudsen

c     Conversion of wet diameter from micrometer to meter, with a
c     threshold of 1.d-3 micrometer.
      diameter_m = 1.d-6 * dmax1(wet_diam, 1.d-3)

      knudsen = 4.d0 * diffusivity / (velocity * diameter_m)

      if (knudsen.le.1.d-2) then
c     Continuous regime.
         rate = 6.283185307178d0 * diffusivity * diameter_m
      else if (knudsen.ge.5.d0) then
c     Free molecular regime.
         rate = 3.141592653589d0 * 2.5d-1 * velocity * accomodation
     $        * diameter_m * diameter_m
      else
c     Transition regime. See Dahneke, 1983.
         rate = 6.283185307178d0 * diffusivity * diameter_m
     $        * (1.d0 + knudsen)
     $        / (1.d0 + 2.d0 * (1.d0 + knudsen)
     $        * knudsen / accomodation)
      end if

      end


c     Function: compute_quadratic_mean_velocity
c
c     Computes condensation transfer rate for a given species.
c
c     Parameters:
c     temperature - Temperature (K).
c     weight - Molecular weight (g/mol).
c
c     Returns:
c     velocity - Quadratic mean molecular velocity (m/s).
      subroutine ssh_compute_quadratic_mean_velocity(temperature, 
     $     weight, velocity)

      double precision temperature
      double precision weight

      double precision velocity

      velocity = dsqrt(2.11714271498563d4 * temperature / weight)

      end


c     Function: compute_saturation_concentration
c
c     Computes saturation concentration of an organic species. The
c     reference temperature is 298 K.
c
c     Parameters:
c     temperature - Temperature (K).
c     weight - Molecular weight (g/mol).
c     enthalpy - Enthalpy of vaporization (J/mol).
c     saturation_pressure - Saturation vapor pressure (Pa).
c
c     Returns:
c     concentration - Saturation concentration (microgram/m^3).
      subroutine ssh_compute_saturation_concentration(temperature, 
     $     weight,enthalpy, saturation_pressure, concentration)

      double precision temperature
      double precision weight
      double precision enthalpy
      double precision saturation_pressure

      double precision concentration

      double precision ratio

c     Constant "0.120279047389945" is equal to 1 / R where R is the
c     molar (perfect) gas constant (8.314 J/K/mol).
      ratio = 0.120279047389945d0 / temperature

c     Constant "4.03620964395788d-4" is equal to 1 / (298. * R) where
c     298 K is the reference temperature.
      concentration = saturation_pressure * 1.d6 * weight * ratio
     $     * dexp(-enthalpy * (ratio - 4.03620964395788d-4))

      end


c     Function: compute_kelvin_coefficient
c
c     Computes the correction factor due to the Kelvin effect on the
c     aerosol surface concentration.
c
c     Parameters:
c     temperature - Temperature (K).
c     weight - Molecular weight (g/mol).
c     surface_tension - Aerosol surface tension (N/m).
c     wet_diameter - Aerosol wet diameter (micrometer).
c     density - Aerosol mass density (kg/m^3).
c
c     Returns:
c     coefficient - Kelvin effect coefficient.
      subroutine ssh_compute_kelvin_coefficient(temp, weight,
     $     surf_tension, wet_diam, density, coefficient)

      double precision temp
      double precision weight
      double precision surf_tension
      double precision wet_diam
      double precision density

      double precision coefficient

      coefficient = dexp(4.d3 * surf_tension * weight
     $     / (8.314d0 * temp * density * wet_diam))
      

      end
