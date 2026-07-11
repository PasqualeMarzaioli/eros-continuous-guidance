This folder contains the kernel of the possible transfer trajectory for the LUMOS mission, updated at July 8th, 2025. 
The kernel has an approximate duration of 92 days.

The files are organised as follows:

* Transfer_D_6.bsp is a kernel of the WSB transfer trajectory computed in the ephemeris model.
	Spice ID          : -1320464
	Start of Interval : 2032 NOV 14 05:16:23.639 (TDB)
	End of Interval   : 2033 FEB 14 07:48:08.552 (TDB)


All the trajectories have been generated using the following kernels: 
(downloaded from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/)

1) Leap Seconds Kernels (LSK)
      a) Text LSK (.tls)
            - naif0012.tls

2) Planetary Constant Kernels (PCK)     
      a) Text PCK (.tpc)
            - gm_de440.tpc
            - pck00011.tpc
      b) Binary PCK (.bpc)
	    - earth_200101_990825_predict.bpc
            - moon_pa_de440_200625.bpc

3) Kernels for ephemeris of vehicles, planets, satellites, comets, asteroids (SPK)
      a) Binary SPK (.bsp)
            - de440s.bsp

The dynamics include the attraction of the following bodies:
	'MERCURY BARYCENTER'
	'VENUS BARYCENTER'
        'EARTH'
        'MARS BARYCENTER'
        'JUPITER BARYCENTER'
        'SATURN BARYCENTER'
        'URANUS BARYCENTER'
        'NEPTUNE BARYCENTER'
        'PLUTO BARYCENTER
        'MOON'
        'SUN'
In addition, the solar radiation pressure (SRP) is considered using the cannonball model.

Finally, the following constants have been used:
- m = 400 kg				[Spacecraft mass]
- A = 3.65 m^s				[Spacecraft area]
- eps = 0.5                             [Spacecraft reflectivity coefficient]
- Psun = 1371 W/m^2			[Solar flux @ 1AU]
