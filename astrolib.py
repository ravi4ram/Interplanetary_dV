
import numpy as np
import pandas as pd
from skyfield.api import load, Loader
from skyfield.elementslib import osculating_elements_of

import lambert as lb

# format for np array printing
float_formatter = "{:9.2f}".format
np.set_printoptions(formatter={'float_kind':float_formatter})

# ----------------------------------------------------------------
# class to provide functions for kernal loading,
# utc time objects, planet name validation and
# planets state vector, classical orbital elements estimation.
# hohmann interplanetary transfer time
# lambert estimation
# ----------------------------------------------------------------
# Author: ravi_ram
# ----------------------------------------------------------------

class astrolib:
    # skyfield init
    jpl_ephemeris_path = r'./data'
    # open the JPL ephemeris DE421
    jpl_ephemeris      = r'de421.bsp'
    
    # constructor
    def __init__(self): #, jpl_ephemeris_path, jpl_ephemeris):
        self.path   = self.jpl_ephemeris_path
        self.kernel = self.jpl_ephemeris

        # avoids multiple copies of large files
        load = Loader(self.path)
        self.planets = load(self.kernel)
        # return
        return
    
    # get current system time
    def get_current_time(self):
        # Create a timescale
        ts = load.timescale(builtin=True) # load.timescale()
        # and ask the current time.
        t = ts.now()
        # return it
        return t
    
    # get time object
    def get_utc_time(self, d, m, yyyy):
        # Create a timescale
        ts = load.timescale(builtin=True) # load.timescale()
        # and set the time.
        t = ts.utc(yyyy, m, d)
        # return it        
        return t
    
    # verify planet name (remove empty spaces and capitalizes the first character)
    def check_planet_name_old(self, planet_name):
        # planets list  
        planets_list = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter',
                   'Saturn', 'Uranus', 'Neptune', 'Pluto']
        try:
            # remove empty spaces and capitalizes the first character
            planet_name = planet_name.strip().capitalize()
            # find the index
            #ind = planets.index(planet_name)
            #return ind
            if planet_name in planets_list:
                ind = planets_list.index(planet_name)
                return planet_name, ind
            else:
                raise ValueError()
            # end try
        except (ValueError, IndexError):
            print('error: Planet \'',planet_name,'\' not in list.')
            exit(2501)
            
    def get_planet_data(self, planet_name):
        
        # https://ssd.jpl.nasa.gov/astro_par.html
        # https://www.jpl.nasa.gov/edu/pdfs/scaless_reference.pdf
        
        # Mass (in kg) and radius (in km), distance from sun (km), GM (in km**3/s**2)
        body = {'Sun': (1.988e30, 6.955e5, 0.0, 132712440041.279419),
                'Mercury': (3.301e23, 2440.0, 57900000.0, 22031.868551),
                'Venus': (4.867e+24, 6052.0, 108200000.0, 324858.592000),
                'Earth': (5.972e24, 6371.0, 149600000.0, 398600.435507),
                'Mars': (6.417e23, 3390.0, 227900000.0, 42828.375816),
                'Jupiter': (1.899e27, 69911.0, 778600000.0, 126712764.100000),
                'Saturn': (5.685e26, 58232.0, 1433500000.0, 37940584.841800),
                'Uranus': (8.682e25, 25362.0, 2872500000.0, 5794556.400000),
                'Neptune': (1.024e26, 24622.0, 4495100000.0, 6836527.100580),
                'Pluto': (0.01303e24, 1188.0, 5906380000.0, 975.500000)
               }

        try:
            # remove empty spaces and capitalizes the first character
            planet_name = planet_name.strip().capitalize()

            if planet_name in body: 
                mass, radius, distance, GM = body[planet_name]
                return planet_name, mass, radius, distance, GM
            else:
                raise ValueError()
            # end try
        except (ValueError, IndexError):
            print('error: Planet \'',planet_name,'\' not in list.')
            exit(2501)
        # dummy    
        return            

    # extract planet's state vector and classical orbital elements
    # planet_name  - valid name from Mercury through Pluto
    # time         - single/list of times 
    # returns
    #     mu, state vector, orbital elements
    def get_planet_ephemeris(self, planet_name, time):
        planet = planet_name + " Barycenter"
        position = self.planets[planet].at(time)
        r_vec = position.position.km
        v_vec = position.velocity.km_per_s
        r    = np.linalg.norm(r_vec)
        v    = np.linalg.norm(v_vec)        
        state_vector = np.array([r_vec, v_vec])
        
        # J2000.0 equatorial plane
        elements = osculating_elements_of(position)
        # GM_sun
        mu   = elements._mu
        # Semimajor axis
        a    = elements.semi_major_axis.km
        # orbital inclination
        i    = elements.inclination.radians
        # eccentricity
        e    = elements.eccentricity
        # angular momentum :
        h_vec = elements._h_vec
        # Specific angular momentum
        h    = np.sqrt( (h_vec * h_vec).sum(axis=0) )

        # Ω (Omega) = Right Ascension of the Ascending Node
        raan = elements.longitude_of_ascending_node.radians  
        # ω (omega) = Argument of Periapsis
        aop  = elements.argument_of_periapsis.radians
        # ν (nu) = True Anomaly
        nu   = elements.true_anomaly.radians
        # period_in_days
        T    = elements.period_in_days
        
        # classical orbital elements
        coe  = np.array( [h, e, i, raan, aop, nu, a, T] )
        
        # return      
        return mu, state_vector, coe
    
    # Calculate classical orbital elements from the state vector
    def coe_from_sv(self, R, V, mu):
        deg = np.pi/180 # multiple this factor radians    
        #eps = 0;
        eps = 1.e-6;
        
        # calculate the distance - norm(R)
        r    = np.linalg.norm(R) 
        # calculate the speed - norm(V)
        v    = np.linalg.norm(V)
        # calculate the radial velocity
        # if vr > 0, the spacecraft is flying away from perigee.
        # if vr < 0, it is flying toward perigee.
        vr   = np.dot(R,V)/r 
        # specific angular momentum h = R X V
        H    = np.cross(R,V) 
        # magnitude of the specific angular momentum, h = sqrt(h.h)
        h    = np.linalg.norm(H) # first orbital element

        # calculate the inclination (equation 4.7) - i = cos-inverse(hz/h)
        # If 90° < i <= 180°, the angular momentum h points in a southerly direction.
        # In that case, the orbit is retrograde
        incl = np.arccos(H[2]/h)     # second orbital element


        # calculate the node line (equation 4.8) N = K-hat X h
        N    = np.cross(np.array([0, 0, 1]).T,H)
        # alculate the magnitude of N
        n    = np.linalg.norm(N)

        # calculate the right ascension of the ascending node, Ω = cos-inv (Nx/N)
        # equation 4.9 (incorporating the case incl = 0):
        if incl != 0:       # inclined orbit
            RA = np.arccos(N[0]/n)
            # if Ny < 0 then 180° <= Ω < 360°. RA = 360° - RA
            if N[1] < 0:
                RA = 2*np.pi - RA
        else:                # equatorial orbit
            RA = 0    

        # calculate the eccentricity vector (equation 4.10)
        E = 1/mu*((v**2 - mu/r)*R - r*vr*V)
        e = np.linalg.norm(E)            # fourth orbital element

        # calculate the argument of perigee (equation 4.12)
        # fifth orbital element
        if incl != 0:                    # inclined orbit
            if e > eps:                  # non-circular orbit
                w = np.arccos(np.dot(N/n,E/e))
                if E[2] < 0:             # eZ < 0 -> 180° < ω < 360°
                    w = 2*np.pi - w
            else:                        # circular orbit
                w = 0 
        else:                            # equatorial orbit
            if e > eps:                  # non-circular orbit
                w = np.arccos(E[0]/e)
                if E[1] < 0:
                    w = 2*pi - w
            else:                        # circular orbit 
                w = 0
        
        # calculate the true anomaly (equation 4.13a)
        # (incorporating the cases incl = 0 and e = 0)
        if incl != 0:                    # inclined orbit
            if e > eps:                  # non-circular orbit
                TA = np.arccos(np.dot(E/e,R/r))
                if vr < 0:
                    TA = 2*np.pi - TA
            else:                        # circular orbit
                TA = np.arccos(dot(N/n,R/r))
                if R[2] < 0:
                    TA = 2*np.pi - TA
        else:                            # equatorial orbit
            if e > eps:                  # non-circular orbit
                TA = np.arccos(dot(E/e,R/r) )
                if vr < 0:
                    TA = 2*np.pi - TA 
            else:                        # circular orbit
                TA = np.arccos(R[0]/r)
                if R[1] < 0:
                    TA = 2*np.pi - TA

        # perigee and apogee radii
        rp = (h**2/mu) * (1/(1 + e * np.cos(np.radians(0))) )
        ra = (h**2/mu) * (1/(1 + e * np.cos(np.radians(180))) )
        
        # semimajor axis of the ellipse
        a = 1/2 * (rp + ra)
        
        # or other form of the eqn
        
        # semimajor axis of the ellipse, equation 4.62 (a < 0 for a hyperbola)   
        #a = h**2/mu/(1 - e**2)
        
        # period
        T = (2 * np.pi * np.power(a, 1.5)) / np.sqrt(mu)
        T = T / 86400.0 # in days
        
        # pack the data - orbital elemnts (h, e, i, raan, aop, true-anomaly, a, T)
        coe = [h, e, incl, RA, w, TA, a, T] #[h, e, RA, incl, w, TA, a]
        
        return coe    

    # print the planets state vector (r,v) 
    # and orbital elemnts (h, e, i, raan, aop, nu, a, T)
    def print_coe(self, planet, state_vector, coe, dt_utc):
        r_vec = state_vector[0]
        v_vec = state_vector[1]
        [h, e, i, raan, aop, nu, a, T] = coe
        
        print('-' * 50)
        print('Planet              :', planet, ' @Date [', dt_utc.utc_strftime('%d-%m-%Y'), ']')
        print('-' * 50)
        print('Position          p :', r_vec, 'km')
        print('Velocity          v :', v_vec, 'km/s')        
        print('Inclination       i : {0:.2f} degrees'.format( np.degrees(i)) )
        print('Eccentricity      e : {0:.5f}'.format(e))
        print('Angular Momentum  h : {0:.2f}'.format(h))
        print('Semimajor axis    a : {0:.0f} km'.format(a))
        print('RAAN              Ω : {0:.0f} degrees'.format( np.degrees(raan)) )
        print('AoP               ω : {0:.2f} degrees'.format( np.degrees(aop)) )
        print('True Anomaly      ν : {0:.2f} degrees'.format( np.degrees(nu)) )
        print('Period            T : {0:.2f} days'.format(T) )
        print()
        # end print_coe
        return

    # calculate the time of interplanetary transfer
    # r1 : semimajor axis of the departing planet
    # r2 : semimajor axis of the arrival planet
    def get_hohmann_tof(self, r1, r2, mu):
        # semi major axis
        a_t = (r1 + r2) / 2.0
        # travel time in seconds
        #t_12_secs = np.pi / np.sqrt(mu) * a_t**(3/2)
        t_12_secs = np.pi * np.sqrt( a_t**3 / mu )
        t_12_days = t_12_secs / 86400
        # return tof
        return t_12_secs, t_12_days #t_12_secs

    # departure phase angle (transfer from 1 to 2)
    # planet period in secs
    # time of flight in secs
    def get_departure_phase_angle(self, T1_sec, T2_sec, TOF_sec, mu):
        # mean motions of the planet, in rad/s
        n1 = 2*np.pi/T1_sec
        n2 = 2*np.pi/T2_sec 
        
        # phase angle 
        gamma_1 = (np.pi - n2 * TOF_sec) % (2 * np.pi)
        gamma_2 = (np.pi - n1 * TOF_sec) % (2 * np.pi)
        
        # synodic period (time between successive launch windows is known as a synodic period)
        Tsyn = 2*np.pi / abs(n2 - n1)
        
        # return dep, arr phase angles
        return gamma_1, gamma_2, Tsyn

    # calculate c3 and delv values from lambert solution
    def get_lambert_estimates(self, mu, v_dep, v_arr, r1, r2, flight_time_secs,
                               orb_type, M, path):
        # multiple-solution for m>0 cases. not considered here
        # def solve(mu, r1, r2, t_sec, orb_type, path, m)
        v1_list, v2_list = lb.solve(mu, r1, r2, flight_time_secs, orb_type, path, M)
        v1, v2 = v1_list[0], v2_list[0]
        # Solve the problem
        #v1, v2 = lb.solve(mu, r1, r2, flight_time_secs, M=M, prograde=True, low_path=low_path)   
        # compute v_inf for departure and arrival (subtract planet velocities)
        v_inf_dep = np.linalg.norm(v_dep - v1) 
        v_inf_arr = np.linalg.norm(v_arr - v2)
        # characteristic energy. v_inf = orbital velocity when the
        # orbital distance tends to infinity.
        c3_dep = v_inf_dep**2
        c3_arr = v_inf_arr**2
        # ∆V = Vplanet1(t1) − VT(t1) + Vplanet2(t2) − VT (t2)
        #delv_total = v_inf_dep + v_inf_arr    
        #return 
        return [c3_dep, c3_arr, v_inf_dep, v_inf_arr ]
    # end class