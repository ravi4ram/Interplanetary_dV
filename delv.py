# ----------------------------------------------------------------
# Interplanetary transfer dV estimation
# ----------------------------------------------------------------
# SPICE kernel used for ephemeris estimation 'de421.bsp' will be
#   downloaded by Skyfield’s load() routine for the first time on
#   the current directory.
#
# estimate planet1 [r1, v1, a1] and planet2 [r2, v2, a2] at start date
# departure from planet1 [r_dep, v_dep] = [r1, v1]
# estimate time of flight for hohmann transfer using a1 and a2
# estimate end date ( start date + time of flight)
# estimate planet1 [r1, v1, a1] and planet2 [r2, v2, a2] at end date
# arrival to planet2 [r_arr, v_arr] = [r2, v2]
#
# estimate orbit [v1, v2] using lamberts solver for the given
#   position vectors and time of flight [r_dep, r_arr, tof].
#
# estimate v_inf (asymptotic velocity at infinite distance)
#   for departure and arrival (subtract planet velocities)
#   v_inf_dep = |v_dep - v1| and  v_inf_arr = |v_arr - v2|
#   ∆V_total  = |Vplanet1(t1) − VT(t1)| + |Vplanet2(t2) − VT (t2)|
#   ∆V_total  = v_inf_dep + v_inf_arr
# estimate characteristic energy C3 (measure of the excess
#   specific energy over that required to just barely escape from a massive body)
#   c3_dep = v_inf_dep**2
#   c3_arr = v_inf_arr**2
# ----------------------------------------------------------------
# Author: ravi_ram
# ----------------------------------------------------------------

import numpy as np
from datetime import datetime

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
class astrolib:
    # skyfield init
    jpl_ephemeris_path = r'./'
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
    def check_planet_name(self, planet_name):
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
        t_12_secs = np.pi * np.sqrt( a_t**3/mu )
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


# test all functions
def estimate_delv(from_planet, to_planet, start_date_string):
        
    # create object
    alib = astrolib()
    
    # create datetime object from string
    dt_in = 0
    try:
        dt_in = datetime.strptime(start_date_string, '%d-%m-%Y')
    except ValueError:
        print ('error : wrong date format. [verify as dd-mm-yyyy]')

    # get start time object
    t1 = alib.get_utc_time(dt_in.day, dt_in.month, dt_in.year )

    # starting from - planet 1    
    planet_1, id = alib.check_planet_name(from_planet)
    mu, sv_dep, coe_dep = alib.get_planet_ephemeris(planet_1, t1)
    r_dep, v_dep = sv_dep[0], sv_dep[1]
    a1 = coe_dep[6]
    T1_sec = coe_dep[7]*86400
    
    # print sv and orbital elements    
    alib.print_coe(planet_1, sv_dep, coe_dep, t1)

    # arrival at - planet 2
    planet_2, id = alib.check_planet_name(to_planet)
    mu, sv, coe = alib.get_planet_ephemeris(planet_2, t1)
    r_vec, v_vec = sv[0], sv[1]
    a2 = coe[6]
    T2_sec = coe[7]*86400

    # print sv and orbital elements     
    alib.print_coe(planet_2, sv_dep, coe_dep, t1)

    # time of flight
    tof_secs, tof_days = alib.get_hohmann_tof(a1, a2, mu)
    #print('time of flight  tof : {0:.2f} days'.format(tof_days) )
    
    # final planets position at 
    t2 = t1 + tof_days

    # final position of planet 1
    mu, sv, coe = alib.get_planet_ephemeris(planet_1, t2)
    r_vec, v_vec = sv[0], sv[1]
    # print sv and orbital elements     
    alib.print_coe(planet_1, sv, coe, t2)
    
    # final position of planet 2
    mu, sv_arr, coe_arr = alib.get_planet_ephemeris(planet_2, t2)
    r_arr, v_arr = sv_arr[0], sv_arr[1]
    # print sv and orbital elements     
    alib.print_coe(planet_2, sv_arr, coe_arr, t2)
    print('.' * 50) 
    
    # lambert estimation (type-I)
    orb_type, M, low_path = 'prograde', 0.0, 'low'
    c3_and_delv = alib.get_lambert_estimates(mu, v_dep, v_arr, r_dep, r_arr,
                                           tof_secs, orb_type, M, low_path)
    c3_dep, c3_arr, v_inf_dep, v_inf_arr = c3_and_delv
    # ∆V = Vplanet1(t1) − VT(t1) + Vplanet2(t2) − VT (t2)
    delv_total = v_inf_dep + v_inf_arr
    
    print()    
    print(t1.utc_strftime('Departure : %d-%m-%Y'), ' from ', planet_1)
    print('Position          p :', r_dep, 'km')
    print('Velocity          v :', v_dep, 'km/s')
    print(t2.utc_strftime('Arrival   : %d-%m-%Y'), ' at ', planet_2)
    print('Position          p :', r_arr, 'km')
    print('Velocity          v :', v_arr, 'km/s')
    print()
    print('time of flight  tof : {0:.2f} days'.format(tof_days) )
    # dep arr phase angles
    gamma1, gamma2, Tsyn = alib.get_departure_phase_angle(T1_sec, T2_sec, tof_secs, mu)
    print('dep phase angle  γ1 : {0:.2f} degrees'.format( np.degrees(gamma1)) )
    print('arr phase angle  γ2 : {0:.2f} degrees'.format( np.degrees(gamma2)) )
    print('Synodic period Tsyn : {0:.2f} days'.format(Tsyn/86400) )        
    print('.' * 50) 
    print('v_inf_dep, v_inf_arr:', round(v_inf_dep, 3), round(v_inf_arr, 3) )
    print('delv_total          :', round(delv_total, 3)  )
    print('c3_dep, c3_arr      :', round(c3_dep, 3), round(c3_arr, 3))
    print()
    # return
    return v_inf_dep, v_inf_arr, c3_dep, c3_arr


# main function
if __name__ == "__main__":
    
    from_planet = 'earth'
    to_planet   = 'venus'
    start_date  = '30-03-2028'    #(format - dd-mm-yyyy)
    estimate_delv(from_planet, to_planet, start_date)
    
