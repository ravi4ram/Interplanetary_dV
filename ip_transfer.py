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
import pandas as pd
from datetime import datetime

from plot import plot
from astrolib import astrolib

# format for np array printing
float_formatter = "{:9.2f}".format
np.set_printoptions(formatter={'float_kind':float_formatter})


# function to find a optimal arrival date. Selecting a set of
# days around hohmann tof day, and applying minimum of 
# weighted mean average on c3, delta-v and phase angle data points
# to a optimal arrival date.
def find_optimal_date(from_planet, to_planet, start_date_string):
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
    planet_1, mass_1, radius_1, distance_1, GM_1 = alib.get_planet_data(from_planet)
    mu, sv_dep, coe_dep = alib.get_planet_ephemeris(planet_1, t1)
    r_dep, v_dep = sv_dep[0], sv_dep[1]
    a1 = coe_dep[6]
    T1_sec = coe_dep[7]*86400
    
    # print sv and orbital elements
    alib.print_coe(planet_1, sv_dep, coe_dep, t1)

    # arrival at - planet 2
    planet_2, mass_2, radius_2, distance_2, GM_2 = alib.get_planet_data(to_planet)
    mu, sv, coe = alib.get_planet_ephemeris(planet_2, t1)
    r_vec, v_vec = sv[0], sv[1]
    a2 = coe[6]
    T2_sec = coe[7]*86400

    # print sv and orbital elements     
    alib.print_coe(planet_2, sv_dep, coe_dep, t1)
    
    # hohmann tof
    tof_secs, tof_days = alib.get_hohmann_tof(a1, a2, mu)
    hohmann_tof = tof_days

    # check for multiple days from minimum hohmann tof
    data_out   = []

    tof_days   = tof_days - int(tof_days/4.0) #40
    tof_days_i = int(tof_days / 4.0)
    tof_days_f = int(tof_days * (3.0 / 4.0) )
    dt         = int(tof_days_f - tof_days_i)
    tof_days   = int(tof_days - tof_days_i)
    
    # scan through the range of dates
    for dt in range(0, dt, 1):
        tof_days = tof_days + 1
        tof_secs = tof_days * 86400
        t2       = t1 + tof_days
        
        # final position of planet 1
        mu, sv, coe = alib.get_planet_ephemeris(planet_1, t2)
        r_vec, v_vec = sv[0], sv[1]
        # print sv and orbital elements     
        #alib.print_coe(planet_1, sv, coe, t2)
        
        # final position of planet 2
        mu, sv_arr, coe_arr = alib.get_planet_ephemeris(planet_2, t2)
        r_arr, v_arr = sv_arr[0], sv_arr[1]
        # print sv and orbital elements     
        #alib.print_coe(planet_2, sv_arr, coe_arr, t2)
        #print('.' * 50) 
        
        # lambert estimation (type-I)
        orb_type, M, low_path = 'prograde', 0.0, 'low'
        c3_and_delv = alib.get_lambert_estimates(mu, v_dep, v_arr, r_dep, r_arr,
                                               tof_secs, orb_type, M, low_path)
        c3_dep, c3_arr, v_inf_dep, v_inf_arr = c3_and_delv
        # ∆V = Vplanet1(t1) − VT(t1) + Vplanet2(t2) − VT (t2)
        delv_total = v_inf_dep + v_inf_arr
        
        # dep arr phase angles
        gamma1, gamma2, Tsyn = alib.get_departure_phase_angle(T1_sec, T2_sec, tof_secs, mu)

        # pack result            
        res = [mu,
               planet_1, radius_1, distance_1, GM_1,
               planet_2, radius_2, distance_2, GM_2,
               v_inf_dep, v_inf_arr]
        
        # calculate delta-v
        dVe_1, dVe_2 = calculate_dV(*res)
        
        # table start, end, duration, gamma1, gamma2, dep_dVe, arr_dVe, total_dv
        out = np.array([t1.utc_strftime('%d-%m-%Y'), t2.utc_strftime('%d-%m-%Y'),
                        round(tof_days, 3),
                        #round(v_inf_dep, 3), round(v_inf_arr, 3),
                        round(c3_dep, 3), round(c3_arr, 3),
                        round(np.degrees(gamma1), 3), round(np.degrees(gamma2), 3) ,
                        round(dVe_1, 3), round(dVe_2, 3), round( (dVe_1+dVe_2), 3) ])
        # data for different tof_days 
        data_out.append(out)
    # end for loop
    
    # create dataframe 
    df = pd.DataFrame(data_out)
    df.columns =['Departure Date', 'Arrival Date', 'Days',
                 'C3 Departure', 'C3 Arrival',
                 'Departure Phase', 'Arrival Phase',
                 'Departure delV', 'Arrival delV', 'Total delV']

    # set column datatypes
    df['Departure Date'] = pd.to_datetime(df['Departure Date'],
                                          dayfirst=True, exact=True, format='%d-%m-%Y')
    df['Arrival Date']   = pd.to_datetime(df['Arrival Date'],
                                          dayfirst=True, exact=True, format='%d-%m-%Y')
    df['C3 Departure']   = df['C3 Departure'].astype(float)
    df['C3 Arrival']     = df['C3 Arrival'].astype(float)
    df['Departure Phase']= df['Departure Phase'].astype(float)
    df['Arrival Phase']  = df['Arrival Phase'].astype(float)
    df['Departure delV'] = df['Departure delV'].astype(float)
    df['Arrival delV']   = df['Arrival delV'].astype(float)
    df['Total delV']     = df['Total delV'].astype(float)

    # find optimal date

    # columns used for weighted mean average 
    cols = ['C3 Departure', 'C3 Arrival',
            'Departure Phase', 'Arrival Phase',
            'Departure delV', 'Arrival delV']
    
    # copy relevant fields
    df1 = df[cols].copy(deep=True) 
    # inverse for finding min
    df1['Departure Phase'] = 360 - df1['Departure Phase']
    # diff
    df1 = df1 - df1.min()
    # cols weight based on importance
    wt = [10, 10, 5, 13, 50, 50]
    # df1['weighted_mean'] using numpy average function
    df1['weighted_mean'] = np.average(df1[cols], weights=wt, axis=1)
    # get index of the the lowest weighted mean
    index = df1[['weighted_mean']].idxmin()
    
    # get optimal values pointed by the index
    optimal_row = df.iloc[index]
    optimal_date = df.iloc[index]['Arrival Date'].astype('str')

    # convert date string in '%d-%m-%Y' format
    opt_arr_date = datetime.strptime(optimal_date.to_list()[0], '%Y-%m-%d').strftime('%d-%m-%Y')

    # tof days
    tof_days_optimal  = df.iloc[index]['Days'].astype('float').to_list()[0]
    t2 = t1 + tof_days_optimal

    print()
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
    
    # print optimal del-v
    planet_1 = planet_1.ljust(8)
    planet_2 = planet_2.ljust(8)
   
    dVe_1 = df.iloc[index]['Departure delV'].astype('float').to_list()[0]
    dVe_2 = df.iloc[index]['Arrival delV'].astype('float').to_list()[0]
    dVe_total = df.iloc[index]['Total delV'].astype('float').to_list()[0]
    
    print('.' * 50)    
    print('Hohmann TOF         : {0:.2f} days'.format(hohmann_tof) )
    print('Optimal TOF         : {0:.2f} days'.format(tof_days_optimal) )    
    print('Optimal arrival date:', opt_arr_date )
    print('dVe at',planet_1,'TOI :', round(dVe_1, 3)  )
    print('dVe at',planet_2,'TOI :', round(dVe_2, 3)  )
    print('Total dVe           :', round( (dVe_total), 3)  )
    print('.' * 50)

    # return dataframe, optimal_row, optimal arrival date
    return df, optimal_row, opt_arr_date    

# function to calculate delta-v from an assumed orbit of
# [170 X 36000] km at TOI of departure planet and an arrival of
# [500 X 60000] km at TOI of arrival planet.
def calculate_dV(mu,
                 planet_1, radius_1, distance_1, GM_1,
                 planet_2, radius_2, distance_2, GM_2,
                 v_inf_dep, v_inf_arr):

    rp_1 = radius_1 + 170.00    # km
    ra_1 = radius_1 + 36000.00  # km
  
    # perigee of the hyperbola
    vp_1 = np.sqrt((2*GM_1/rp_1) + v_inf_dep**2)
    # Δv-value of the boost to be performed at transfer orbit injection (TOI) maneuver
    # circular orbital velocity = np.sqrt(GM_1/rp_1)
    dVe_1 = vp_1 - np.sqrt( GM_1 * (2*ra_1/(rp_1*(ra_1 + rp_1 ))) ) 
    #print('dVe at',planet_1,'TOI :', round(dVe_1, 3)  )

    # turning angle, angle between the transverse axis of the departure hyperbola
    # and the direction of the velocity vector of the Earth at the moment of departure is
    alpha_1 = np.arccos(1/((rp_1 * v_inf_dep**2 / GM_1) + 1))
    # In case of a prograde orbit, the perigee is placed on the dark side.

    rp_2 = radius_2 + 500.00    # km
    ra_2 = radius_2 + 60000.00   # km
    # perigee of the hyperbola
    vp_2 = np.sqrt((2*GM_2/rp_2) + v_inf_arr**2)
    # Δv-value of the boost to be performed at TOI maneuver
    # circular orbit orbital velocity = np.sqrt(GM_2/rp_2)
    dVe_2 = vp_2 - np.sqrt( GM_2 * (2*ra_2/(rp_2*(ra_2 + rp_2 ))) ) 
    #print('dVe at',planet_2,'TOI :', round(dVe_2, 3)  )
    #print('Total dVe           :', round( (dVe_1+dVe_2), 3)  )

    # The angle between the transverse axis and the asymptote of the arrival hyperbola
    alpha_2 = np.arccos(1/((rp_2 * v_inf_arr**2 / GM_2) + 1))
    
    # return
    return dVe_1, dVe_2

# function to create data for pork-chop plot using the calculated optimal
# arrival date
def get_pork_data(from_planet, start_date, to_planet, arrival_date):
    
    # 'delv_plot' or 'c3_plot'
    plot_type = 'delv_plot'
    
    # create object
    alib = astrolib()  
    
    # starting from - planet 1    
    planet_1, mass_1, radius_1, distance_1, GM_1 = alib.get_planet_data(from_planet)
    # arrival at - planet 2
    planet_2, mass_2, radius_2, distance_2, GM_2 = alib.get_planet_data(to_planet)    
    
    # create datetime object from string
    dt_dep, dt_arr = 0, 0
    try:
        dt_dep = datetime.strptime(start_date, '%d-%m-%Y')
        dt_arr = datetime.strptime(start_date, '%d-%m-%Y')
    except ValueError:
        print ('error : wrong date format. [verify as dd-mm-yyyy]')

    dep_window, arr_window = 40, 30
    # dep_start
    # start_date  = '29-03-2028'    #(format - dd-mm-yyyy)
    sd = pd.to_datetime(start_date, dayfirst=True, exact=True, format='%d-%m-%Y')
    smin = sd - pd.Timedelta(days=dep_window) 
    smax = sd + pd.Timedelta(days=dep_window) 
    dep_dates = pd.date_range(smin, smax, freq='D').strftime('%Y-%m-%d')
    
    # arr_date  = '19-07-2028'    #(format - dd-mm-yyyy)
    ad = pd.to_datetime(arrival_date, dayfirst=True, exact=True, format='%d-%m-%Y')
    amin = ad - pd.Timedelta(days=arr_window) 
    amax = ad + pd.Timedelta(days=arr_window) 
    arr_dates = pd.date_range(amin, amax, freq='D').strftime('%Y-%m-%d')
    
    # dataframe matrix with row-column names as dates
    df = pd.DataFrame(np.ones((len(arr_dates), len(dep_dates)), int), arr_dates, dep_dates, dtype='object')

    # generate result matrix
    rows, cols = df.shape[0], df.shape[1]
    
    # create matrix of data between arrival(ix,rows) and departure (iy, cols)
    for ix in range(0, rows):
        # arrival
        row_name = arr_dates.to_list()[ix]
        arr_date = datetime.strptime(row_name, '%Y-%m-%d')
        # t2 UTC
        t2 = alib.get_utc_time(arr_date.day, arr_date.month, arr_date.year )
        
        for iy in range(0, cols):
            # departure
            col_name = dep_dates.to_list()[iy]            
            dep_date = datetime.strptime(col_name, '%Y-%m-%d')
            # t1 UTC
            t1 = alib.get_utc_time(dep_date.day, dep_date.month, dep_date.year )

            # difference between times, in days
            tof_days = t2 - t1
            tof_secs = tof_days*86400

            # departure conditions
            mu, sv_dep, coe_dep = alib.get_planet_ephemeris(planet_1, t1)
            r_dep, v_dep = sv_dep[0], sv_dep[1]
            a1 = coe_dep[6]
            T1_sec = coe_dep[7]*86400
            
            # arrival conditions
            mu, sv_arr, coe_arr = alib.get_planet_ephemeris(planet_2, t2)
            r_arr, v_arr = sv_arr[0], sv_arr[1]
            a2 = coe_arr[6]
            T2_sec = coe_arr[7]*86400
            
            # lambert estimation (type-I)
            orb_type, M, low_path = 'prograde', 0.0, 'low'
            c3_and_delv = alib.get_lambert_estimates(mu, v_dep, v_arr, r_dep, r_arr,
                                                   tof_secs, orb_type, M, low_path)
            c3_dep, c3_arr, v_inf_dep, v_inf_arr = c3_and_delv

            # ∆V = Vplanet1(t1) − VT(t1) + Vplanet2(t2) − VT (t2)
            v_inf_total = v_inf_dep + v_inf_arr
            
            # dep-arr phase angles
            gamma1, gamma2, Tsyn = alib.get_departure_phase_angle(T1_sec, T2_sec, tof_secs, mu)

            # pack result            
            res = [mu,
                   planet_1, radius_1, distance_1, GM_1,
                   planet_2, radius_2, distance_2, GM_2,
                   v_inf_dep, v_inf_arr]
            
            # delta-v estimations
            dVe_1, dVe_2 = calculate_dV(*res)
            # total delta-v
            dV_total = dVe_1 + dVe_2
            
            # matrix data [ix][iy]
            data = [dep_date, arr_date, tof_days,
                    c3_dep, c3_arr,
                    v_inf_dep, v_inf_arr, v_inf_total,
                    gamma1, gamma2,
                    dVe_1, dVe_2, dV_total]

            # set list into cell as object
            df.at[row_name, col_name] = data
        # end for-iy
    # end for-ix
    
    #print(df)
    # arrival-ix-rows,  departure-iy-cols
    rows, cols = df.shape[0], df.shape[1]

    # initialize contour lists
    dep_date_str_list = np.empty(cols, dtype=object)
    arr_date_str_list = np.empty(rows, dtype=object)
    # 2d array shape   
    contour_shape = (rows, cols)
    tof_days_list = np.zeros( contour_shape, dtype=np.float64)
    c3_dep_1_list = np.zeros( contour_shape, dtype=np.float64)
    c3_dep_2_list = np.zeros( contour_shape, dtype=np.float64)
    delv_t_1_list = np.zeros( contour_shape, dtype=np.float64)
    delv_t_2_list = np.zeros( contour_shape, dtype=np.float64)
    dV_1_list = np.zeros( contour_shape, dtype=np.float64)
    dV_2_list = np.zeros( contour_shape, dtype=np.float64)    
    # generate result matrix
    for ix in range(0, rows):      # arrival
        #row_name = arr_dates.to_list()[ix]
        for iy in range(0, cols):  # departure
            col_name = dep_dates.to_list()[iy]
            data = df.iloc[ix][col_name]

            # pack the list
            dep_date_str_list[iy] = data[0]
            arr_date_str_list[ix] = data[1]            
            tof_days_list[ix][iy] = data[2]
            c3_dep_1_list[ix][iy] = data[3]
            c3_dep_2_list[ix][iy] = data[4]
            delv_t_1_list[ix][iy] = data[5]
            delv_t_2_list[ix][iy] = data[6]
            dV_1_list[ix][iy]     = data[10]
            dV_2_list[ix][iy]     = data[11]
        # end-for-iy
    # end-for-ix 
    
    # contour levels for c3 and delv   
    c3_levels = [1, 1.2, 1.5, 1.7, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 19, 20, 30, 50]
    # t_levels for transfer dates 
    t_levels  = [50, 100, 110, 120, 150, 170, 200, 220, 250, 280, 300, 350, 400, 450, 500, 550]
    
    # plot    
    if plot_type == 'delv_plot':
        title = 'Porkchop plot (∆V Total = ||$ΔV_{' +planet_1+'}|| + ||ΔV_{'+ planet_2 +'}$||)' + '\n'
        out = [title, start_date, arrival_date, dep_date_str_list, arr_date_str_list,
                         dV_1_list, dV_2_list, c3_levels,
                         tof_days_list, t_levels]
    elif plot_type == 'c3_plot':       
        title = 'Porkchop plot (C3-characteristic energy = $v_{id}^{2}$)' + '\n'
        out = [title, start_date, arrival_date, dep_date_str_list, arr_date_str_list,
                         c3_dep_1_list, c3_dep_2_list, c3_levels,
                         tof_days_list, t_levels]
    else: raise Exception("error: plot types should be c3_plot or delv_plot")   
    
    # out = [title, start_date, arrival_date, dep_date_str_list, arr_date_str_list,
    #         c3_dep_1_list, c3_dep_2_list, c3_levels,
    #         tof_days_list, t_levels]
    return out

# main function
if __name__ == "__main__":
    # Inputs
    from_planet = 'earth'
    to_planet   = 'venus'
    start_date  = '29-03-2028'    #(format - dd-mm-yyyy)
    
    # find optimal arrival date
    df, optimal_row, optimal_arrival_date = find_optimal_date(from_planet, to_planet, start_date)
    # create data for porkchop plot
    out = get_pork_data(from_planet, start_date, to_planet, optimal_arrival_date)
    # combined plot
    plot(df, *out)
