# Interplanetary transfer dV estimation
> Program to estimate interplanetary dV for Type-I trajectories using lambert's solver   

## Table of contents
* [General info](#general-info)
* [Results](#Results)
* [Setup](#setup)
* [How to run ](#how)
* [Updates](#updates)
* [To-do list](#to-do)

## NOTE :   
SPICE kernel 'de421.bsp' [~18 MB] used for ephemeris estimation will be downloaded by Skyfield’s load() routine **once** on the current directory.  It will not connect to internet again on future runs.   

## Results

###Input :    
```
from_planet = 'earth'
to_planet   = 'venus'
start_date  = '30-03-2028'
```   
* Arrival date is taken as calculated through the hohmann transfer time-of-flight.  For a different arrival date, edit the line which says  
```t2 = t1 + tof_days```  
accordingly.   
* As satellite orbital velocity is a variable with altitude, latitude etc, I did not include that part. That will be a part of patched conic estimation.   
*  

###Output:   
```
--------------------------------------------------
Planet              : Earth  @Date [ 30-03-2028 ]
--------------------------------------------------
Position          p : [-147231343.69 -22837213.32 -9891752.40] km
Velocity          v : [     4.40    -27.06    -11.73] km/s
Inclination       i : 23.44 degrees
Eccentricity      e : 0.01952
Angular Momentum  h : 4451234226.94
Semimajor axis    a : 149153145 km
RAAN              Ω : 0 degrees
AoP               ω : 95.19 degrees
True Anomaly      ν : 94.40 degrees
Period            T : 363.39 days

--------------------------------------------------
Planet              : Venus  @Date [ 30-03-2028 ]
--------------------------------------------------
Position          p : [-147231343.69 -22837213.32 -9891752.40] km
Velocity          v : [     4.40    -27.06    -11.73] km/s
Inclination       i : 23.44 degrees
Eccentricity      e : 0.01952
Angular Momentum  h : 4451234226.94
Semimajor axis    a : 149153145 km
RAAN              Ω : 0 degrees
AoP               ω : 95.19 degrees
True Anomaly      ν : 94.40 degrees
Period            T : 363.39 days

--------------------------------------------------
Planet              : Earth  @Date [ 21-08-2028 ]
--------------------------------------------------
Position          p : [129706037.59 -72055777.40 -31229648.53] km
Velocity          v : [    14.92     23.30     10.10] km/s
Inclination       i : 23.43 degrees
Eccentricity      e : 0.01682
Angular Momentum  h : 4465644573.58
Semimajor axis    a : 150105689 km
RAAN              Ω : 360 degrees
AoP               ω : 96.66 degrees
True Anomaly      ν : 232.15 degrees
Period            T : 366.87 days

--------------------------------------------------
Planet              : Venus  @Date [ 21-08-2028 ]
--------------------------------------------------
Position          p : [101914515.48 36027873.26 9785980.14] km
Velocity          v : [   -12.16     29.58     14.08] km/s
Inclination       i : 24.44 degrees
Eccentricity      e : 0.00983
Angular Momentum  h : 3792448638.19
Semimajor axis    a : 108239910 km
RAAN              Ω : 8 degrees
AoP               ω : 119.35 degrees
True Anomaly      ν : 253.24 degrees
Period            T : 224.65 days

..................................................

Departure : 30-03-2028  from  Earth
Position          p : [-147231343.69 -22837213.32 -9891752.40] km
Velocity          v : [     4.40    -27.06    -11.73] km/s
Arrival   : 21-08-2028  at  Venus
Position          p : [101914515.48 36027873.26 9785980.14] km
Velocity          v : [   -12.16     29.58     14.08] km/s

time of flight  tof : 144.83 days
dep phase angle  γ1 : 304.85 degrees
arr phase angle  γ2 : 36.52 degrees
Synodic period Tsyn : 568.75 days
..................................................
v_inf_dep, v_inf_arr: 8.227 11.295
delv_total          : 19.522
c3_dep, c3_arr      : 67.681 127.573

```


&nbsp;         

## General info
SPICE kernel 'de421.bsp' used for ephemeris estimation will be downloaded by Skyfield’s load() routine for the first time on the current directory.     
&nbsp;    
estimate planet1 [r1, v1, a1] and planet2 [r2, v2, a2] at start date     
departure from planet1 [r_dep, v_dep] = [r1, v1]    
estimate time of flight for hohmann transfer using a1 and a2     
estimate end date ( start date + time of flight)    
estimate planet1 [r1, v1, a1] and planet2 [r2, v2, a2] at end date    
arrival to planet2 [r_arr, v_arr] = [r2, v2]     
&nbsp;     
estimate orbit [v1, v2] using lamberts solver for the given position vectors and time of flight [r_dep, r_arr, tof].     
&nbsp;   
estimate v_inf (asymptotic velocity at infinite distance) for departure and arrival (subtract planet velocities)   
  v_inf_dep = |v_dep - v1| and  v_inf_arr = |v_arr - v2|    
  ∆V_total  = |Vplanet1(t1) − VT(t1)| + |Vplanet2(t2) − VT (t2)|    
  ∆V_total  = v_inf_dep + v_inf_arr     
&nbsp;    
estimate characteristic energy C3 (measure of the excess specific energy over that required to just barely escape from a massive body)    
  c3_dep = v_inf_dep<sup>2</sup>    
  c3_arr = v_inf_arr<sup>2</sup>    


&nbsp;    


## Reference    

1.  On the Nature of Earth-Mars Porkchop Plots  
[ https://trs.jpl.nasa.gov/bitstream/handle/2014/44336/13-0679_A1b.pdf ]  

2.  Interplanetary Mission Design Handbook: Earth-to-Mars Mission Opportunities 2026 to 2045  
[ https://ntrs.nasa.gov/api/citations/20100037210/downloads/20100037210.pdf ]     

## Setup
Script is written with python (Version: 3.10.12) on linux. Additional modules required :   

* numpy  (tested with Version: 1.21.5 )
* skyfield  (tested with Version: 1.4.5 )

## How to run   
* Verify and install required modules 
* run `python delv.py`. 


## Updates   
*   
*   *  

## To-do list
* 

