import numpy as np

# Gooding algorithm for solving Lambert's problem 
# On the Solution of Lambert's Orbital Boundary-Value Problem
#     [ https://apps.dtic.mil/sti/pdfs/ADA200383.pdf ]
# based on the fortran code
#     [ https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit ]
#
# Author: ravi_ram
#

# solve lambert
# inputs
# orb_type:  'prograde' or 'retrograde'
# m:         number of revolutions
# r1:        position vector 1
# r2:        position vector 2
# t_sec:     time of flight in secs
# orb_type:  'prograde' or 'retrograde'
# path:      'low' or 'high'
# m:         number of revolutions
def solve(mu, r1, r2, t_sec, orb_type, path, m):
    # compute angular momentum h, transfer angle ψ0
    h, psi = get_transfer_angle(orb_type, m, r1, r2)
    # multi-revolution
    psi = m * 2.0*np.pi + psi   
    # get list of radial and tangential velocity components
    n, vr1_list, vt1_list, vr2_list, vt2_list  = vlamb(mu, r1, r2, psi, t_sec)
    # loop through and calculate velocities  
    v1_list, v2_list = [], []    
    for i, (vr1, vt1, vr2, vt2) in enumerate( zip(vr1_list, vt1_list, vr2_list, vt2_list)):
        v1, v2 = calculate_velocities(h, r1, r2, vr1, vt1, vr2, vt2)
        v1_list.append(v1)
        v2_list.append(v2)
    # end - for
    
    # return list of velocities
    return v1_list, v2_list

# calculate transfer angle
# ----- as explained on this paper -----
# The Superior Lambert Algorithm
# (https://amostech.com/TechnicalPapers/2011/Poster/DER.pdf )
# -----------------------------------------
# inputs
# orb_type:  'prograde' or 'retrograde'
# m:         number of revolutions
# r1:        position vector 1
# r2:        position vector 2
# outputs
# h:         angular momentum
# psi:       transfer angle
# 
def get_transfer_angle(orb_type, m, r1, r2):
    # calculate magnitudes r1_n and r2_n
    r1_n = np.linalg.norm(r1)
    r2_n = np.linalg.norm(r2)    
    # angular momentum of the transfer orbit
    h = np.cross(r1, r2) / np.linalg.norm(np.cross(r1, r2))
    # α, which is third component of the angular momentum of the transfer orbit
    alpha  = np.dot(np.array([0, 0, 1]), h)
    # transfer angle between the two given position vectors
    psi0 = np.arccos(np.dot(r1, r2) / (r1_n * r2_n))
    
    # For direct transfer orbit: ψ = 2π − ψ0 , if α < 0
    # For retrograde transfer orbit: ψ = 2π − ψ0 , if α > 0
    if alpha < 0:
        if orb_type == 'prograde':
            psi = 2*np.pi - psi0
            h = -h
        else:  psi = psi0
    # end - if
            
    if alpha > 0:
        if orb_type == 'retrograde':
            psi = 2*np.pi - psi0
            h = -h
        else:  psi = psi0
    # end - if
            
    # return h and psi
    return h, psi

# velocities from radial and tangential components
# inputs
# h:      angular momentum
# r1:     position vector 1
# r2:     position vector 2
# vr1:    radial velocity of the initial position vector
# vt1:    tangential velocity of the initial position vector
# vr2:    radial velocity of the final position vector
# vt2:    tangential velocity of the final position vector
# output
# v1:     velocity at position vector 1
# v2:     velocity at position vector 2
#
def calculate_velocities(h, r1, r2, vr1, vt1, vr2, vt2):
    # calculate magnitudes r1_n and r2_n
    r1_n = np.linalg.norm(r1)
    r2_n = np.linalg.norm(r2)
    # tangential unit vectors
    hXr1   = np.cross(h, r1)
    hXr1_n = np.linalg.norm(hXr1)
    hXr2   = np.cross(h, r2)
    hXr2_n = np.linalg.norm(hXr2)
    # final velocities from components
    v1 = vr1 * (r1/r1_n) + vt1 * ( hXr1 / hXr1_n )
    v2 = vr2 * (r2/r2_n) + vt2 * ( hXr2 / hXr2_n )
    # return velocities
    return v1, v2

# --------------------------------------------
# from appendix C: Fortran-77 Subroutine VLAMB
# --------------------------------------------
# !  Gooding support routine - vlamb(gm,r1,r2,th,tdelt,n,vri,vti,vrf,vtf)
# inputs
# mu:     gravitational parameter (GM)
# r1:     position vector 1
# r2:     position vector 2
# t_sec   time of flight in secs
# outputs
# vr1:    list of radial velocity of the initial position vector
# vt1:    list of tangential velocity of the initial position vector
# vr2:    list of radial velocity of the final position vector
# vt2:    list of tangential velocity of the final position vector
#
def vlamb(mu, r1, r2, th, t_sec): 
    # !the following yields m = 0 when th = 2 pi exactly
    # ! neither this nor the original code works for th < 0.0
    thr2 = th
    m = 0
    while thr2 > 2 * np.pi:
        thr2 = thr2 - 2 * np.pi
        m = m + 1
    thr2 = thr2 / 2.0
    
    # calculate magnitudes r1_n and r2_n
    r1_n = np.linalg.norm(r1)
    r2_n = np.linalg.norm(r2)
    
    dr     = r1_n - r2_n
    r1r2   = r1_n * r2_n
    r1r2th = 4.0 * r1r2 * np.sin(thr2) ** 2
    csq    = dr**2 + r1r2th
    c      = np.sqrt(csq)
    s      = (r1_n + r2_n + c) / 2.0
    gms    = np.sqrt(mu * s / 2.0) 
    qsqfm1 = c / s
    q      = np.sqrt(r1r2) * np.cos(thr2) / s

    if c != 0.0:
        rho = dr / c
        sig = r1r2th / csq
    else:
        rho = 0.0
        sig = 1.0

    t = 4.0 * gms * t_sec / s**2
    
    # CALL XLAMB (M, Q, QSQFMI, T, N, XI, X2)
    # PROCEED FOR SINGLE SOLUTION, OR A PAIR
    n, x1, x2 = xlamb( m, q, qsqfm1, t)
    
    vr1_list, vt1_list, vr2_list, vt2_list = [], [], [], []
    # !proceed for single solution, or a pair
    for i in range(n):
        if i == 0: x = x1
        else: x = x2
        
        # compute radial and tangential velocity components
        others, qzminx, qzplx, zplqx = tlamb(m, q, qsqfm1, x, -1)
        
        vt2 = gms * zplqx * np.sqrt(sig)
        vr1 = gms * (qzminx - qzplx * rho) / r1_n
        vt1 = vt2 / r1_n
        vr2 = -gms * (qzminx + qzplx * rho) / r2_n
        vt2 = vt2 / r2_n
        
        vr1_list.append(vr1)
        vt1_list.append(vt1)
        vr2_list.append(vr2)
        vt2_list.append(vt2)
    # end - for
    
    # no of solutions and radial and tangential velocity components list
    return n, vr1_list, vt1_list, vr2_list, vt2_list


# 8th root function, used by xlamb
def d8rt(x):
    return np.sqrt(np.sqrt(np.sqrt(x)))

# !  Gooding support routine - xlamb(m,q,qsqfm1,tin,n,x,xpl)
# inputs
# m:      no of revolutions.
# q:      transfer angle.
# qsqfm1: (1-q^2)
# tin:    time of flight
# outputs
# n:      np of output parameters
# x:      first solution
# xpl:    second solution
#
def xlamb(m, q, qsqfm1, tin):
    # auxiliary parameters.
    tol = 3.0e-7
    xpl = 0
    c0, c1, c2, c3, c41, c42 = 1.7, 0.5, 0.03, 0.15, 1.0, 0.24
    thr2 = np.arctan2(qsqfm1, 2.0 * q) / np.pi

    # set flag for goto statement
    goto3 = False
    
    # !single-rev starter from t (at x = 0) & bilinear (usually)
    if m == 0:
        n = 1
        # calculate time derivatives.
        t0, dt, d2t, d3t = tlamb(m, q, qsqfm1, 0.0, 0)
        tdiff = tin - t0        
        if tdiff <= 0.0:
            # !(-4 is the value of dt, for x = 0)
            x = t0 * tdiff / (-4.0 * tin)
        else:
            x = -tdiff / (tdiff + 4.0)
            w = x + c0 * np.sqrt(2.0 * (1.0 - thr2))

            if w < 0.0:
                x = x - np.sqrt(d8rt(-w)) * (x + np.sqrt(tdiff / (tdiff + 1.5 * t0)))

            w = 4.0 / (4.0 + tdiff)
            x = x * (1.0 + x * (c1 * w - c2 * x * np.sqrt(w)))
        # end - else
    # end - if
    
    # !with multirevs, first get t(min) as basis for starter
    else:
        xm = 1.0 / (1.5 * (m + 0.5) * np.pi)

        if thr2 < 0.5: xm = d8rt(2.0 * thr2) * xm
        if thr2 > 0.5: xm = (2.0 - d8rt(2.0 - 2.0 * thr2)) * xm
        # !(starter for tmin)
        i = 1
        while i < 12:
            # calculate time derivatives.
            tmin, dt, d2t, d3t = tlamb(m, q, qsqfm1, xm, 3)

            if d2t == 0.0: break

            xmold = xm
            xm = xm - dt * d2t / (d2t * d2t - dt * d3t / 2.0)
            xtest = np.abs(xmold/xm - 1.0)
            if (xtest<=tol): break
            i = i + 1
            
        # !(break off & exit if tmin not located - should never happen)
        if (i>12):
            # !now proceed from t(min) to full starter
            n = -1
            return n, x, None,
        
        tdiffm = tin - tmin
        
        # !(exit if no solution with this m)
        if tdiffm < 0.0:
            raise ValueError("error : no solution with this m")
        elif tdiffm == 0.0:
            x = xm
            n = 1
            # !(exit if unique solution already from x(tmin))
            return n, x, None,
        else:
            n = 3
            if d2t == 0.0:   d2t = 6.0 * m * np.pi

            x = np.sqrt(tdiffm / (d2t / 2.0 + tdiffm / (1.0 - xm) ** 2))
            w = xm + x
            w = w*4.0 / (4.0 + tdiffm) + (1.0 - w)**2
            x = ( x*(1.0 - (1.0 + m + c41*(thr2 - 0.5)) \
                     / (1.0 + c3*m) * x * (c1*w + c2*x*np.sqrt(w)) ) + xm )
            d2t2 = d2t / 2.0

            if x >= 1.0:
                n = 1
                goto3 = True # goto 3 statement in fortran
            # !(no finite solution with x > xm)
        # end - else
    # end - else

    # !(now have a starter, so proceed by halley)    
    while True:          # 5   continue statement     
        if goto3 is False:        
            for i in range(3):
                t, dt, d2t, d3t = tlamb(m, q, qsqfm1, x, 2)           
                t = tin - t
                if dt != 0.0:   x = x + t * dt / (dt * dt + t * d2t / 2.0)
                    
                # !(exit if only one solution, normally when m = 0)
                if n != 3:
                    return n, x, xpl,                
            n = 2
            xpl = x           
        # Update the goto condition
        else:  goto3 = False            
        # !(second multi-rev starter)
        t0, dt, d2t, d3t = tlamb(m, q, qsqfm1, 0.0, 0)
        tdiff0 = t0 - tmin
        tdiff = tin - t0

        if tdiff <= 0.0:
            x = xm - np.sqrt( tdiffm / (d2t2 - tdiffm*(d2t2 / tdiff0 - 1.0 / xm**2)) )
        else:
            x = -tdiff / (tdiff + 4.0)
            w = x + c0*np.sqrt(2.0*(1.0 - thr2))

            if w < 0.0:
                x = x - np.sqrt(d8rt(-w)) * ( x + np.sqrt(tdiff / (tdiff + 1.5 * t0)) )

            w = 4.0 / (4.0 + tdiff)
            x = x * (1.0+(1.0+m+ c42*(thr2-0.5)) / (1.0+c3*m)*x*(c1*w-c2*x*np.sqrt(w)) )

            if x <= -1.0:
                n = n - 1
                # !(no finite solution with x < xm)
                if n == 1:   x = xpl
            # end - if
        # end - else
    # end while - goto 5

    # return n, first sol, second sol
    return n, x, xpl,

# !  Gooding support routine - tlamb(m,q,qsqfm1,x,n,t,dt,d2t,d3t)
# inputs
# m:      no of revolutions
# q:      transfer angle
# qsqfm1: (1-q^2)
# x:      independent variable
# n:      np of output parameters
# outputs
# t:      time at x
# dt:     first derivative of the non-dimensional time evaluated at x
# d2t:    second derivative of the non-dimensional time evaluated at x
# d3t:    third derivative of the non-dimensional time evaluated at x
#
def tlamb(m, q, qsqfm1, x, n):
    # define necessary parameters.
    sw = 0.4
    lm1 = n == -1
    l1 = n >= 1
    l2 = n >= 2
    l3 = n == 3
    qsq = q * q
    xsq = x * x
    u = (1.0 - x) * (1.0 + x)

    # !(needed if series, and otherwise useful when z = 0)
    if not lm1:
        dt, d2t, d3t = 0.0, 0.0, 0.0

    # !direct computation (not series)
    if lm1 or m > 0 or x < 0.0 or np.abs(u) > sw:
        y = np.sqrt(np.abs(u))
        z = np.sqrt(qsqfm1 + qsq * xsq)
        qx = q * x

        if qx <= 0.0:
            a = z - qx
            b = q * z - x
            
        if qx <= 0.0 and lm1:
            aa = qsqfm1 / a
            bb = qsqfm1 * (qsq * u - xsq) / b

        if (qx == 0.0 and lm1) or (qx > 0.0):
            aa = z + qx
            bb = q * z + x

        if qx > 0.0:
            a = qsqfm1 / aa
            b = qsqfm1 * (qsq * u - xsq) / bb
        
        if lm1:
            t, dt, d2t, d3t = 0, b, bb, aa
        else:
            if qx * u >= 0.0:  g = x * z + q * u
            else:              g = (xsq - qsq * u) / (x * z - q * u)

            f = a * y

            if x <= 1.0:       t = m * np.pi + np.arctan2(f, g)
            else:
                if f > sw:     t = np.log(f + g)
                else:
                    # !(continue looping for inverse tanh)
                    fg1 = f / (g + 1.0)
                    term = 2.0 * fg1
                    fg1sq = fg1 * fg1
                    t = term
                    twoi1 = 1.0

                    # do-loop implementation.
                    twoi1 = twoi1 + 2.0
                    term = term * fg1sq
                    told = t
                    t = t + term / twoi1
                    # if (t/=told) cycle else exit
                    while t != told:
                        twoi1 = twoi1 + 2.0
                        term = term * fg1sq
                        told = t
                        t = t + term / twoi1
                    # end - while
                # end - else
            # end - else
            t = 2.0 * (t / y + b) / u

            if l1 and z != 0.0:
                qz = q / z
                qz2 = qz * qz
                qz = qz * qz2
                dt = (3.0*x*t - 4.0*(a + qx * qsqfm1) / z) / u
                if l2:    d2t = (3.0*t + 5.0*x*dt + 4.0*qz*qsqfm1) / u
                if l3:    d3t = ( 8.0*dt + 7.0*x*d2t - 12.0*qz*qz2*x*qsqfm1 ) / u
            # end - if
        # end - else
    else:
        # !compute by series
        u0i = 1.0

        if l1:    u1i = 1.0
        if l2:    u2i = 1.0
        if l3:    u3i = 1.0

        term = 4.0
        tq = q * qsqfm1
        i = 0

        if q < 0.5:   tqsum = 1.0 - q * qsq
        if q >= 0.5:  tqsum = (1.0 / (1.0 + q) + q) * qsqfm1

        ttmold = term / 3.0
        t = ttmold * tqsum

        # do-loop implementation
        i = i + 1
        p = i
        u0i = u0i * u

        if l1 and i > 1:   u1i = u1i*u
        if l2 and i > 2:   u2i = u2i*u
        if l3 and i > 3:   u3i = u3i*u

        term = term * (p - 0.5) / p
        tq = tq * qsq
        tqsum = tqsum + tq
        told = t
        tterm = term / (2.0 * p + 3.0)
        tqterm = tterm * tqsum
        t = t - u0i*( (1.5*p + 0.25)*tqterm / (p*p - 0.25) - ttmold*tq )
        ttmold = tterm
        tqterm = tqterm*p

        if l1:    dt = dt + tqterm*u1i
        if l2:   d2t = d2t + tqterm*u2i*(p - 1.0)
        if l3:   d3t = d3t + tqterm*u3i*(p - 1.0)*(p - 2.0)

        # if (i<n .or. t/=told) cycle else exit.
        while i < n or t != told:
            i = i + 1
            p = i
            u0i = u0i * u

            if l1 and i > 1:   u1i = u1i*u
            if l2 and i > 2:   u2i = u2i*u
            if l3 and i > 3:   u3i = u3i*u

            term = term * (p - 0.5) / p
            tq = tq * qsq
            tqsum = tqsum + tq
            told = t
            tterm = term / (2.0 * p + 3.0)
            tqterm = tterm * tqsum
            t = t - u0i*( (1.5*p + 0.25)*tqterm / (p*p - 0.25) - ttmold*tq )
            ttmold = tterm
            tqterm = tqterm * p

            if l1:   dt = dt + tqterm*u1i
            if l2:  d2t = d2t + tqterm*u2i*(p - 1.0)
            if l3:  d3t = d3t + tqterm*u3i*(p - 1.0)*(p - 2.0)
        # end - while
        if l3:   d3t = 8.0*x*(1.5*d2t - xsq*d3t)
        if l2:   d2t = 2.0*(2.0*xsq*d2t - dt)
        if l1:    dt = -2.0*x*dt

        t = t / xsq
    # return time derivatives
    return t, dt, d2t, d3t


    
   