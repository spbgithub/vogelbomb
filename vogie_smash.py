# -*- coding: utf-8 -*-
"""
Sean Boyd did this shit
Do what you want with it, but give him a little credit at least. 
"""

import math
import numpy as np
import matplotlib.pyplot as mp


###############################################################
#Needed constants and utility functions
C_drag   = 0.4     #          coefficient of drag
C_lift   = 0.30    #          coefficient of lift (Magnus force)
rho      = 1.23    #kg/m^3  - assumed density of air at sea level
A        = 0.00426 #m^2     - cross sectional area of baseball
r        = A/(2*math.pi) #radius of baseball
romega   = 32.1
tau      = 25.0

m        = 0.145   #kg      - mass of baseball
g        = 9.8     #m/s^2   - acceleration of gravity
dt       = 0.01    #seconds - increment of time used in RK

                   #coefficient used to scale abs(v)v for drag computation

                    #coefficient used for lift in Magnus force
startx   = 0.0
starty   = 2.0/3.0

###############################################################
v_tol    = 0.001   #target difference in velocities

#takes meters per second, converts to miles per hour
def mps_to_mph(mps):
    return 2.236936 * mps

#mph to meters per second
def mph_to_mps(mph):
    return mph/2.236936

#miles to meters
def miles_to_meters(miles):
    return 1609.344 * miles

#meters to miles
def meters_to_miles(meters):
    return meters/1609.344

#degrees to radians
def dtr(deg):
    return math.pi * deg /180.

#radians to degrees
def rtd(rad):
    return 180. * rad / math.pi

def m2f(meters):
    return 3.28084 * meters

def f2m(feet):
    return feet/3.28084

targetx = f2m(375.0)
targety = f2m(75.0)

#################################################################
    
#the vectors we perform Runge-Kutta on are four-vectors (x,y,x',y'),
#where x',y' are the time derivatives of x,y respectively.
#this is needed to cast our differential equation into first order form,
#which is required for Runge-Kutta
#
#we accomplish this by using numpy arrays
#here, we create some utility functions to achieve this
#
#we will do our computations in metric, convert to and from english
#for presentation purposes

def dist(rk_vec):
    return math.sqrt(rk_vec[0]**2 + rk_vec[1]**2)
    
def speed(rk_vec):
    return math.sqrt(rk_vec[2]**2 + rk_vec[3]**2)
    
def direction(rk_vec):
    s = speed(rk_vec)
    return np.array([0.0, 0.0, rk_vec[2] / s, rk_vec[3] / s])
                   
#this is as in the first order ODE y' = f(t, y)
def f(t, rk_vec):
    kappa    = 0.5 * c_d(S(t, speed(rk_vec))) * rho * A / m
    c = -kappa * speed(rk_vec)
    return np.array([rk_vec[2], rk_vec[3], c * rk_vec[2], c * rk_vec[3] - g])

def f2(t, rk_vec):
    sp       = speed(rk_vec)
    kappa    = 0.5 * c_d(S(t, m2f(sp))) * rho * A / m
    c        = -kappa * sp
    mu       = 0.5 * c_l(S(t, m2f(sp))) * rho * A / m    
    d        = mu * sp
    return np.array([rk_vec[2], rk_vec[3], c * rk_vec[2] - d * rk_vec[3], c * rk_vec[3] + d * rk_vec[2] - g])

def S(t, v):
    return romega/v * math.exp(-v*t/(146.7*tau))

#computes the drag coefficient, given the spin ratio
def c_d(s):
    return 0.385*(1.0 + 0.2017 * s**2)

def c_l(s):
    return s/(2.32*s + 0.4)

#this is a quick and dirty implementation of the standard 4th order
#explicit Runge-Kutta with a twist; we proceed until a specified condition on 
#the first coordinate of x is achieved
#this corresponds to seeing how high the baseball is at 375 feet from home plate
def rungekutta_with_cond(start_pt, start_t, h, fn, cutoff_x):
    cur_pt = start_pt
    cur_t  = start_t
    datalist = [cur_pt]
    while cur_pt[0] < cutoff_x:
        k1 = h * fn(cur_t,       cur_pt)
        k2 = h * fn(cur_t + h/2, cur_pt + 0.5 * k1)
        k3 = h * fn(cur_t + h/2, cur_pt + 0.5 * k2)
        k4 = h * fn(cur_t + h,   cur_pt + k3)
        
        cur_pt = cur_pt + (k1 + 2 * k2 + 2 * k3 + k4)/6
        datalist.append(cur_pt)
        cur_t = cur_t + h
    return datalist

#compute the points for the graph of the relation between exit velocity and launch angle
def es_to_langle(fn, dt, tx, larray):
    ploty        = []
    for l in larray:
        #start with ridiculously high exit speed, step down until we've found
        v_cur = int(mph_to_mps(200.0))+1
        x     = np.array([startx, starty, v_cur * math.cos(dtr(l)), v_cur * math.sin(dtr(l))])
        rval  = rungekutta_with_cond(x, 0, dt, fn, tx)
        y     = rval[-1][1]

        while y > targety and y > 0:
            v_cur = v_cur - 2
            x     = np.array([startx, starty, v_cur * math.cos(dtr(l)), v_cur * math.sin(dtr(l))])
            rval  = rungekutta_with_cond(x, 0, dt, fn, tx)
            y     = rval[-1][1]

        #we've found our initial bracket - now we use bisection method to narrow down on velocity
        #needed to hit the target
        if y > 0:
            v_lo = v_cur
            v_hi = v_cur + 5

            while v_hi - v_lo > v_tol:
                v_cur  = 0.5 * (v_hi + v_lo)
                x      = np.array([startx, starty, v_cur * math.cos(dtr(l)), v_cur * math.sin(dtr(l))])
                rval   = rungekutta_with_cond(x, 0, dt, fn, tx)
                y      = rval[-1][1]
                if y < targety:
                    v_lo = v_cur
                else:
                    v_hi = v_cur

            ploty.append(mps_to_mph(v_cur))

        else:
            ploty.append(0)

    return ploty



if __name__ == "__main__":
    #l = 40.0
    #v = 56.0
    #x = np.array([0, 1, v * math.cos(dtr(l)), v * math.sin(dtr(l))])
    #r = rungekutta_with_cond(x, 0, dt, f, targetx)
    #print(vec_m2f(r[1][-1]))
    
    launch_angle = [20 + i for i in range(0,20)]
    
    print("Begin, using Magnus force computation")
    ploty2 = es_to_langle(f2, dt, targetx, launch_angle)
    mp.plot(launch_angle, ploty2)
    mp.xlabel("Launch angle (degrees)")
    mp.ylabel("Exit velocity (mph)")
    mp.show()
    
    print(ploty2)
    val, idx = min((val, idx) for (idx, val) in enumerate(ploty2))
    
    l = launch_angle[idx]
    v = mph_to_mps(ploty2[idx])
    print(l, v)
    u    = np.array([startx, starty, v * math.cos(dtr(l)), v * math.sin(dtr(l))])
    rval = rungekutta_with_cond(u, 0, dt, f2, f2m(600))
    
    x = []
    y = []
    
    for r in rval:
        if r[1] < 0: break
        x.append(m2f(r[0]))
        y.append(m2f(r[1]))
        
    mp.plot(x,y)
    mp.xlabel("Distance from home plate (feet)")
    mp.ylabel("Height (feet)")
    mp.show()
    print(x[-1])
    