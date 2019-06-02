#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 18:19:07 2019

@author: seanp
"""
import math
import matplotlib.pyplot as mp


###############################################################
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

g        = 9.8
startx   = 0.0
starty   = 2.0/3.0
targetx  = f2m(375.0)
targety  = f2m(75.0)


def f(l):
    return math.sqrt(g * targetx**2/(2*math.cos(l)**2*(targetx*math.tan(l) - targety + starty)))

ls = [20 + 0.5*i for i in range(0, 29)]
v = [mps_to_mph(f(dtr(l))) for l in ls]
mp.plot(ls, v)
mp.xlabel("Launch angle (degrees)")
mp.ylabel("Exit velocity (mph)")
mp.show()


