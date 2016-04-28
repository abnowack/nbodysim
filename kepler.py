# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 23:47:54 2016

@author: Aaron
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def mag(vec):
    return np.sqrt(np.dot(vec, vec))
    
def kepler_to_cart(a, e, i, O, w, M, mass=0.0, G=39.4876393, convert_degree=False):
    if convert_degree:
        i = np.deg2rad(i)
        O = np.deg2rad(O)
        w = np.deg2rad(w)
        M = np.deg2rad(M)
    
    n = np.sqrt(G * (1 + mass)) * pow(a, -1.5)
    
    niterations = 20
    E = np.ones(np.size(M)) * 0.5
    for itervar in xrange(niterations):
        E -= (E - e * np.sin(E) - M)
    
    q1 = a * (np.cos(E) - e)
    q2 = a * np.sqrt(1 - e * e) * np.sin(E)
    q1dot = -n * a * np.sin(E) / (1 - e * np.cos(E))
    q2dot = n * a * np.sqrt(1 - e*e) * np.cos(E) / (1 - e * np.cos(E))
    
    r = np.zeros((3, np.size(q1)))
    rdot = np.zeros((3, np.size(q1)))
    
    r[0, :] = \
        (np.cos(O) * np.cos(w) - np.sin(O) * np.cos(i) * np.sin(w)) * q1 + \
        (-np.cos(O) * np.sin(w) - np.sin(O) * np.cos(i) * np.cos(w)) * q2
    r[1, :] = \
        (np.sin(O) * np.cos(w) + np.cos(O) * np.cos(i) * np.sin(w)) * q1 + \
        (-np.sin(O) * np.sin(w) + np.cos(O) * np.cos(i) * np.cos(w)) * q2
    r[2, :] = (np.sin(i) * np.sin(w)) * q1 + (np.sin(i) * np.cos(w)) * q2
    
    rdot[0, :] = \
        (np.cos(O) * np.cos(w) - np.sin(O) * np.cos(i) * np.sin(w)) * q1dot + \
        (-np.cos(O) * np.sin(w) - np.sin(O) * np.cos(i) * np.cos(w)) * q2dot
    rdot[1, :] = \
        (np.sin(O) * np.cos(w) + np.cos(O) * np.cos(i) * np.sin(w)) * q1dot + \
        (-np.sin(O) * np.sin(w) + np.cos(O) * np.cos(i) * np.cos(w)) * q2dot
    rdot[2, :] = \
        (np.sin(i) * np.sin(w)) * q1dot + (np.sin(i) * np.cos(w)) * q2dot
    
    return r, rdot

def cart_to_kepler(r, rdot, mass=0.0, G=39.4876393, convert_degree=False): 
    GM = G * (1+mass)
    
    h = np.cross(r, rdot)
    
    e_vec = np.cross(rdot, h) / G - r / mag(r)
    
    n = np.array([-h[1], h[0], 0])
    
    nu = np.dot(e_vec, r) / (mag(e_vec) * mag(r))
    if nu > 1:
        nu = 1
    nu = np.arccos(nu)
    if np.dot(r, rdot) < 0:
        nu = 2 * np.pi - nu
    
    i = np.arccos(h[2] / mag(h))    
    
    e = mag(e_vec)
    
    E = np.tan(nu / 2) / np.sqrt((1 + e) / (1 - e))
    E = 2 * np.arctan(E)
    
    Omega = np.arccos(n[0] / mag(n))
    if n[1] < 0:
        Omega = 2 * np.pi - Omega
    
    w = np.arccos(np.dot(n, e_vec) / mag(e_vec) / mag(n))
    if e_vec[2] < 0:
        w = 2 * np.pi - w
    
    M = E - e * np.sin(E)
    
    a = 2 / mag(r) - mag(rdot) * mag(rdot) / GM
    a = 1 / a
    
    return (a, e, np.rad2deg(i), np.rad2deg(Omega), np.rad2deg(w), M)

def cart_to_kepler2(r, rdot):
    G = 39.4876393 # in AU^3 / M_solar / year^2    
    GM = G * 1
    
    h = np.cross(r, rdot)
    
    e_vec = np.cross(rdot, h) / G - r / mag(r)
    
    n = np.array([-h[1], h[0], 0])
    
    nu = np.dot(e_vec, r) / (mag(e_vec) * mag(r))
    if nu > 1:
        nu = 1
    nu = np.arccos(nu)
    if np.dot(r, rdot) < 0:
        nu = 2 * np.pi - nu
    
    i = np.arccos(h[2] / mag(h))    
    
    e = mag(e_vec)
    
    E = np.tan(nu / 2) / np.sqrt((1 + e) / (1 - e))
    E = 2 * np.arctan(E)
    
    Omega = np.arccos(n[0] / mag(n))
    if n[1] < 0:
        Omega = 2 * np.pi - Omega
    
    w = np.arccos(np.dot(n, e_vec) / mag(e_vec) / mag(n))
    if e_vec[2] < 0:
        w = 2 * np.pi - w
    
    M = E - e * np.sin(E)
    
    a = 2 / mag(r) - mag(rdot) * mag(rdot) / GM
    a = 1 / a
    
    return (a, e, np.rad2deg(i), np.rad2deg(Omega), np.rad2deg(w), M)
    

def orbital_trace(r, rdot, n=50, mass=0.0, G=39.4876393):
    (a, e, i, Omega, w, M) = cart_to_kepler(r, rdot, mass=0.0, G=39.4876393, convert_degree=False)
    M = np.linspace(0., 360, n+1)
    t, tdot = kepler_to_cart(a, e, i, Omega, w, M, convert_degree=True)  
    return t, tdot

if __name__ == '__main__':
    jupiter = {'a': 5.20336301, 'e': 0.04839266, 'i': 1.30530, 'Omega': 100.55615, 'omega': 14.75385 - 100.55615}
    saturn = {'a': 9.53707032, 'e': 0.05415060, 'i': 2.48446, 'Omega': 113.71504, 'omega': 92.43194 - 113.71504}
    uranus = {'a': 19.19126393, 'e': 0.04716771, 'i': 0.76986, 'Omega': 74.22988, 'omega': 170.96424 - 74.22988}
    neptune = {'a': 30.06896348, 'e': 0.00858587, 'i': 1.76917, 'Omega': 131.72169, 'omega': 44.97135 - 131.72169}
    pluto = {'a': 39.48168677, 'e': 0.24880766, 'i': 17.14175, 'Omega': 110.30347, 'omega': 224.06676 - 110.30347}
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot([0], [0], [0], 'o', color='red')
    
    for p in [jupiter, saturn, uranus, neptune, pluto]:
        r, rdot = kepler_to_cart(p['a'], p['e'], p['i'], p['Omega'], p['omega'], np.linspace(0., 2 * np.pi, 51), convert_degree=True)
        ax.plot(r[0], r[1], r[2])
    
    ax.set_xlim3d(-50, 50)
    ax.set_ylim3d(-50, 50)
    ax.set_zlim3d(-50, 50)