# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:02:08 2016

@author: Aaron
"""

import kepler as kp
import numpy as np

def gen_object(n):
    return np.zeros((n), dtype=[('a', 'f8'), ('e', 'f8'), ('i', 'f8'), 
                                ('omega', 'f8'), ('w', 'f8'), ('M', 'f8'), 
                                ('mass', 'f8'), ('name', 'S10')])

def generate_solar_system():
    planets = gen_object(5)
    
    planets['name'] = ('p_jupiter', 'p_saturn', 'p_uranus', 'p_neptune', 'p_pluto')
    planets['mass'] = (0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 1./1.3e8)
    planets['a'] = (5.20336301, 9.53707032, 19.19126393, 30.06896348, 39.48168677)
    planets['e'] = (0.04839266, 0.05415060, 0.04716771, 0.00858587, 0.24880766)
    planets['i'] = (1.30530, 2.48446, 0.76986, 1.76917, 17.14175)
    planets['omega'] = (100.55615, 113.71504, 74.22988, 131.72169, 110.30347)
    planets['w'] = (274.19769287, 338.71688843, 96.73435974, 273.24966431, 113.76329041)
    planets['M'] = (19.65054321, 317.51239014, 142.26794434, 259.90869141, 14.86204529)
    
    return planets

def generate_kuiper_objects(n=400):
    kbos = gen_object(n)
    
    for i, kbo in enumerate(kbos):
        kbo['name'] = "a_%d" % i
    
    kbos['mass'] = 0.
    kbos['a'] = np.random.uniform(50, 550, n)
    kbos['e'] = 1. - np.random.uniform(30, 50, n) / kbos['a']
    kbos['i'] = 0.
    kbos['omega'] = np.random.uniform(0., 360, n)
    kbos['w'] = np.random.uniform(0., 360, n)
    kbos['M'] = np.random.uniform(0., 360, n)
    
    return kbos

def generate_perturber(mass=10):
    earth_mass = 3.003e-6
    p = gen_object(1)
    
    p['name'] = "p_new"
    p['mass'] = earth_mass * mass
    p['a'] = 500
    p['e'] = 0.6
    p['i'] = 0.0
    p['omega'] = 0.
    p['w'] = 150
    p['M'] = 0
    
    return p

# Just do a simple CSV, write name, position[3], velocity[3]
def output_objects(objects, filename='objects.csv'):
    with open(filename, 'w') as csv_file:
        # write Sun
        csv_file.write('s_sun,1.00000597682,0.0,0.0,0.0,0.0,0.0,0.0\n')
        for obj in objects:
            r, v = kp.kepler_to_cart(obj['a'], obj['e'], obj['i'], obj['omega'], 
                                     obj['w'], obj['M'], convert_degree=True)
            data = '%s,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f' % \
                (obj['name'], obj['mass'], r[0], r[1], r[2], v[0], v[1], v[2])
            csv_file.write(data + '\n')

if __name__ == '__main__':
    planets = generate_solar_system()
    p = generate_perturber()
    kbos = generate_kuiper_objects(50)
    
    all_objects = np.concatenate((planets, p, kbos))
    
    output_objects(all_objects, r'H:\nbody\objects.csv')
#    output_objects(planets, r'H:\nbody\objects.csv')