# -*- coding: utf-8 -*-

# we need to import the time module from astropy
from astropy.time import Time
# import some coordinate things from astropy
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric_posvel
from astropy import units as u
from astropy import constants as c
import matplotlib.pyplot as plt
import numpy as np


def grav(r,const):
    #Calculates the gravitational acceleration using the location of the planet
    #in AU
    #Calculate the distance between the planet and the sun
    dist = (np.sqrt(np.sum((r)**2)))
    #Calculate the gravitational acceleration using F = ma
    a = (-const*r/(dist**3))
    return a

    
def leapfrog(x0,v0,dt, N, const):
    #Determines the location and velocity of an object in orbit using the leapfrog algorithm
    #Kick the algorithm with a Gaussian half step for the velocity
    a = grav(x0, const)
    halfstep = 0.5*dt
    vhalf = v0 + halfstep*grav(x0,const)
    x=np.zeros((N,3))
    v=np.zeros((N,3))
    x[0,:] = x0
    v[0,:] = vhalf
    
    for i in range(1,N):
        #Take N-1 steps of the algorithm
        a = grav(x[i-1,:],const)
        #Alternate x_i and v_i+1/2
        x_i = x[i-1,:] + dt*v[i-1,:]
        v_i = v[i-1,:] + dt*a
        x[i,:] = x_i
        v[i,:] = v_i
        
    return x, v

# pick a time (please use either this or the current time)
t = Time("2021-12-07 10:00")
# initialize the planets;
with solar_system_ephemeris.set('jpl'):
    mars = get_body_barycentric_posvel('mars', t)
    earth = get_body_barycentric_posvel('earth', t)
    sun = get_body_barycentric_posvel('sun', t)
    mercury = get_body_barycentric_posvel('mercury', t)
    venus = get_body_barycentric_posvel('venus', t)
    jupiter = get_body_barycentric_posvel('jupiter', t)
    saturn = get_body_barycentric_posvel('saturn', t)
    neptune = get_body_barycentric_posvel('neptune', t)
    uranus = get_body_barycentric_posvel('uranus', t)


#Cast the planets' positions and velocities to a usable form
marsposition = mars[0].xyz.to_value(u.AU)
marsvelocity = mars[1].xyz.to_value(u.AU/u.d)
earthposition = earth[0].xyz.to_value(u.AU)
earthvelocity = earth[1].xyz.to_value(u.AU/u.d)
sunposition = sun[0].xyz.to_value(u.AU)
#Assume the sun is stationary for simplicity
sunvelocity = [0,0,0]
mercuryposition = mercury[0].xyz.to_value(u.AU)
mercuryvelocity = mercury[1].xyz.to_value(u.AU/u.d)
venusposition = venus[0].xyz.to_value(u.AU)
venusvelocity = venus[1].xyz.to_value(u.AU/u.d)
jupiterposition = jupiter[0].xyz.to_value(u.AU)
jupitervelocity = jupiter[1].xyz.to_value(u.AU/u.d)
saturnposition = saturn[0].xyz.to_value(u.AU)
saturnvelocity = saturn[1].xyz.to_value(u.AU/u.d)
neptuneposition = neptune[0].xyz.to_value(u.AU)
neptunevelocity = neptune[1].xyz.to_value(u.AU/u.d)
uranusposition = uranus[0].xyz.to_value(u.AU)
uranusvelocity = uranus[1].xyz.to_value(u.AU/u.d)

positions = np.array([sunposition,mercuryposition,venusposition,earthposition,
                      marsposition, jupiterposition,saturnposition,uranusposition,neptuneposition])
velocities = np.array([sunvelocity,mercuryvelocity,venusvelocity,earthvelocity,marsvelocity,
                       jupitervelocity,saturnvelocity,uranusvelocity,neptunevelocity])
names = np.array(['Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune'])

#Create a plot of the initial conditions
fig, axs = plt.subplots(1,2)
for i in range(len(names)):
    axs[0].plot(positions[i,0],positions[i,1],label=names[i],marker='o')
    axs[1].plot(positions[i,0],positions[i,2],marker='o')
axs[0].set_xlabel('x-position (AU)')
axs[0].set_ylabel('y-position (AU)')
axs[1].set_xlabel('x-position (AU)')
axs[1].set_ylabel('z-position (AU)')
fig.suptitle('Initial positions of the solar system')
handles, labels = axs[0].get_legend_handles_labels()
fig.tight_layout()
fig.subplots_adjust(top=0.8)
fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.05,0.95),ncols=int(len(names)/2+1))
fig.savefig('./plots/initials.pdf')
plt.close()


iterations = int(200*365.25*2)
time = np.arange(0,200*365.25,0.5)
#Create a plot of the orbits
fig, axs = plt.subplots(2,1,height_ratios=[2,1])
constants = (c.G*c.M_sun).to(u.AU**3/u.d**2).value
for i in range(1,len(names)):
    x,v = leapfrog(positions[i,:],velocities[i,:],0.5,iterations,constants)
    axs[0].plot(x[:,0],x[:,1],label=names[i])
    axs[1].plot(time,x[:,2])
axs[0].set_xlabel('x-position (AU)')
axs[0].set_ylabel('y-position (AU)')
axs[1].set_xlabel('time (days)')
axs[1].set_ylabel('z-position (AU)')
fig.suptitle('Simulated orbits of the solar system')
fig.tight_layout()
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right')
fig.savefig('./plots/orbits.pdf')
plt.close()