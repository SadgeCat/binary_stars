from vpython import *

scene = canvas(width=600, height=600)

#constants
G = 6.67e-11
M0 = 1.989e30       # solar mass
Rs = 6.957e8        # radius of sun

mA = 3.4*M0
mB = 0.8*M0
q = mB/mA           # mass ratio

x1 = mB/(mA+mB)
x2 = 1-x1

# 2 stars for binary star system
starA = sphere(pos = vector(-(1.5e11), 0, 0), radius = 2.7*Rs, color = color.yellow)
starB = sphere(pos = vector(1.5e11, 0, 0), radius = 3.4*Rs, color = color.blue)

# Please note that I made the radii of the earth and the Sun much too large, just so they're more visible. 
# All other quantities are realistic.

starA.mass = mA
starB.mass = mB

starA.velocity = vector(0,0,0)
starB.velocity = vector(0,0,0)

starA.acc = vector(0,0,0)
starB.acc = vector(0,0,0)


def gravity(star, satellite):
    rad = satellite.pos - star.pos
    return -G*star.mass*satellite.mass*hat(rad)/(mag(rad)**2)
    
def potential(cordx, cordy):
    sep_dist = starA.pos - starB.pos
    w_squared = G * starA.mass * (1 + q)/(sep_dist ** 3)
    
    

t=0; dt=3600
    
while((starA.pos-starB.pos).mag>(starA.radius+starB.radius)):
    rate(1000)
  
    starA.acc = gravity(starB,starA)/starA.mass
    starB.acc = gravity(starA,starB)/starB.mass
    
    starA.velocity = starA.velocity + starA.acc*dt
    starB.velocity = starB.velocity + starB.acc*dt
    
    starA.pos = starA.pos + starA.velocity*dt
    starB.pos = starB.pos + starB.velocity*dt
     
    
    t = t+dt