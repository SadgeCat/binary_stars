from vpython import *
scene = canvas(width=600, height=600)

#constants
G = 6.67e-11
M0 = 1.989e30       # solar mass
Rs = 6.957e8        # radius of sun

# 2 stars for binary star system
starA = sphere(pos = vector(-(1.5e11), 0, 0), radius = 2.7*Rs, color = color.yellow)
starB = sphere(pos = vector(1.5e11, 0, 0), radius = 3.4*Rs, color = color.blue)

starA.mass = 3.4*M0
starB.mass = 0.8*M0
q = starB.mass/starA.mass       # mass ratio

# Please note that I made the radii of the earth and the Sun much too large, just so they're more visible. 
# All other quantities are realistic.

def gravity(star, satellite):
    rad = satellite.pos - star.pos
    return -G*star.mass*satellite.mass*hat(rad)/(mag(rad)**2)
    
def potential(cordx, cordy):
    sep_dist = starA.pos - starB.pos
    w_squared = G * starA.mass * (1 + q)/(sep_dist ** 3)