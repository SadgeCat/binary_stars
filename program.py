from vpython import *

scene = canvas(width=600, height=600)

#constants
G = 6.67e-11
M0 = 1.989e30       # solar mass
Rs = 6.957e8        # radius of sun

mA = 3.4*M0
mB = 0.8*M0
q = mB/mA           # mass ratio
dist = 1e11

x1 = -dist*mB/(mA+mB)
x2 = dist*mA/(mA+mB)

#lobe_rad = 

# 2 stars for binary star system
starA = sphere(pos = vector(x1, 0, 0), radius = 2.7*Rs, color = color.yellow)
starB = sphere(pos = vector(x2, 0, 0), radius = 3.4*Rs, color = color.blue)

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
    

sep_dist = mag(starA.pos - starB.pos)
def potential(x, y, z):
#    sep_dist = starA.pos - starB.pos
    w_squared = G * starA.mass * (1 + q)/(sep_dist ** 3)
    
    r1 = sqrt((starA.pos.x - x) ** 2 + y ** 2 + z ** 2)
    r2 = sqrt((starB.pos.x - x) ** 2 + y ** 2 + z ** 2)
    r3_squared = x**2 + y**2 + z**2
    
#    W = 1/r1 + 1/r2 + 0.5*(1+q)*x1**2
    W = G*mA/r1 + G*mB/r2 + 0.5*w_squared*r3_squared
    return W
    
# find x value of lagrangian pt

best_x = 0
least_force = 1e100

for i in range(200):
    x = x1 + (x2-x1)*i/200
    r1 = abs(x1-x)
    r2 = abs(x2-x)
    w_squared = G * starA.mass * (1 + q)/(sep_dist ** 3)
    f = abs(G*mA/r1**2 - G*mB/r2**2 - w_squared*x)
    
    if f<least_force:
        least_force = f
        best_x = x
        
equipotential = potential(best_x, 0, 0)
    
# calculate points close enough to the value of equipotential
grid_size = 2 * dist
step = 2e9
pts_list = []
for x in arange(-grid_size, grid_size, step):
    for y in arange(-grid_size, grid_size, step):
        W = potential(x, y, 0)
        if abs(W - equipotential) < abs(equipotential)*1e-2:
            pts_list.append(vector(x,y,0))

for p in pts_list:
    print(p)
    sphere(pos = p, radius = 1.2e9, color = color.white)

print("done drawing equipotential")
    
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
