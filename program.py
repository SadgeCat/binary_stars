from vpython import *

scene = canvas(width=600, height=600)

#constants
G = 6.67e-11
M0 = 1.989e30       # solar mass
Rs = 6.957e8        # radius of sun

mA = 3.4
mB = 0.8
q = mB/mA           # mass ratio
dist = 1e11

x1 = -dist*mB/(mA+mB)
x2 = dist*mA/(mA+mB)


#lobe_rad = 

# 2 stars for binary star system
starA = sphere(pos = vector(x1, 0, 0), radius = 2.7*Rs, color = color.yellow)
starB = sphere(pos = vector(x2, 0, 0), radius = 3.4*Rs, color = color.blue)

sep_dist = mag(starA.pos - starB.pos)

# Please note that I made the radii of the earth and the Sun much too large, just so they're more visible. 
# All other quantities are realistic.


######################## USER INTERFACES ########################
# button to start/pause simulation
running = False
button(text="Click to Run", pos=scene.title_anchor, bind=Run)
def Run(b):
    global running
    running = not running
    if running: 
        b.text = "Click to Pause"
        print("running")
    else: 
        b.text = "Click to Run"
        print("not running")

# sliders to adjust dist, mA, mB

def update_system():
    global x1, x2, dist, sep_dist, q
    q = mB/mA
    x1 = -dist*mB/(mA + mB)
    x2 = dist*mA/(mA + mB)

    starA.pos.x = x1
    starB.pos.x = x2

    sep_dist = mag(starA.pos - starB.pos)

#######################
# DISTANCE SLIDER
#######################
def changeDistSlider(evt):
    global dist
    dist = evt.value
    dist_text.text = '{:1.2f}'.format(evt.value)
    update_system()
changeDist = slider(bind=changeDistSlider, min=0.5*1e11, max=2*1e11, value=dist)
dist_text = wtext(text='{:1.2f}'.format(changeDist.value))


#######################
# Mass A SLIDER
#######################
def change_mASlider(evt):
    global mA
    mA = evt.value
    mA_text.text = '{:1.2f}'.format(evt.value)
    update_system()
change_mA = slider(bind=change_mASlider, min=0.1, max=5, value=mA)
mA_text = wtext(text='{:1.2f}'.format(change_mA.value))

    
#######################
# Mass B SLIDER
#######################
def change_mBSlider(evt):
    global mB
    mB = evt.value
    mB_text.text = '{:1.2f}'.format(evt.value)
    update_system()
change_mB = slider(bind=change_mBSlider, min=0.1, max=5, value=mB)
mB_text = wtext(text='{:1.2f}'.format(change_mB.value))
################################################################# 

starA.mass = mA * M0
starB.mass = mB * M0

starA.velocity = vector(0,0,0)
starB.velocity = vector(0,0,0)

starA.acc = vector(0,0,0)
starB.acc = vector(0,0,0)


def gravity(star, satellite):
    rad = satellite.pos - star.pos
    return -G*star.mass*satellite.mass*hat(rad)/(mag(rad)**2)
    

def potential(x, y, z):
#    sep_dist = starA.pos - starB.pos
    w_squared = G * starA.mass * (1 + q)/(sep_dist ** 3)
    
    r1 = sqrt((starA.pos.x - x) ** 2 + y ** 2 + z ** 2)
    r2 = sqrt((starB.pos.x - x) ** 2 + y ** 2 + z ** 2)
    r3_squared = x**2 + y**2 + z**2
    
#    W = 1/r1 + 1/r2 + 0.5*(1+q)*x1**2
    W = G*starA.mass/r1 + G*starB.mass/r2 + 0.5*w_squared*r3_squared
    return W
    
# find x value of lagrangian pt

best_x = 0
least_force = 1e100

for i in range(200):
    x = x1 + (x2-x1)*i/200
    r1 = abs(x1-x)
    r2 = abs(x2-x)
    w_squared = G * starA.mass * (1 + q)/(sep_dist ** 3)
    f = abs(G*starA.mass/r1**2 - G*starB.mass/r2**2 - w_squared*x)
    
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
    if running:
        starA.acc = gravity(starB,starA)/starA.mass
        starB.acc = gravity(starA,starB)/starB.mass
        
        starA.velocity = starA.velocity + starA.acc*dt
        starB.velocity = starB.velocity + starB.acc*dt
        
        starA.pos = starA.pos + starA.velocity*dt
        starB.pos = starB.pos + starB.velocity*dt
         
        
        t = t+dt