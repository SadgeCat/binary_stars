Web VPython 3.2

scene = canvas(width=600, height=600)

#constants
G = 6.67e-11
M0 = 1.989e30       # solar mass
Rs = 6.957e8        # radius of sun

mA = 3.4
mB = 0.8
q = mB/mA           # mass ratio
q2 = 1/q
dist = 1e11

x1 = -dist*mB/(mA+mB)
x2 = dist*mA/(mA+mB)

#k = 1.6/(2.7*Rs)
k = 1e-10


# 2 stars for binary star system
starA = sphere(pos = vector(x1, 0, 0), radius = 2.7*Rs, color = color.yellow)
starB = sphere(pos = vector(x2, 0, 0), radius = 3.4*Rs, color = color.blue)

sep_dist = mag(starA.pos - starB.pos)

lobe_rad = sep_dist * (0.49 * q2 ** .6666667) / (0.6 * q2 ** .6666667 + log(1 + q2 ** .3333333))
lobe_rad2 = sep_dist * (0.49 * q ** .6666667) / (0.6 * q ** .6666667 + log(1 + q ** .3333333))

# make star overflow
starA.radius = lobe_rad * (1 + 1e-2)


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
        print(scene.range)
    else: 
        b.text = "Click to Run"
        print("not running")
        
button(text="Reset", pos=scene.title_anchor, bind=reset)
scene.append_to_title("\n\n")
def reset():
    global mA, mB, x1, x2, dist, sep_dist, q, q2, lobe_rad, lobe_rad2, sum_mass, reduced_mass, momentum
    mA = 3.4
    mB = 0.8
    q = mB/mA
    q2 = 1/q
    x1 = -dist*mB/(mA + mB)
    x2 = dist*mA/(mA + mB)

    starA.pos = vector(x1, 0, 0)
    starB.pos = vector(x2, 0, 0)

    sep_dist = mag(starA.pos - starB.pos)
    lobe_rad = sep_dist * (0.49 * q2 ** .6666667) / (0.6 * q2 ** .6666667 + log(1 + q2 ** .3333333))
    lobe_rad2 = sep_dist * (0.49 * q ** .6666667) / (0.6 * q ** .6666667 + log(1 + q ** .3333333))
    
    starA.mass = mA * M0
    starB.mass = mB * M0
    
    starA.radius = lobe_rad * (1 + 1e-2)
    starB.radius = radius_from_mass(starB.mass, C_b, k)
    
    starA.velocity = (sqrt(G * starB.mass * abs(starA.pos.x))/sep_dist)*vector(0, -1, 0)
    starB.velocity = (sqrt(G * starA.mass * abs(starB.pos.x))/sep_dist)*vector(0, 1, 0)
    
    starA.acc = vector(0,0,0)
    starB.acc = vector(0,0,0)
    
    sum_mass = starA.mass + starB.mass
    reduced_mass = starA.mass * starB.mass / sum_mass
    
    momentum = reduced_mass * sep_dist * sqrt(G * sum_mass / sep_dist)
#    changeDistSlider(sep_dist)
#    change_mASlider(mA)
#    change_mBSlider(mB)
    dist_text.text = f"{sep_dist:.3e} m"
    mA_text.text = f"{mA:.3f} M☉"
    mB_text.text = f"{mB:.3f} M☉"
    
    draw_potential()
    
    q_text.text = f"{q:.3f}"
    roche_text.text = f"{lobe_rad:.3e} m"
    
    scene.autoscale = False
    scene.range = scale
    scene.autoscale = True
        

# sliders to adjust dist, mA, mB
def update_system():
    global x1, x2, dist, sep_dist, q, q2, lobe_rad, lobe_rad2, sum_mass, reduced_mass, momentum
    q = mB/mA
    q2 = 1/q
    x1 = -dist*mB/(mA + mB)
    x2 = dist*mA/(mA + mB)

    starA.pos.x = x1
    starB.pos.x = x2

    sep_dist = mag(starA.pos - starB.pos)
    lobe_rad = sep_dist * (0.49 * q2 ** .6666667) / (0.6 * q2 ** .6666667 + log(1 + q2 ** .3333333))
    lobe_rad2 = sep_dist * (0.49 * q ** .6666667) / (0.6 * q ** .6666667 + log(1 + q ** .3333333))
    
    starA.velocity = (sqrt(G * starB.mass * abs(starA.pos.x))/sep_dist)*vector(0, -1, 0)
    starB.velocity = (sqrt(G * starA.mass * abs(starB.pos.x))/sep_dist)*vector(0, 1, 0)
    
    sum_mass = starA.mass + starB.mass
    reduced_mass = starA.mass * starB.mass / sum_mass
    momentum = reduced_mass * sep_dist * sqrt(G * sum_mass / sep_dist)
    
    draw_potential()
    
    q_text.text = f"{q:.3f}"
    roche_text.text = f"{lobe_rad:.3e} m"
    
    
#######################
# DISPLAY INFO
######################
scene.append_to_caption("\nMass Ratio q = ")
q_text = wtext(text=f"{q:.3f}")
scene.append_to_caption("\nRobe Lobe Radius = ")
roche_text = wtext(text=f"{lobe_rad:.3e} m")


#######################
# DISTANCE SLIDER
#######################
scene.append_to_caption("\n\n")
scene.append_to_caption("Separation Dist. (m): ")
scene.append_to_caption("\n")
def changeDistSlider(evt):
    global dist
    dist = evt.value
    dist_text.text = f"{evt.value:.3e} m"
    update_system()
changeDist = slider(bind=changeDistSlider, min=0.5*1e11, max=2*1e11, value=dist, length=300)
dist_text = wtext(text=f"{dist:.3e} m")


#######################
# Mass A SLIDER
#######################
scene.append_to_caption("\n\n")
scene.append_to_caption("Star A Mass (Solar Masses): ")
scene.append_to_caption("\n")
def change_mASlider(evt):
    global mA, C_a
    mA = evt.value
    starA.mass = mA * M0
    C_a = starA.mass / mass_from_radius(starA.radius, 1, k)
    mA_text.text = f"{evt.value:.3f} M☉"
    update_system()
change_mA = slider(bind=change_mASlider, min=0.1, max=5, value=mA, length=300)
mA_text = wtext(text=f"{mA:.3f} M☉")

    
#######################
# Mass B SLIDER
#######################
scene.append_to_caption("\n\n")
scene.append_to_caption("Star B Mass (Solar Masses): ")
scene.append_to_caption("\n")
def change_mBSlider(evt):
    global mB, C_b
    mB = evt.value
    starB.mass = mB * M0
    C_b = starB.mass / mass_from_radius(starB.radius, 1, k)
    mB_text.text = f"{evt.value:.3f} M☉"
    update_system()
change_mB = slider(bind=change_mBSlider, min=0.1, max=5, value=mB, length=300)
mB_text = wtext(text=f"{mB:.3f} M☉")


#######################
# Redraw potential btn 
#######################
scene.append_to_caption("\n\n")
button(text="Draw Equipotential", bind=draw_potential)



scene.append_to_caption("\n\n\n")
################################################################# 

starA.mass = mA * M0
starB.mass = mB * M0

starA.velocity = (sqrt(G * starB.mass * abs(starA.pos.x))/sep_dist)*vector(0, -1, 0)
starB.velocity = (sqrt(G * starA.mass * abs(starB.pos.x))/sep_dist)*vector(0, 1, 0)

starA.acc = vector(0,0,0)
starB.acc = vector(0,0,0)

sum_mass = starA.mass + starB.mass
reduced_mass = starA.mass * starB.mass / sum_mass

momentum = reduced_mass * sep_dist * sqrt(G * sum_mass / sep_dist)
#momentum = starA.mass * sep_dist * mag(starA.velocity)


def mass_from_radius(R,C,k):
    term = 2/k**3 - exp(-k*R) * (R**2/k + 2*R/k**2 + 2/k**3)
    return 4*pi*C*term
    
C_a = starA.mass / mass_from_radius(starA.radius, 1, k)
C_b = starB.mass / mass_from_radius(starB.radius, 1, k)

def radius_from_mass(new_mass, C, k):
    left = 0
    right = 1e11
    for i in range(50):
        mid = (left + right)/2
        if mass_from_radius(mid, C, k) > new_mass:
            right = mid
        else:
            left = mid
    return (left + right)/2
    

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
    
# find x value of lagrange pt
def find_x():
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
    return best_x
        
equipotential = potential(find_x(), 0, 0)
    
# calculate points close enough to the value of equipotential
spheres_list = []
def draw_potential(b):
    global spheres_list
    equipotential = potential(find_x(), 0, 0)
    for sph in spheres_list:
        sph.visible = False
    spheres_list = []
    grid_size = 2 * dist
    step = 2e9
    pts_list = []
    for x in arange(-grid_size, grid_size, step):
        for y in arange(-grid_size, grid_size, step):
            W = potential(x, y, 0)
            if abs(W - equipotential) < abs(equipotential)*1e-2:
                pts_list.append(vector(x,y,0))
    
#    p_curve = curve(pos=pts_list, color=color.cyan, radius=0.5*Rs)

    for p in pts_list:
        s = sphere(pos=p, radius=Rs, color=color.white, opacity=0.9)
        spheres_list.append(s)
        
#    print("done drawing equipotential")
    
draw_potential()
scale = scene.range
    
t=0; dt=3600
# custom inc in rad for now
rad_inc_rate = 5e-6 * Rs
transfer_rate = 1e-7

graph_w = 480
graph_h = 360
vel_graph = graph(title='Velocity over time', xtitle='Time', ytitle='Velocity', width=graph_w, height=graph_h, align='left')
av_graph = gcurve(color=color.red, label='StarA Velocity')
bv_graph = gcurve(color=color.blue, label='StarB Velocity')

rad_graph = graph(title='Radius over time', xtitle='Time', ytitle='Radius', width=graph_w, height=graph_h, align='left')
ar_graph = gcurve(color=color.red, label='StarA Radius')
br_graph = gcurve(color=color.blue, label='StarB Radius')

mass_graph = graph(title='Mass over time', xtitle='Time', ytitle='Mass', width=graph_w, height=graph_h, align='left')
am_graph = gcurve(color=color.red, label='StarA Mass')
bm_graph = gcurve(color=color.blue, label='StarB Mass')

type = "detached"

while((starA.pos-starB.pos).mag>(starA.radius+starB.radius)):
    rate(1000)
    if running:
        starA.acc = gravity(starB,starA)/starA.mass
        starB.acc = gravity(starA,starB)/starB.mass
        
        starA.velocity = starA.velocity + starA.acc*dt
        starB.velocity = starB.velocity + starB.acc*dt
        
        starA.pos = starA.pos + starA.velocity*dt
        starB.pos = starB.pos + starB.velocity*dt
        
        P = starA.mass*starA.velocity + starB.mass*starB.velocity
        
        A_overflow = starA.radius >= lobe_rad
        B_overflow = starB.radius >= lobe_rad2
        if A_overflow and B_overflow:
            type = "contact"
        elif A_overflow or B_overflow:
            type = "semi-detached"
        else:
            type = "detached"
        
        if type == "semi-detached":
            if A_overflow:
                overflow = starA.radius - lobe_rad
#                dm = 5e23 * (overflow/lobe_rad) ** .3333333 * dt
                dm = mass_from_radius(starA.radius, C_a, k) - mass_from_radius(lobe_rad, C_a, k)
                starA.mass -= dm * transfer_rate * dt
                starB.mass += dm * transfer_rate * dt
            elif B_overflow:
                overflow = starB.radius - lobe_rad2
#                dm = 5e23 * (overflow/lobe_rad2) ** .3333333 * dt
                dm = mass_from_radius(starB.radius, C_b, k) - mass_from_radius(lobe_rad2, C_b, k)
                starB.mass -= dm * transfer_rate * dt
                starA.mass += dm * transfer_rate * dt
            
            starA.radius = radius_from_mass(starA.mass, C_a, k)
            starB.radius = radius_from_mass(starB.mass, C_b, k)
            
            # preserve linear momentum so COM doesn't move
            P = starA.mass*starA.velocity + starB.mass*starB.velocity
            v_cm = P/sum_mass
            starA.velocity -= v_cm
            starB.velocity -= v_cm
        
            # updating variables
            q = starB.mass/starA.mass
            q2 = 1/q
            reduced_mass = starA.mass * starB.mass / sum_mass
            
            # using conservation of angular momentum to calc vel
#            sep_dist = mag(starB.pos - starA.pos)
#            v_rel = momentum / (reduced_mass * sep_dist)
#            r_hat = hat(starB.pos - starA.pos)
#            t_hat = vector(-r_hat.y, r_hat.x, 0)
#            vA = v_rel * starB.mass / (starA.mass + starB.mass)
#            vB = v_rel * starA.mass / (starA.mass + starB.mass)
#            starA.velocity = -vA * t_hat
#            starB.velocity =  vB * t_hat
            
            # not much mass is transferred and stop really fast
#            sep_dist = mag(starA.pos - starB.pos)

            # change separation dist using conservation of momentum but only reflected in lobe_rad not actual dist between stars
            sep_dist = (momentum ** 2) / (reduced_mass ** 2 * G * sum_mass)
#            angle1 = atan2(starA.pos.y / starA.pos.x)
#            angle2 = atan2(starB.pos.y / starB.pos.x)
#            x1 = -sep_dist * starB.mass / sum_mass
#            x2 = sep_dist * starA.mass / sum_mass
#            newpos1 = vector(x1 * cos(angle1), x1 * sin(angle1), 0)
#            newpos2 = vector(x2 * cos(angle2), x2 * sin(angle2), 0)
#            starA.pos = newpos1
#            starB.pos = newpos2

            lobe_rad = sep_dist * (0.49 * q2 ** .6666667) / (0.6 * q2 ** .6666667 + log(1 + q2 ** .3333333))
            lobe_rad2 = sep_dist * (0.49 * q ** .6666667) / (0.6 * q ** .6666667 + log(1 + q ** .3333333))
            
        if t % (3600 * 1000) == 0:
            print("radius of star B: " + starB.radius)
            print("total momentum: " + mag(P))
        if t % (3600 * 10) == 0:
            draw_potential()
        
        q_text.text = f"{q:.3f}"
        roche_text.text = f"{lobe_rad:.3e} m"
        dist_text.text = f"{sep_dist:.3e} m"
        mA_text.text = f"{starA.mass/M0:.3f} M☉"
        mB_text.text = f"{starB.mass/M0:.3f} M☉"
        
        vel_graph.select()
        av_graph.plot(t, mag(starA.velocity))
        bv_graph.plot(t, mag(starB.velocity))
        
        rad_graph.select()
        ar_graph.plot(t, starA.radius)
        br_graph.plot(t, starB.radius)
        
        mass_graph.select()
        am_graph.plot(t, starA.mass)
        bm_graph.plot(t, starB.mass)
         
        
        t = t+dt
