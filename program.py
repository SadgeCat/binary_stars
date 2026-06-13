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


bkg = box(pos = vector(0, 0, -5*Rs), size=vec(1000*Rs, 1000*Rs, 1), texture="https://i.imgur.com/UiEicbd.jpeg")
bkg.visible = True


# 2 stars for binary star system
starA = sphere(pos = vector(x1, 0, 0), radius = 2.7*Rs, color = color.white, texture="https://i.imgur.com/sDNTWNi.jpeg")
starB = sphere(pos = vector(x2, 0, 0), radius = 3.4*Rs, color = color.white, texture="https://i.imgur.com/y0MpbvN.png")

mergedStar = sphere(pos = vector(0, 0, 0), radius = (starA.radius + starB.radius) * 0.95, color = color.white, texture="https://i.imgur.com/mS27qdf.jpeg")
mergedStar.visible = False

sep_dist = mag(starA.pos - starB.pos)

lobe_rad = sep_dist * (0.49 * q2 ** .6666667) / (0.6 * q2 ** .6666667 + log(1 + q2 ** .3333333))
lobe_rad2 = sep_dist * (0.49 * q ** .6666667) / (0.6 * q ** .6666667 + log(1 + q ** .3333333))

# make star overflow
def setA_radius():
    starA.radius = lobe_rad * (1 + 1e-1)
    update_system()
def setB_radius():
    starB.radius = lobe_rad2 * (1 + 1e-1)
    update_system()
    
#starA.radius = lobe_rad * (1 + 1e-2)
#starB.radius = lobe_rad2 * (1 + 1e-2)

#test = sphere(pos = vector(0, 0, 0), radius = lobe_rad, color = color.yellow)


# Please note that I made the radii of the earth and the Sun much too large, just so they're more visible. 
# All other quantities are realistic.


######################## USER INTERFACES ########################
# button to start/pause simulation
running = False
button(text="Click to Run", pos=scene.title_anchor, bind=Run)
def Run(b):
    global running, has_reset
    running = not running
    if running: 
        b.text = "Click to Pause"
        change_mA.disabled = True
        change_mB.disabled = True
        changeDist.disabled = True
        A_over.disabled = True
        B_over.disabled = True
#        print("running")
#        print(scene.range)
    else: 
        b.text = "Click to Run"
#        print("not running")
#        if has_reset:
#            change_mA.disabled = False
#            change_mB.disabled = False
#            changeDist.disabled = False
#            A_over.disabled = False
#            B_over.disabled = False
#            has_reset = False

has_reset = False
button(text="Reset", pos=scene.title_anchor, bind=reset)
scene.append_to_title("\n\n")
def reset():
    global mA, mB, x1, x2, dist, sep_dist, q, q2, lobe_rad, lobe_rad2, sum_mass, reduced_mass, momentum, C_a, C_b, keep_running, has_reset, t
    t = 0
    has_reset = True
    if not running:
        change_mA.disabled = False
        change_mB.disabled = False
        changeDist.disabled = False
        A_over.disabled = False
        B_over.disabled = False
        has_reset = False
    if not drawP:
        change_draw()
    keep_running = True
    starA.visible = True
    starB.visible = True
    mergedStar.visible = False
    mA = 3.4
    mB = 0.8
    q = mB/mA
    q2 = 1/q
    dist = 1e11
    x1 = -dist*mB/(mA + mB)
    x2 = dist*mA/(mA + mB)

    starA.pos = vector(x1, 0, 0)
    starB.pos = vector(x2, 0, 0)

    sep_dist = mag(starA.pos - starB.pos)
    lobe_rad = sep_dist * (0.49 * q2 ** .6666667) / (0.6 * q2 ** .6666667 + log(1 + q2 ** .3333333))
    lobe_rad2 = sep_dist * (0.49 * q ** .6666667) / (0.6 * q ** .6666667 + log(1 + q ** .3333333))
    
    starA.mass = mA * M0
    starB.mass = mB * M0
    
    starA.radius = 2.7 * Rs
    starB.radius = 3.4 * Rs
#    starB.radius = radius_from_mass(starB.mass, C_b, k)
    
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
    dist_text.text = f"{dist/Rs:.3f} R☉"
    mA_text.text = f"{mA:.3f} M☉"
    mB_text.text = f"{mB:.3f} M☉"
    
    changeDist.value = dist
    change_mA.value = mA
    change_mB.value = mB
    
    C_a = starA.mass / mass_from_radius(starA.radius, 1, k)
    C_b = starB.mass / mass_from_radius(starB.radius, 1, k)
    Ca_text.text = f"{C_a:.3f}"
    Cb_text.text = f"{C_b:.3f}"
    
    draw_potential()
    
    q_text.text = f"{q:.3f}"
    rocheA_text.text = f"{lobe_rad/Rs:.3f} R☉"
    rocheB_text.text = f"{lobe_rad2/Rs:.3f} R☉"
    
    type_text.text = "detached"
    
    ar_graph.data = []
    br_graph.data = []
    alobe_graph.data = []
    blobe_graph.data = []
    am_graph.data = []
    bm_graph.data = []
    av_graph.data = []
    bv_graph.data = []
    q_graph.data = []
    p_graph.data = []
    
    scene.autoscale = False
    scene.range = scale
    scene.autoscale = True
    
        

# sliders to adjust dist, mA, mB
def update_system():
    global x1, x2, dist, sep_dist, q, q2, lobe_rad, lobe_rad2, sum_mass, reduced_mass, momentum, C_a, C_b
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
    
    C_a = starA.mass / mass_from_radius(starA.radius, 1, k)
    C_b = starB.mass / mass_from_radius(starB.radius, 1, k)
    Ca_text.text = f"{C_a:.3f}"
    Cb_text.text = f"{C_b:.3f}"
    
    draw_potential()
    
    q_text.text = f"{q:.3f}"
    rocheA_text.text = f"{lobe_rad/Rs:.3f} R☉"
    rocheB_text.text = f"{lobe_rad2/Rs:.3f} R☉"
    
    A_overflow = starA.radius >= lobe_rad
    B_overflow = starB.radius >= lobe_rad2
    if A_overflow and B_overflow:
        type = "contact"
    elif A_overflow or B_overflow:
        type = "semi-detached"
    else:
        type = "detached"
    
    type_text.text = type
        

scene.append_to_caption("Upon loading the program, the <b>scene</b> will display a binary star system with their <b>Roche lobe equipotential</b> drawn.\n")
scene.append_to_caption("The <b>separation distance</b>, <b>starA mass</b>, and <b>starB mass sliders</b> can be used to adjust those values respectively and the equipotential\n")
scene.append_to_caption("is redrawn to reflect those changes.\n")
scene.append_to_caption("The <b>Star A Overflow</b> and <b>Star B Overflow buttons</b> enlarge the radius of the star beyond its Roche lobe to test <b>semi-detached</b> and <b>contact binaries</b>.\n")
scene.append_to_caption("The <b>Equipotential button</b> toggles continuous Roche lobe updates. Click the <b>run button</b> to begin the simulation.\n")
scene.append_to_caption("Click the <b>reset button</b> to reset the inputs to their initial values.\n")
scene.append_to_caption("All the graphs display properties/info of the stars over time like their velocities, mass, and radius.\n\n")
scene.append_to_caption("We approximated the <b>density per unit volume</b> of the stars using an exponential function \\(e^{-kt}\\) for some constant \\(k\\), which\n")
scene.append_to_caption("we will call the <i>mass density constant</i> in this project. This is meant to give the viewer an idea of the relative densities of the\n")
scene.append_to_caption("stars compared to each other.\n\n")
    
    
#######################
# DISPLAY INFO
######################
scene.append_to_caption("\nSystem Type = ")
type_text = wtext(text="detached")

scene.append_to_caption("\nMass Ratio = ")
q_text = wtext(text=f"{q:.3f}")

scene.append_to_caption("\nRobe Lobe A Radius = ")
rocheA_text = wtext(text=f"{lobe_rad2/Rs:.3f} R☉")

scene.append_to_caption("\nRobe Lobe B Radius = ")
rocheB_text = wtext(text=f"{lobe_rad2/Rs:.3f} R☉")

scene.append_to_caption("\nMass Density Const of A = ")
Ca_text = wtext(text=f"{C_a:.3f}")
scene.append_to_caption("\nMass Density Const of B = ")
Cb_text = wtext(text=f"{C_b:.3f}")
MathJax.Hub.Queue(["Typeset",MathJax.Hub])

#######################
# CHANGE RADIUS
#######################
scene.append_to_caption("\n\n")
drawP_btn = button(text="Stop Drawing Equipotential", bind=change_draw)
scene.append_to_caption("\n\n")
A_over = button(text="StarA Overflow", bind=setA_radius)
B_over = button(text="StarB Overflow", bind=setB_radius)


#######################
# DISTANCE SLIDER
#######################
scene.append_to_caption("\n\n")
scene.append_to_caption("Separation Dist. (Solar Radii): ")
scene.append_to_caption("\n")
def changeDistSlider(evt):
    global dist
    dist = evt.value
    dist_text.text = f"{evt.value/Rs:.3f} R☉"
    update_system()
changeDist = slider(bind=changeDistSlider, min=0.5*1e11, max=2*1e11, value=dist, length=300)
dist_text = wtext(text=f"{dist/Rs:.3f} R☉")


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
#scene.append_to_caption("\n\n")
#drawP_btn = button(text="Stop Drawing Equipotential", bind=change_draw)



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
Ca_text.text = f"{C_a:.3f}"
Cb_text.text = f"{C_b:.3f}"

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

drawP = True
def change_draw():
    global drawP
    drawP = not drawP
    if drawP:
        drawP_btn.text = "Stop Drawing Equipotential"
    else:
        drawP_btn.text = "Draw Equipotential"

scale = scene.range
    
t=0; dt=3600
# custom inc in rad for now
rad_inc_rate = 5e-6 * Rs
transfer_rate = 1e-7

graph_w = 480
graph_h = 360

rad_graph = graph(title='Radius Over Time', xtitle='Time (days)', ytitle='Radius (R☉)', width=graph_w, height=graph_h, align='left')
ar_graph = gcurve(color=color.red, label='StarA Radius')
br_graph = gcurve(color=color.blue, label='StarB Radius')

#roche_graph = graph(title='Roche Lobe Radius Over Time', xtitle='Time (days)', ytitle='Lobe Radius (R☉)', width=graph_w, height=graph_h, align='left')
alobe_graph = gcurve(color=color.orange, label='Roche Lobe A Radius')
blobe_graph = gcurve(color=color.purple, label='Roche Lobe B Radius')

mass_graph = graph(title='Mass Over Time', xtitle='Time (days)', ytitle='Mass (M☉)', width=graph_w, height=graph_h, align='left')
am_graph = gcurve(color=color.red, label='StarA Mass')
bm_graph = gcurve(color=color.blue, label='StarB Mass')

vel_graph = graph(title='Velocity Over Time', xtitle='Time (days)', ytitle='Velocity (m/s)', width=graph_w, height=graph_h, align='left')
av_graph = gcurve(color=color.red, label='StarA Velocity')
bv_graph = gcurve(color=color.blue, label='StarB Velocity')

ratio_graph = graph(title='Mass Ratio of StarB to StarA Over Time', xtitle='Time (days)', ytitle='Mass Ratio', width=graph_w, height=graph_h, align='left')
q_graph = gcurve(color=color.green, label='Mass Ratio of B to A')

period_graph = graph(title='Orbital Period Over Time', xtitle='Time (days)', ytitle='Period (days)', width=graph_w, height=graph_h, align='left')
p_graph = gcurve(color=color.cyan, label='Orbital Period')

type = "detached"

keep_running = True

while(True):
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
        
        type_text.text = type
        
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
        elif type == "contact" and keep_running:
            if (starA.pos-starB.pos).mag < (starA.radius+starB.radius):
                keep_running = False
                if drawP:
                    change_draw()
                starA.visible = False
                starB.visible = False
                mergedStar.radius = (starA.radius + starB.radius) * 0.95
                mergedStar.mass = starA.mass + starB.mass
                mergedStar.visible = True
            else:
                overflow1 = starA.radius - lobe_rad
                overflow2 = starB.radius - lobe_rad2
                dm1 = mass_from_radius(starA.radius, C_a, k) - mass_from_radius(lobe_rad, C_a, k)
                dm2 = mass_from_radius(starB.radius, C_b, k) - mass_from_radius(lobe_rad2, C_b, k)
                starA.mass += (dm2 - dm1) * transfer_rate * dt
                starB.mass += (dm1 - dm2) * transfer_rate * dt
                
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
                sep_dist = mag(starB.pos - starA.pos)
                v_rel = momentum / (reduced_mass * sep_dist)
                r_hat = hat(starB.pos - starA.pos)
                t_hat = vector(-r_hat.y, r_hat.x, 0)
                vA = v_rel * starB.mass / (starA.mass + starB.mass)
                vB = v_rel * starA.mass / (starA.mass + starB.mass)
                starA.velocity = -vA * t_hat
                starB.velocity =  vB * t_hat
                
                lobe_rad = sep_dist * (0.49 * q2 ** .6666667) / (0.6 * q2 ** .6666667 + log(1 + q2 ** .3333333))
                lobe_rad2 = sep_dist * (0.49 * q ** .6666667) / (0.6 * q ** .6666667 + log(1 + q ** .3333333))
            
#        if t % (3600 * 1000) == 0:
#            print("radius of star B: " + starB.radius)
#            print("total momentum: " + mag(P))
        if t % (3600 * 10) == 0 and (starA.pos-starB.pos).mag > (starA.radius+starB.radius) and drawP:
            draw_potential()
            
        if keep_running:
            q_text.text = f"{q:.3f}"
            rocheA_text.text = f"{lobe_rad/Rs:.3f} R☉"
            rocheB_text.text = f"{lobe_rad2/Rs:.3f} R☉"
            dist_text.text = f"{sep_dist/Rs:.3f} R☉"
            mA_text.text = f"{starA.mass/M0:.3f} M☉"
            mB_text.text = f"{starB.mass/M0:.3f} M☉"
            
            
            rad_graph.select()
            ar_graph.plot(t/(3600*24), starA.radius/Rs)
            br_graph.plot(t/(3600*24), starB.radius/Rs)
            
    #        roche_graph.select()
            alobe_graph.plot(t/(3600*24), lobe_rad/Rs)
            blobe_graph.plot(t/(3600*24), lobe_rad2/Rs)
            
            mass_graph.select()
            am_graph.plot(t/(3600*24), starA.mass/M0)
            bm_graph.plot(t/(3600*24), starB.mass/M0)
            
            vel_graph.select()
            av_graph.plot(t/(3600*24), mag(starA.velocity))
            bv_graph.plot(t/(3600*24), mag(starB.velocity))
            
            ratio_graph.select()
            q_graph.plot(t/(3600*24), q)
            
            period_graph.select()
            T = 2*pi*sqrt(sep_dist**3 / (G*sum_mass))
            p_graph.plot(t/(3600*24), T/(3600*24))
                 
        
        t = t+dt