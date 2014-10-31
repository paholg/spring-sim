#!/usr/bin/python3
import sys
from pylab import *
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib import animation

rc('font', family='serif', serif ='Computer Modern')
rc('text', usetex=True)

params = ['show-anim', 'save-anim', 'show-plots', 'save-plots']
for arg in sys.argv[1:]:
    if arg not in params:
        print('The argument %s is not an allowed parameter.' %arg)
        print('Please run with any of:')
        for p in params:
            print('\t',p)
        print('But not both show-anim and save-anim.')
        exit(1)
    if 'show-anim' in sys.argv and 'save-anim' in sys.argv:
        print('Please don\'t run with both save-anim and show-anim.')
        exit(1)

# system properties
m1 = 55000 # kg
m2 = 41000 # kg

k1 = 120e3 # N/m
k2 = 4100e3 # N/m
kC = 3e9 # N/m

c2 = 0.05*2*sqrt(m2*k2)

g = 0 #-9.81 # m/s^2

dt = .001

# initial conditions
vL = -0.2 # m/s
x0 = array([2, 4.5*.0254, 0.]) # m
v = array([vL, vL, 0.])
a = array([0., 0., 0.])

# animation parameters
dim = .1
spring_width = .05
wind1 = 20
wind2 = 5
ground = -.5
Time = 8 # s


x = zeros_like(x0)
x[:] = x0

x1 = array([x[1]])
x2 = array([x[2]])

I = zeros(1)

def impact(x):
    a = zeros(3)
    if x[1] <= x[2]:
        a[1] = -kC*(x[1] - x[2])/m1
        a[2] = -kC*(x[2] - x[1])/m2
    return a

def accel(x, v):
    a = zeros(3)
    a[1] = -k1*(x[1] - x[0] - (x0[1] - x0[0]))/m1 + g
    a[2] = -k2*x[2]/m2 - c2*v[2]/m2 + g
    a += impact(x)
    return a

def rk4(x, v, dt):
        x1 = x
        v1 = v
        a1 = accel(x1, v1)

        x2 = x + 0.5*v1*dt
        v2 = v + 0.5*a1*dt
        a2 = accel(x2, v2)

        x3 = x + 0.5*v2*dt
        v3 = v + 0.5*a2*dt
        a3 = accel(x3, v3)

        x4 = x + v3*dt
        v4 = v + a3*dt
        a4 = accel(x4, v4)

        xf = x + (dt/6)*(v1 + 2*v2 + 2*v3 + v4)
        vf = v + (dt/6)*(a1 + 2*a2 + 2*a3 + a4)
        return xf, vf

def draw_spring(ymin, ymax, windings, width, end = .05, curl = .025):
    dt = (ymax-ymin)/windings/100
    t = arange(ymin, ymax, dt)
    x = zeros_like(t)
    middle = (t > ymin + end) & (t < ymax - end - curl)
    if len(t[middle]) == 0:
        return x, t
    coil_len = t[middle][-1] - t[middle][0]
    xcoil = width*sin(2*pi*(t[middle] - t[middle][0])/coil_len*windings)
    ycoil = t[middle] - curl*cos(2*pi*(t[middle] - t[middle][0])/coil_len*windings) + curl
    x[middle] = xcoil
    y = zeros_like(t)
    y[:] = t
    y[middle] = ycoil
    return x, y

def init():
    ax.add_patch(patch0)
    ax.add_patch(patch1)
    ax.add_patch(patch2)

    ax.add_line(line1)
    ax.add_line(line2)
    return patch0, patch1, patch2, line1, line2

skip = 30
time = 0
def animate(i):
    global x, v, x1, x2, I, time
    for i in range(skip):
        x, v = rk4(x, v, dt)
        x1 = append(x1, x[1])
        x2 = append(x2, x[2])
        I = append(I, impact(x)[1]*m1)
        time += dt
    if 'show-anim' in sys.argv or 'save-anim' in sys.argv:
        patch0.xy = (-dim, x[0])
        patch1.xy = (-dim/2, x[1])
        patch2.xy = (-dim/2, x[2]-dim)

        line1.set_data(*draw_spring(x[1] + dim, x[0], wind1, spring_width))
        line2.set_data(*draw_spring(ground, x[2] - dim, wind2, spring_width))
        ax.set_title('$t = %.2f$ s' %time)
    #i += 1
    print('t:', time)

    return patch0, patch1, patch2, line1, line2

# set up animation stuff
patch0 = Rectangle((-dim, x[0]), 2*dim, .1*dim, color='black')
patch1 = Rectangle((-dim/2, x[1]), dim, dim, color='blue')
patch2 = Rectangle((-dim/2, x[2]-dim), dim, dim, color='red')

dy = .001
y1 = arange(x[1] + dim, x[0] + dy/2, dy)
y2 = arange(ground, x[2] - dim + dy/2, dy)
len1 = x[0] - x[1] - dim
len2 = x[2] - dim - ground
line1 = Line2D(*draw_spring(x[1]+dim, x[0], wind1, spring_width), color='k')
line2 = Line2D(*draw_spring(ground, x[2]-dim, wind2, spring_width), color='k')


animfig = figure(figsize=(4,6))
ax = animfig.add_subplot(111)
ax.set_aspect('equal')
ax.set_ylim(-.5, 1.0)
ax.set_xlim(-3*dim, 3*dim)
ax.set_xticks([])
ax.set_ylabel('$x$ (m)')

if 'show-anim' in sys.argv or 'save-anim' in sys.argv:
    anim = animation.FuncAnimation(animfig, animate, init_func=init,
                                   frames=int(round(Time/(skip*dt))), interval = dt)
    if 'show-anim' in sys.argv:
        show()
    if 'save-anim' in sys.argv:
        anim.save('spring-anim.gif', writer='imagemagick', fps=((1/(skip*dt))))
else:
    i = 0
    while time < 20:
        animate(i)

# Plots in metric
# posfig = figure()
# t = arange(0, len(x1)*dt, dt)
# plot(t, x1, 'b', label='$m_1$')
# plot(t, x2, 'r', label='$m_2$')
# xlabel('$t$ (s)')
# ylabel('$x$ (m)')
# xlim(0, int(t[-1]))
# title('Positions of the masses ($x_1$ initial is %.2f m)' %x0[1])
# legend()
# if 'save' in sys.argv:
#     savefig('xs-nog.png')

# forcefig = figure()
# plot(t, I, label='$F_\\mathrm{impact}$')
# xlabel('$t$ (s)')
# ylabel('$F$ (N)')
# xlim(0, int(t[-1]))
# title('Impact Force')
# if 'save' in sys.argv:
#     savefig('F-nog.png')


# Plots in inches and pounds
posfig = figure()
t = arange(0, len(x1)*dt, dt)
plot(t, -x1/.0254, 'b', label='$m_1$')
plot(t, -x2/.0254, 'r', label='$m_2$')
xlabel('$t$ (s)')
ylabel('$x$ (in)')
xlim(0, int(t[-1]))
ylim(-2, 5)
title('Positions of the masses ($x_1$ initial is %.2f in)' %(x0[1]/.0254))
legend(loc='best')
if 'save-plots' in sys.argv:
    savefig('positions.png')

forcefig = figure()
plot(t, I*.225, label='$F_\\mathrm{impact}$')
xlabel('$t$ (s)')
ylabel('$F$ (lbs)')
xlim(0, int(t[-1]))
title('Impact Force')
if 'save-plots' in sys.argv:
    savefig('impact-force.png')

if 'show-plots' in sys.argv:
    show()
