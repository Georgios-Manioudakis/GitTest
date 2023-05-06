import numpy as np


#Uppgift 1


#Emils försök på uppgift 1

def xt(v0, t):
    vx = v0*np.cos(2*np.pi/9)
    return (((-0.001**2)*vx-0.020-vx)*np.exp(-0.001*t)-np.cos(t)*0.020*(0.001**2)-0.020*0.001*np.sin(t)+((0.001**2)+1)*(0.020+vx))/(((0.001**2)+1)*0.001)

def yt(v0 , t):
    vy = v0*np.sin(2*np.pi/9)
    return ((-vy*0.001-9.82)*np.exp(-0.001*t)+(-9.82*t+vy)*0.001+9.82)/(0.001**2)

def eN(xtn, xn, ytn, yn):
    fel = (xtn - xn, ytn -yn)
    return np.abs(fel)

def q(en1, en2, k1, k2):
    return np.log10(en1/en2)/np.log10(k1/k2)

def F(t,u):
    x  = u[0]
    y  = u[1]
    vx = u[2]
    vy = u[3]
    return np.array([vx, vy, x, y])

g= 9.82
k = [1, 0.5, 0.25, 0.125]
theta = 2*np.pi/9
T = 4
N1 = T/k[0]
N2 = T/k[1]
N3 = T/k[2]
N4 = T/k[3]
v0 = 40
c = 0.001
a = 0.020

#beg data
x0 = 0
y0 = 0
x0_t = v0*np.cos(theta)
y0_t = v0*np.sin(theta)
vx0 = x0_t
vy0 = y0_t

f = np.array([x0,y0,vx0,vy0])

t1 = np.linspace(0, 4, int(N1+1))
t2 = np.linspace(0, 4, int(N2+1))
t3 = np.linspace(0, 4, int(N3+1))
t4 = np.linspace(0, 4, int(N4+1))

u1 = np.zeros((int(N1+1),4),dtype = float)
u2 = np.zeros((int(N2+1),4),dtype = float)
u3 = np.zeros((int(N3+1),4),dtype = float)
u4 = np.zeros((int(N4+1),4),dtype = float)

u1[0,:] = f
u2[0,:] = f
u3[0,:] = f
u4[0,:] = f

def RK4(k_, i):
    f = np.array([0,0,40*np.cos(2*np.pi/9),40*np.sin(2*np.pi/9)])
    N = 4/k_[i]
    k_ = k_[i]
    N = int(N)
    t = np.linspace(0, 4, int(N+1))
    u = np.zeros((int(N+1),4),dtype = float)
    u[0,:] = f
    for n in range(N):
        w1=F(t[n],u[n,:])
        w2=F(t[n]+k_/2,u[n,:]+k_/2*w1)
        w3=F(t[n]+k_/2,u[n,:]+k_/2*w2)
        w4=F(t[n+1],u[n,:]+k_*w3)
        u[n+1,:] = u[n,:] + k_/6*(w1+2*w2+2*w3+w4)
    return u

yrk4 = RK4(k, 1)[0:, 1:2]
xrk4 = RK4(k, 1)[0:, 0:1]


print(RK4(k, 0))
print(xrk4)
print(xt(40, t2))
print(yrk4)
print(yt(40, t2))

#slut på emils försök

def v(t):
    return np.array([[x_t - w(t)], [y_t]])


def Fd(a, t):
    return -a * np.abs(v(t)) * v(t)


def g(y):
    return 3.986 * 10**14 / (6.371 * 10**6 + y)**2


def c(c0, y):
    return c0 * np.exp(-1 * 10**(-4) * y)


def x_tt(y, t):
    return -c(y) * np.abs(v(t)) * (x_t - w(t))


def y_tt(y, t):
    return -c(y) * np.abs(v(t)) * y_t - g(y)





x0 = 0
y0 = 0
