""" Uppgift 2 """



import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt




#Relativa hastigheten av kanonkulan
def v(t, x_t, y_t):
    return np.array([[x_t - w(t)], [y_t]])



#Gravitationen beroende av höjden y
def g(y):
    return 3.986 * 10**14 / (6.371 * 10**6 + y)**2


#Luftmotståndskoefficienten beroende av c0 och höjden y
def c(c0, y):
    return c0 * np.exp(-1 * 10**(-4) * y)


#Vinden (uppg 2-6)
def w(t):
    return -20*np.exp(-((t-10)/5)**2)



#Funktion för omskrivning av variabler
def F(t, u):
    x, y, x_t, y_t = u
    Vx = x_t    #Inför nya obekanta
    Vy = y_t    #Inför nya obekanta
    Vx_x = -c(c0, y) * np.linalg.norm(v(t, Vx, Vy)) * (Vx - w(t))
    Vy_y = -c(c0, y) * np.linalg.norm(v(t, Vx, Vy)) * Vy - g(y)
    return np.array([Vx, Vy, Vx_x, Vy_y])







#Begynnelsevillkor
x0 = 0
y0 = 0
c0 = 4.518 * 10**(-4)
v0 = 400
θ = np.pi/4
vx = v0 * np.cos(θ)
vy = v0 * np.sin(θ)
x_t0 = vx
y_t0 = vy





def RK4_2(k):
    T = 35
    f = np.array([0, 0, vx, vy])
    N = int(T/k)
    t = np.linspace(0, T, N+1)
    u = np.zeros((N+1, 4),dtype = float)
    u[0] = f
    for n in range(N):
        w1 = F(t[n], u[n])
        w2 = F(t[n] + k/2, u[n] + k/2 * w1)
        w3 = F(t[n] + k/2, u[n] + k/2 * w2)
        w4 = F(t[n+1], u[n] + k * w3)
        u[n+1] = u[n] + k/6 * (w1 + 2 * w2 + 2 * w3 + w4)
        if u[n+1][1] < 0.05:
            #print(n)
            # print(u[n+1][1])
            u = np.delete(u, np.s_[n+2::], 0)
            #print(u.shape)
            #print(round((n+2)*k))
            break
    return u


k1 = 0.1
k2 = 0.01
k3 = 0.001
k4 = 0.0001

print('k = 0.1: \n', RK4_2(k1), '\n\n') 
print('k = 0.01: \n', RK4_2(k2), '\n\n') 
print('k = 0.001: \n', RK4_2(k3), '\n\n') 
print('k = 0.0001: \n', RK4_2(k4), '\n\n') 

print('skillnad mellan 0.1 och 0.01: ', abs(RK4_2(k1)[-1][0]-RK4_2(k2)[-1][0]))






""" Uppgift 3 """






k = k3
th = np.linspace(0,90,91)
th2 = np.linspace(14,15,100)
th3 = np.linspace(58,59,100)

def RK4_3(k, theta):
    T = 10000
    vx = v0 * np.cos(theta*np.pi/180)
    vy = v0 * np.sin(theta*np.pi/180)
    f = np.array([0, 0, vx, vy])
    N = int(T/k)
    t = np.linspace(0, T, N+1)
    u = np.zeros((N+1, 4),dtype = float)
    u[0] = f
    for n in range(N):
        w1 = F(t[n], u[n])
        w2 = F(t[n] + k/2, u[n] + k/2 * w1)
        w3 = F(t[n] + k/2, u[n] + k/2 * w2)
        w4 = F(t[n+1], u[n] + k * w3)
        u[n+1] = u[n] + k/6 * (w1 + 2 * w2 + 2 * w3 + w4)
        if u[n+1][1] < 0.05:
            # print(u[n+1][1])
            u = np.delete(u, np.s_[n+2::], 0)
            print('matrisens form:', u.shape)
            print('Antalet tidssteg', n+2)
            print('Tiden:', (n+2)*k)
            break
    return u

#print(RK4_3(k3, 45))

#for i in range(th.size):
#    if np.abs(RK4_3(k3, th[i])[0]-2700) < 70:
#        print(th[i], RK4_3(k3, th[i]))

u3 = RK4_3(k3, 58.04441)
#v58 = RK4_3(k3, 58)
#v59 = RK4_3(k3, 59)
#v58_5 = RK4_3(k3, 58.5)
#v58_1 = RK4_3(k3, 58.1)
#v58_05 = RK4_3(k3, 58.05)
v58_04441 = u3[-1]

#print('Vinkel: 58 ger', v58, '\nAvstånd från x=0:', v58[0], ' \n|x^N - 2700|=', abs(v58[0]-2700), '\n')
#print('Vinkel: 59 ger', v59, '\nAvstånd från x=0:', v59[0], ' \n|x^N - 2700|=', abs(v59[0]-2700), '\n')
#print('Vinkel: 58.5 ger', v58_5, '\nAvstånd från x=0:', v58_5[0], ' \n|x^N - 2700|=', abs(v58_5[0]-2700), '\n')
#print('Vinkel: 58.1 ger', v58_1, '\nAvstånd från x=0:', v58_1[0], ' \n|x^N - 2700|=', abs(v58_1[0]-2700), '\n')
#print('Vinkel: 58.05 ger', v58_05, '\nAvstånd från x=0:', v58_05[0], ' \n|x^N - 2700|=', abs(v58_05[0]-2700), '\n')
print('Vinkel: 58.04441 ger', v58_04441, '\nAvstånd från x=0:', v58_04441[0], ' \n|x^N - 2700|=', abs(v58_04441[0]-2700), '\n')
print('Tiden för nedslaget:', 38509*k3)


#for i in range(th2.size):
#    if np.abs(RK4_3(k3, th2[i])[0]-2700) < 0.05:
#        print(th2[i], RK4_3(k3, th2[i]))


#for i in range(th3.size):
#    if np.abs(RK4_3(k3, th3[i])[0]-2700) < 0.05:
#        print(th3[i], RK4_3(k3, th3[i]))




def I(u):
    # Vx = v0 * np.cos(theta*np.pi/180)
    # Vy = v0 * np.sin(theta*np.pi/180)
    Vx = u[:,2]
    Vy = u[:,3]
    return np.sqrt(Vx**2 + Vy**2)


tN = 38.509
k = 0.001
N3 = int(tN/k)
t = np.arange(0, 10000, k)[:u3.shape[0]]
theta = 58.04441
S1 = simpson(I(u3), t)
print(S1)


u4 = RK4_3(2*k3, 58.04441)
t4 = np.arange(0, 10000, 2*k)[:u4.shape[0]]
S2 = simpson(I(u4), t4)
print(S2)



Fel = (S1 - S2)/15
print(Fel)

fig, ax =plt.subplots()
ax.set(xlabel= 'x[m]', ylabel='y[m]', title='Kulans bana')
plt.plot(u3[:,0],u3[:,1])
plt.show()
l1 = u3[-1][0]

c0 = 0
u3_2 = RK4_3(k3, 58.04441)
l2 = u3_2[-1][0]

fig, ax =plt.subplots()
ax.set(xlabel= 'x[m]', ylabel='y[m]', title='Kulans bana utan luftmotstånd')
plt.plot(u3_2[:,0],u3_2[:,1])
plt.show()

print(u3_2[-1])
print('Utan luftmotstånd flyger kulan', l2-l1, 'meter längre.')
t_ = np.arange(0, 10000, k)[:u3_2.shape[0]]
S3 = simpson(I(u3_2), t_)

print('Utan luftmotstånd är kulans banas längd', S3-S1, 'meter längre.')
