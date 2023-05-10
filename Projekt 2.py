import numpy as np



#Analytiska lösningen för positionen i x-led
def xt(v0, t):
    vx = v0*np.cos(2*np.pi/9)
    return (((-(0.001)**2)*vx-0.020-vx)*np.exp(-0.001*t)-np.cos(t)*0.020* \
            (0.001**2)-0.020*0.001*np.sin(t)+((0.001**2)+1)*(0.020+vx)) / \
                (((0.001**2)+1)*0.001)



#Analytiska lösningen för positionen i y-led
def yt(v0 , t):
    vy = v0*np.sin(2*np.pi/9)
    return ((-vy*0.001-9.82)*np.exp(-0.001*t)+(-9.82*t+vy)*0.001+9.82) / \
        (0.001**2)


#Funktion för L2-felen
def eN(xtn, xn, ytn, yn):
    return (np.sqrt((xtn - xn)**2 + (ytn -yn)**2))


#Funktion för konvergensen
def q(en1, en2, k1, k2):
    return np.log10(en1/en2) / np.log10(k1/k2)


#Funktion för omskrivning av modellen. Returnerar matris med hastigheterna och accelerationerna i x-led och y-led
def F(t,u):
    x  = u[0]
    y  = u[1]
    vx = u[2]
    vy = u[3]
    return np.array([vx, vy, -c*vx + a*np.sin(t), -c*vy-g])


#Begynnelsevillkor
g= 9.82
k = [1, 0.5, 0.25, 0.125]
theta = 2*np.pi/9
T = 4
v0 = 40
c = 0.001
a = 0.020
x0 = 0
y0 = 0
x0_t = v0*np.cos(theta)
y0_t = v0*np.sin(theta)
vx0 = x0_t
vy0 = y0_t

#Antal tidssteg
N1 = int(T/k[0])
N2 = int(T/k[1])
N3 = int(T/k[2])
N4 = int(T/k[3])



f = np.array([x0,y0,vx0,vy0])

#Tiden
t1 = np.linspace(0, 4, N1+1)
t2 = np.linspace(0, 4, N2+1)
t3 = np.linspace(0, 4, N3+1)
t4 = np.linspace(0, 4, N4+1)


u1 = np.zeros((N1+1,4),dtype = float)
u2 = np.zeros((N2+1,4),dtype = float)
u3 = np.zeros((N3+1,4),dtype = float)
u4 = np.zeros((N4+1,4),dtype = float)

u1[0,:] = f
u2[0,:] = f
u3[0,:] = f
u4[0,:] = f


#Egendefinierad RK4
def RK4(k_):
    f = np.array([0,0,40*np.cos(2*np.pi/9),40*np.sin(2*np.pi/9)])
    N = int(4/k_)
    t = np.linspace(0, 4, N+1)
    u = np.zeros((N+1, 4),dtype = float)
    u[0,:] = f
    for n in range(N):
        w1=F(t[n],u[n,:])
        w2=F(t[n]+k_/2,u[n,:]+k_/2*w1)
        w3=F(t[n]+k_/2,u[n,:]+k_/2*w2)
        w4=F(t[n+1],u[n,:]+k_*w3)
        u[n+1,:] = u[n,:] + k_/6*(w1+2*w2+2*w3+w4)
    return u

# yrk4 = RK4(k[3])[0:, 1:2]
# xrk4 = RK4(k[1])[0:, 0:1]

#print('RK4(k[0]): \n', RK4(k[0]), '\n')
#print('RK4(k[1]): \n', RK4(k[1]), '\n')
# print('xrk4: \n', xrk4, '\n')
# print('xt(40, t2): \n', xt(40, t2), '\n')
# print('yrk4: \n', yrk4, '\n\n')
# print('yt(40, t4): \n', yt(40, t4), '\n')

#RK4-lösning för x
xN1 = RK4(k[0])[0:, 0:1][-1]
xN2 = RK4(k[1])[0:, 0:1][-1]
xN3 = RK4(k[2])[0:, 0:1][-1]
xN4 = RK4(k[3])[0:, 0:1][-1]


#RK4-lösning för y
yN1 = RK4(k[0])[0:, 1:2][-1]
yN2 = RK4(k[1])[0:, 1:2][-1]
yN3 = RK4(k[2])[0:, 1:2][-1]
yN4 = RK4(k[3])[0:, 1:2][-1]


#Analytisk lösning för x
AnX1 = xt(40,t1)[-1]
AnX2 = xt(40,t2)[-1]
AnX3 = xt(40,t3)[-1]
AnX4 = xt(40,t4)[-1]


#Analytisk lösning för y
AnY1 = yt(40, t1)[-1]
AnY2 = yt(40, t2)[-1]
AnY3 = yt(40, t3)[-1]
AnY4 = yt(40, t4)[-1]


#L2-felen för dom olika tidsstegen
L2_1 = eN(AnX1, xN1, AnY1, yN1)
L2_2 = eN(AnX2, xN2, AnY2, yN2)
L2_3 = eN(AnX3, xN3, AnY3, yN3)
L2_4 = eN(AnX4, xN4, AnY4, yN4)

print(L2_1)
print(L2_2)
print(L2_3)
print(L2_4, '\n')

#Konvergensen
K_12 = q(L2_1, L2_2, k[0], k[1])
K_23 = q(L2_2, L2_3, k[1], k[2])
K_34 = q(L2_3, L2_4, k[2], k[3])
K_41 = q(L2_4, L2_1, k[3], k[0])

print(K_12)
print(K_23)
print(K_34)
print(K_41)

felnorm1=np.linalg.norm(np.transpose(RK4(k[0])[0:, 0:1]) - xt(40,t1))
felnorm2=np.linalg.norm(np.transpose(RK4(k[1])[0:, 0:1]) - xt(40,t2))
felnorm3=np.linalg.norm(np.transpose(RK4(k[2])[0:, 0:1]) - xt(40,t3))
felnorm4=np.linalg.norm(np.transpose(RK4(k[3])[0:, 0:1]) - xt(40,t4))
print("Felnormen blir {}.".format(felnorm1))
print("Felnormen blir {}.".format(felnorm2))
print("Felnormen blir {}.".format(felnorm3))
print("Felnormen blir {}.".format(felnorm4))




#slut på emils försök

















# #Uppgift 2-6


# #Relativa hastigheten av kanonkulan
# def v(t):
#     return np.array([[x_t - w(t)], [y_t]])


# #Friktionen
# def Fd(a, t):
#     return -a * np.abs(v(t)) * v(t)


# #Gravitationen beroende av höjden y
# def g(y):
#     return 3.986 * 10**14 / (6.371 * 10**6 + y)**2


# #Luftmotståndskoefficienten beroende av c0 och höjden y
# def c(c0, y):
#     return c0 * np.exp(-1 * 10**(-4) * y)


# #Accelerationen i x-led
# def x_tt(y, t):
#     return -c(y) * np.abs(v(t)) * (x_t - w(t))


# #Accelerationen i y-led
# def y_tt(y, t):
#     return -c(y) * np.abs(v(t)) * y_t - g(y)


# #Osäkert vad detta är
# def F(t, u):
#     x, y, x_t, y_t = u
#     Vx = x_t    #Inför nya obekanta
#     Vy = y_t    #Inför nya obekanta
#     return np.array([x], [y], [Vx], [Vy])


# #Vinden (uppg 2-6)
# def w(t):
#     return -20*np.exp(-((t-10)/5)**2)


# #Inför nya obekanta
# Vx_t = x_tt
# Vy_t = y_tt


# #Begynnelsevillkor
# x0 = 0
# y0 = 0
# c0 = 4.518 * 10**(-4)
# v0 = 400
# vx = v0 * np.cos(θ)
# vy = v0 * np.sin(θ)
# x_t0 = vx
# y_t0 = vy