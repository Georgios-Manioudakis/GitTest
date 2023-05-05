import numpy as np
print('din mamma')

#Uppgift 1





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
