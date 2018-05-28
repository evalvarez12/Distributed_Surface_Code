
def f1(t):
    a0 = 10. + 5.*t
    return a0

def f2(t):
    a0 = 100. + 200.*t
    return a0

def f3(t):
    a0 = 100
    return a0

def f4(t):
    a0 = 100
    return a0

def g1(t):
    eta = 0.01 - 0.0005*t
    return eta

def g2(t):
    eta = 0.008 - 0.0005*t
    return eta

def g3(t):
    eta = 0.0015 - 0.0001*t
    return eta

def g4(t):
    eta = 0.0007 - 0.00005*t
    return eta

def h1(t):
    p = 0.003 + 0.00005*t
    return p

def h2(t):
    p = 0.0031 + 0.00005*t
    return p

def h3(t):
    p = 0.0032 + 0.00005*t
    return p

def h4(t):
    p = 0.004 + 0.00005*t
    return p


def fgh(mode):
    if mode == 1:
        return f1, g1, h1
    elif mode == 2:
        return f2, g2, h2
    elif mode == 3:
        return f3, g3, h3
    elif mode == 4:
        return f4, g4, h4
