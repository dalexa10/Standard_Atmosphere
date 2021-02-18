from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np

# Definition of the plots' fonts
plt.rc('font',family='serif')
plt.rc('font', size = 11)
plt.rc('axes', labelsize = 11)

# Definition of plotting function
def t_R2F(T_R):
    """
    Transform a given temperature from R to F degrees
    """
    T_F = T_R - 459.67
    
    return T_F

def temp(h):
    """
    Computes the temperature at a certain altitude

    Paramters
    -----------
    h: float, altitude [ft]

    Returns
    ----------
    T: float, temperature at certain altitude [F]

    """
    if h < 36000:
        To = 518.69
        a = -3.57*(10**-3)
        ho = 0
        T = To + a*(h-ho)
        T_ratio = T/To
    else:
        T = 389.97
        T_ratio = None
    
    return T, T_ratio

def goR_ratio():
    """
    Computes the value of the exponential part of both pressure and density

    """
    go = 32.174
    R = 1716.5
    k = go/R

    return k

def press(h):
    """
    Computes the atmospheric pressure at a certain altitude

    Paramters
    -----------
    h: float, altitude [ft]

    Returns
    ----------
    P: float, pressure at certain altitude [psf]

    """
    Po = 2116.22
    a = -3.57*10**(-3)
    if h < 36000:
        T, T_ratio = temp(h)
        P = Po * T_ratio**(-goR_ratio()/a)
    else:
        T36 = 389.97
        P36 = Po * (T36/518.69)**(-goR_ratio()/a)
        P = P36 * np.exp((-goR_ratio()/T36)*(h-36000))

    return P


def dens(h):
    """
    Computes the air density at a certain altitude

    Paramters
    -----------
    h: float, altitude [ft]

    Returns
    ----------
    rho: float, density at certain altitude [slugs/ft3]

    """
    rho_o = 2.3769*10**-3
    a = -3.57*10**(-3)
    if h < 36000:
        T, T_ratio = temp(h)
        rho = rho_o * T_ratio**(-(goR_ratio()/a)-1)
    else:
        T36 = 389.97
        rho36 = rho_o * (T36/518.69)**(-(goR_ratio()/a)-1)
        rho = rho36 * np.exp((-goR_ratio()/T36)*(h-36000))

    return rho

def air_prop(h):
    """
    Computes the thermodynamic air properties (T, P and rho) for a
    desired altitude

    Paramters
    -----------
    h: float, altitude [ft]

    Returns
    ----------
    T: float, temperature [R]
    P: float, pressure [psf]
    rho: float, density at certain altitude [slugs/ft3]
    """
    T = temp(h)
    P = press(h)
    rho = dens(h)

    return T, P, rho


# ========================================
# ----------- Test Zone ------------------
# ========================================

#Evalute air properties at h = 4000 [ft]

result = air_prop(74123)
T = t_R2F(result[0][0])
P = result[1]
rho = result[2]






