# Code for CHEME 7770 PS6 Problem 2 and Problem 3
# Justin Paek

# Problem 2
visc = 0.001
del_q = -1.602 * 10^-19
r = 2 * 10^-9
E = 100 * 100

t = -6*pi*visc*r/(E*del_q) /3600    #hours


#Problem 3a

#initialize parameters
N1t = 10^5                              # num/cell
N2t = 10^5                              # num/cell
Amax = 10^-6 /10000                     # m^2
kb = 1.38065*(10^-23)                   # m^2*kg/(s^2*K)
T = 37 + 273                            # Kelvin
A1 = (2*10^-6) /10000                   # m^2
A2 = (2*10^-6) /10000                   # m^2
L = (2*10^-6) /100                      # m
tau = (10^-6) /100                      # m
spring_const = 0.1 *(10^-7)             # N/m
KL = 10^-8 /10000                       # m^2

# radical = ((L + (kb*T/(spring_const*tau)))^2 + 4*kb*T/spring_const)^0.5
# S_hat = 0.5*(L + (kb*T/(spring_const*tau)) + radical)
S_hat = 0.3*L + L     # value found in Figure 7 where non-dimensional separation distance
                      # crosses from Region III to Region II

kappa = KL*exp(-0.5*spring_const*(S_hat - L)^2/(kb*T))
z = S_hat*kb*T/(exp(-S_hat/tau)*Amax)     # collection of constants to simplify things
z2 = (A1*A2)/(Amax*kappa)                 # another collection of constants

gamma_star = 0.5*(z*(N1t+N2t+z2) + sqrt(z^2*(N1t+N2t+z2)^2 - 4*N1t*N2t*z^2))