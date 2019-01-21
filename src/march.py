from __future__ import print_function
import numpy as np

# Exact solution for q
def q(t):
    return 50 - 2*np.exp(-0.196*t)

# Exact qdot (from differential equation)
def qd(t):
    return 9.8 - 0.196*q(t)

def residual(qdot, q):
    """
    Residual of the ODE
    """
    return qdot - 9.8 + 0.196*q

# BDF coefficients
alpha1 = np.array([1.0, -1.0])
alpha2 = np.array([3.0/2.0, -4.0/2.0, 1.0/2.0])
alpha3 = np.array([11.0/6.0, -18.0/6.0, 9.0/6.0, -2.0/6.0])
alpha4 = np.array([25.0/12.0, -48.0/12.0, 36.0/12.0, -16.0/12.0, 3.0/12.0])
alpha5 = np.array([137.0/60.0, -300.0/60.0, 300.0/60.0, -200.0/60.0, 75.0/60.0, -12.0/60.0])
alpha6 = np.array([147.0/60.0, -360.0/60.0, 450.0/60.0, -400.0/60.0, 225.0/60.0, -72.0/60.0, 10.0/60.0])

# Make a list
alpha  = [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6]

# Sanity check coeffs
for alp in alpha:
    print(np.sum(alp))
    
# Supported order of accuracy
orders = [1 ,2, 3, 4, 5, 6]

print("original step size")
t = 0.0; h = 1.0; e1 = []
for order in orders:
    # An O(h^p) approximation of qdot using BDF formula
    qdot = 0.0
    for i in xrange(order+1):
        qdot += alpha[order-1][i]*q(t-i*h)
    qdot = qdot/h    
    e1.append(np.abs(qdot -  qd(t)))

print("reducing step size")
t = 0.0; h = h/2.0; e2 = []
for order in orders:
    # An O(h^p) approximation of qdot using BDF formula
    qdot = 0.0
    for i in xrange(order+1):
        qdot += alpha[order-1][i]*q(t-i*h)
    qdot = qdot/h
    e2.append(np.abs(qdot -  qd(t)))
    print(residual(qdot,q(t)), e2[order-1])

# Compare theoretical and numerical accuracies
for order in orders:
    pexp = np.log(e2[order-1]/e1[order-1])/np.log(2)
    print(order, pexp)

# Compute residual
print (residual(qdot, q(t)))
