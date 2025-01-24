import pandas as pd
import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve

# Define basis functions and their derivatives
def phi1(x):
    return 1

def phi2(x):
    return 3 * x - 1

def phi3(x):
    return x ** 2

def phi4(x):
    return x ** 3

def phi1_prime(x):
    return 0

def phi2_prime(x):
    return 3

def phi3_prime(x):
    return 2 * x

def phi4_prime(x):
    return 3 * x ** 2

# Define the integrands for the matrix A
def integrand_A(i, j, x):
    phi_i = [phi1, phi2, phi3, phi4][i]
    phi_j = [phi1, phi2, phi3, phi4][j]
    phi_i_prime = [phi1_prime, phi2_prime, phi3_prime, phi4_prime][i]
    phi_j_prime = [phi1_prime, phi2_prime, phi3_prime, phi4_prime][j]
    
    return (phi_i_prime(x) * phi_j_prime(x) + phi_i(x) * phi_j(x))

# Calculate matrix A
def compute_A():
    A = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            A[i, j], _ = quad(lambda x: integrand_A(i, j, x), 0, 1)
    return A

# Define the integrands for the vector b
def integrand_b(j, x):
    phi_j = [phi1, phi2, phi3, phi4][j]
    return 2 * phi_j(x)

def compute_b():
    b = np.zeros(4)
    for j in range(4):
        b[j], _ = quad(lambda x: integrand_b(j, x), 0, 1)
    
    # Boundary condition adjustments
    phi = [phi1, phi2, phi3, phi4]
    b  - 5 * phi   # Adjust with boundary conditions
    return b

# Compute matrix A and vector b
A = compute_A()
b = compute_b()

# Solve the linear system A * u = b
u = solve(A, b)

# Print the coefficients
print("Coefficients of the basis functions:", u)
