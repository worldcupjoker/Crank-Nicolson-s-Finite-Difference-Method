import numpy as np
import matplotlib.pyplot as plt

# Define Crank-Nicolson scheme for initial boundary value problem.
def pdeCN(a0, b0, D, uLeftBound, uRightBound, uInitial, tStop = 10, dt = 0.25, n = 4):
    a0 *= 1.0
    b0 *= 1.0
    dt *= 1.0
    D *= 1.0
    tStop *= 1.0
    timeSize = int(tStop / dt) + 2
    dx = (b0 - a0) / n
    xSize = n + 2
    wSize = xSize - 3 # For unknows in the mesh
    lambd = D * dt / (2.0 * dx ** 2.0)
    wN = np.empty([wSize, 1], dtype = float)

    # Generate result matrix. The first row is time axis. The first column is x axis.
    result = np.zeros((xSize, timeSize))
    for i in range(1, xSize):
        result[i, 0] = a0 + (i - 1) * dx
        result[i, 1] = uInitial(result[i, 0]) # Initial values

    for j in range(1, timeSize):
        result[0, j] = (j - 1) * dt
        result[1, j] = uLeftBound(result[0, j]) # Left boundary
        result[-1, j] = uRightBound(result[0, j]) # Right Boundary


    # Create A matrix.
    A = np.identity(wSize)
    A[0, 0] = 2.0
    A[0, 1] = -1.0
    A[-1, -1] = 2.0
    A[-1, -2] = -1.0
    for i in range(1, len(A) - 1):
        A[i, i] = 2.0
        A[i, i - 1] = -1.0
        A[i, i + 1] = -1.0


    # Create an identity matrix.
    I = np.identity(wSize)

    # Solve the tridiagonal linear system at each time step.
    leftA = np.add(I, np.dot(lambd, A))
    coeW = np.subtract(I, np.dot(lambd, A))
    for j in range(2, timeSize):

        # Setup the linear system.
        wN[:, 0] = result[2 : -1, j - 1]
        temp1 = np.matmul(coeW, wN)
        temp2 = np.dot(lambd, np.add(bMatrix(uLeftBound, uRightBound, result[0, j - 1], wSize), bMatrix(uLeftBound, uRightBound, result[0, j], wSize)))
        rightB = np.subtract(temp1, temp2)

        # Solve for w.
        result[2 : -1, j] = tridiagonal(leftA, rightB)[:, 0]

    return result


# Essential function for pdeCN
def bMatrix(uLeftBound, uRightBound, time, size):
    b = np.zeros((size, 1))
    b[0, 0] = -1.0 * uLeftBound(time)
    b[-1, 0] = -1.0 * uRightBound(time)
    return b


# Solve a linear system that contains a tridiagonal matrix.
def tridiagonal(A, b):
    n = len(A)

    # Create solution vector.
    z = np.empty([n, 1], dtype = float)
    x = np.empty([n, 1], dtype = float)

    # LU decomposition
    L = np.zeros((n, n))
    U = np.identity(n)
    L[0, 0] = A[0, 0]
    L[1, 0] = A[1, 0]
    U[0, 1] = A[0, 1] / L[0, 0]

    for i in range(1, n - 1):
        L[i, i] = A[i, i] - L[i, i - 1] * U[i - 1, i]
        L[i + 1, i] = A[i + 1, i]
        U[i, i + 1] = A[i, i + 1] / L[i, i]

    L[-1, -1] = A[-1, -1] - L[-1, -2] * U[-2, -1]

    # Forward substitution
    z[0, 0] = b[0, 0] / L[0, 0]
    for i in range(1, n):
        for j in range(i - 1, i):
            z[i, 0] = (b[i, 0] - L[i, j] * z[i - 1, 0]) / L[i, i]



    # Backward substitution
    x[n - 1, 0] = z[n - 1, 0]
    for i in range(n - 2, -1, -1):
        for j in range(i + 1, i, -1):
            x[i, 0] = z[i, 0] - U[i, j] * x[i + 1, 0]


    return x


####
def uT0(x):
    return 0.0


def uLeft(t):
    if (t <= 2.0):
        return 1.0 / 2 * t ** 2.0
    else:
        return 2.0 * np.exp(-(t - 2.0))



def uRight(t):
    if (t <= 2.0):
        return 1.0 / 2 * t ** 2.0
    else:
        return 2.0 * np.exp(-(t - 2.0))



####
alpha = 509.76

# (a)
print("(a)")
print("The solution of h(x, t) at t=2 is plotted below.")
r4a = pdeCN(0, 1100, alpha, uLeft, uRight, uT0, 2, 0.005, 500)
plt.plot(r4a[1 :, 0], r4a[1 :, -1], label = "t = 2 days")
plt.xlabel("x (m)")
plt.ylabel("Approximate solution h(x) (m)")
plt.legend()
plt.grid(True)
plt.show()
print("")

# (b)
print("(b)")
print("The solution of h(x, t) at t=5, 10, 20 are plotted below.")
r4b5 = pdeCN(0, 1100, alpha, uLeft, uRight, uT0, 5, 0.005, 500)
r4b10 = pdeCN(0, 1100, alpha, uLeft, uRight, uT0, 10, 0.005, 500)
r4b20 = pdeCN(0, 1100, alpha, uLeft, uRight, uT0, 15, 0.005, 500)
plt.plot(r4b5[1 :, 0], r4b5[1 :, -1], label = "t = 5 days")
plt.plot(r4b10[1 :, 0], r4b10[1 :, -1], label = "t = 10 days")
plt.plot(r4b20[1 :, 0], r4b20[1 :, -1], label = "t = 20 days")
plt.xlabel("x (m)")
plt.ylabel("Approximate solution h(x) (m)")
plt.legend()
plt.grid(True)
plt.show()
print("")
