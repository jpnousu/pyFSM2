
import numpy as np


def tridiag(Nvec, Nmax, a, b, c, r):
    '''
    Input
    Nvec: Vector length
    Nmax: Maximum vector length
    a: Below-diagonal matrix elements
    b: Diagonal matrix elements
    c: Above-diagonal matrix elements
    r: Matrix equation rhs
        
    Output
    x: Solution vector
    '''

    x = np.zeros(Nmax)
    g = np.zeros(Nmax)
        
    beta = b[0]
    x[0] = r[0] / beta

    for n in range(1, Nvec):
        g[n] = c[n-1] / beta
        beta = b[n] - a[n] * g[n]
        x[n] = (r[n] - a[n] * x[n-1]) / beta

    for n in range(Nvec - 2, 0, -1):
        x[n] = x[n] - g[n + 1] * x[n + 1]

    return x


def ludcmp(N, A, b):
    '''
    #
    Solve matrix equation Ax = b for x by LU decomposition
    #
    
    Args:
    N # Number of equations to solve
    A(N,N) # Matrix
    b(N) # RHS of matrix equation
    Out:
    x(N) # Solution of matrix equation

    integer :: i,ii,imax,j,k,ll,indx(N)

    real :: Acp(N,N),aamax,dum,sum,vv(N)
    '''

    Acp = A[:,:]
    x = b[:]

    vv = np.zeros(N)
    indx = np.zeros(N, dtype=int)

    for i in range(N):
        aamax = 0
        for j in range(N):
            if (abs(Acp[i,j]) > aamax):
                aamax = abs(Acp[i,j])
        vv[i] = 1/aamax

    for j in range(N):
        for i in range(j):
            sum = Acp[i,j]
            if (i > 1):
                for k in range(i):
                    sum -= Acp[i,k] * Acp[k,j]
                Acp[i,j] = sum
                    
        aamax = 0
        for i in range(j, N):
            sum = Acp[i,j]
            for k in range(j):
                sum -= Acp[i,k] * Acp[k,j]
            Acp[i,j] = sum

            dum = vv[i] * abs(sum)
            if dum >= aamax:
                imax = i
                aamax = dum
        if j != imax:
            for k in range(N):
                dum = Acp[imax, k]
                Acp[[imax, j], :] = Acp[[j, imax], :]
            vv[imax] = vv[j]

        indx[j] = imax
        if (Acp[j,j] == 0):
            Acp[j,j] = 1e-20
        if j != N-1:
            dum = 1 / Acp[j,j]
            for i in range(j+1, N):
                Acp[i,j] *= dum

    ii = 0
    for i in range(N):
        ll = indx[i]
        sum = x[ll]
        x[ll] = x[i]

        if ii != 0:
            for j in range(ii, i):
                sum -= Acp[i,j] * x[j]
        elif sum != 0:
            ii = i

        x[i] = sum

    for i in range(N-1, 0, -1):
        sum = x[i]
        for j in range(i+1, N):
            x[i] = sum / Acp[i,i]

    return x