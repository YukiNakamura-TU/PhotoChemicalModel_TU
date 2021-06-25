from scipy.special import gammainc # lower incomplete gamma function / gamma
from scipy.special import gamma
import numpy as np #mathematical library

BJ = 4.35e5 # nT
RJ = 1 

M = 1

D = 2*240*M+1
#D = 10

b=[[0,0] for i in range(D)]
bold=[[0,0] for i in range(D)]
theta=[[0,0] for i in range(D)]
db_dL=[[0,0] for i in range(D)]
f = [[0,0,0] for i in range(D)]
print("start calculation")

for ix in range(D):
    L = 1.0+float(ix)/float(2*M)

    alpha = -2.0/5.0
    x = (L/14.501)**2.5

    uicg_p1 = gamma(alpha+1)*(1.0-gammainc(alpha+1,x))
    uicg = (uicg_p1 - x**alpha * np.exp(-x))/alpha

    Fe = 2.841e4 + 9.199e3 * uicg + 5.4e4/(2.71-2.0) * (1.0/L)**(2.71-2.0)

    bold[ix][0] = L
    bold[ix][1] = 1.0/L

    b[ix][0] = L
    if L <  5.0:
        b[ix][1] = 1.0/L
    if L >= 5.0:
        b[ix][1] = Fe/BJ

    db_dL[ix][0] = L
    if L <  5.0:
        db_dL[ix][1] = -1.0/L/L
    if L >= 5.0:
        db_dL[ix][1] = ( - 9.199e3 * 2.5 * 14.501* L**(-2.0) * (np.exp(-(L/14.501)**2.5)) - 5.4e4 * L**(-1.71) ) / BJ

    theta[ix][0] = L
    if L <  5.0:
        theta[ix][1] = 90-np.degrees(np.arcsin(np.sqrt(1/L)))
    if L >= 5.0:
        theta[ix][1] = 90-np.degrees(np.arcsin(np.sqrt(Fe/BJ/RJ/RJ)))

    f[ix][0] = L
    f[ix][1] = -b[ix][1] /L/L * np.sqrt(1.0-b[ix][1]) * db_dL[ix][1]
    f[ix][2] = np.sqrt(1.0-1/L) / L**5.0

    if ix%10000 == 0:
        print(float(ix)/float(2*240*M) * 100)

    #print(L, theta[ix][1])
    
#print((b[6][1]-b[4][1])*float(M), db_dL[5][1])

out = np.array(b)
np.savetxt("b_function/b.dat", out)

out = np.array(bold)
np.savetxt("b_function/bold.dat", out)

out = np.array(theta)
np.savetxt("b_function/theta.dat", out)

out = np.array(db_dL)
np.savetxt("b_function/db_dL.dat", out)

out = np.array(f)
np.savetxt("b_function/f.dat", out)