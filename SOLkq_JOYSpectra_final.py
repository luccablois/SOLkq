import numpy as np
import scipy as sp
from scipy import special as spe
import sympy.physics.wigner as wgn
import os

'''
Note: In scipy the theta and phi angles are swapped.
that is: theta: (0-2pi) and phi: (0-pi)
'''
#Begin#

#Let's extract the relevant data from the JOYSpectra output and create a temporary xyz file (it will later be deleted)
f = open('jobname.log',mode='r',encoding='utf-8')
content = f.readlines()
f.close()

start = [i for i in range(len(content)) if content[i].find('Mass centered Geometry (in Angstrom):') != -1][0]
stop = [i for i in range(len(content)) if content[i].find(' INERTIA TENSOR (in amu A^2):') != -1][0]

ln = content[start:stop]
ln = [i.split() for i in ln]

for i in ln:
    try:
        if int(i[0]) == 1:
            ln = i[1]
    except:
        pass

start2 = [i for i in range(len(content)) if content[i].find('1st Coordination sphere geometry (Cartesian coordinates in Angstrom):') != -1][0]
stop2 = [i for i in range(len(content)) if content[i].find('1st Coordination sphere geometry (Spherical coordinates): ') != -1][0]

geometry = content[start2:stop2]
geometry = [i.split() for i in geometry]
prov = []

for i in geometry:
    try:
        if isinstance(int(i[0]),int) == True:
            prov.append(i)
    except:
        pass

start3 = [i for i in range(len(content)) if content[i].find('1st Coordination sphere G factors:') != -1][0]
end3 = [i for i in range(len(content)) if content[i].find('1st Coordination sphere effective Alphas: ') != -1][0]

charges = content[start3:end3]
charges = [i.split() for i in charges]
g = []

for i in charges:
    try:
        if isinstance(int(i[0]), int) == True:
            g.append(float(i[4]))
    except:
        pass
g = np.array(g)

f = open('temp.xyz',mode='w')
f.write(f'{1+len(prov)}\n')
f.write('Temporary file for Bkq computation\n')
f.write(f'{ln}\t')
f.write('0\t')
f.write('0\t')
f.write('0\n')
for i in prov:
    f.write(f'{i[1]}\t')
    f.write(f'{i[2]}\t')
    f.write(f'{i[3]}\t')
    f.write(f'{i[4]}\n')
f.close()

#Centering the Geometry around the Ln3+ ion
#Input file: .xyz file with coordinates of Eu and the 1st coordination sphere. Distances in Angstroms, space separated

Input = np.loadtxt('temp.xyz',skiprows=2,usecols=(1,2,3))

gg = [Input[i]-Input[0] for i in range(1,len(Input))]

G = np.array(gg)


#G is the coordinate of the 1st coordination around the Eu3+ atoms (Cartesian)

#Getting the atomic numbers

ipt = open('temp.xyz',mode='r')
lines = ipt.readlines()
ipt.close()
os.remove('temp.xyz')
atheadder = [lines[i+2][0:2].strip() for i in range(len(lines)-2)]

def ZLn(x):
    lanthanides = ['Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb']
    atw = [i for i in range(58,71)]
    for i in range(len(lanthanides)):
        if x == lanthanides[i]:
            return atw[i]

ZLn = ZLn(atheadder[0])
def ZLig(x):
    ligands = ['N','O','F','P','S','Cl','Se','Br','I']
    atw = [7,8,9,15,16,17,34,35,53]
    for i in range(len(ligands)):
        if x == ligands[i]:
            return atw[i]

ZL = [ZLig(atheadder[i+1]) for i in range(len(gg))]

#JOYSpectra g-factors around the atoms. Best-fit, not necessarily real charges.

#Coordinate transformation to spherical system

def r(x):
    return np.sqrt(x[0]**2+x[1]**2)

def R(x):
    return np.sqrt(x[0]**2+x[1]**2+x[2]**2)

def ph(x):
    cos = x[2]/R(x)
    return np.arccos(cos)

def th(x):
    if x[0] == 0 and x[1]>0:
        tan = np.inf
    if x[0] == 0 and x[1]<0:
        tan = -np.inf
    else:
        tan = x[1]/x[0]
    cos = x[0]/r(x)
    sin = x[1]/r(x)
    if tan>0 and sin>0 and cos>0:
        return np.arctan(tan)
    if tan<0 and sin>0 and cos<0:
        return np.arctan(tan)+np.pi
    if tan>0 and sin<0 and cos<0:
        return np.arctan(tan)+np.pi
    if tan<0 and sin<0 and cos>0:
        return np.arctan(tan)
    if tan == 0 and sin == 0 and cos == 1:
        return 0
    if tan == np.inf and sin == 1 and cos == 0:
        return np.pi/2
    if tan == 0 and sin == 0 and cos == -1:
        return np.pi
    if tan == - np.inf and sin == -1 and cos == 0:
        return 3*np.pi/2
    else:
        None

def pol(X):
    return np.array([R(X),th(X),ph(X)])

P = [pol(G[i]) for i in range(len(G))]


#P are the coordinates of the 1st coordination sphere Eu (Spherical)
#Ln-O overlap integrals (ref: Carneiro Neto A.N., Moura Jr R.T. Chem. Phys. Lett. 10.1016/j.cplett.2020.137884)

def rho(Ln,L,X):
    if L == 7:
        c = [[3.022,-4.291,0.582],
            [-15.964,16.732,-4.877],
            [-16.359,17.061,-4.963],
            [0.587,-1.504,0.040],
            [1.159,-1.995,0.151],
            [3.102,-3.328,0.331],
            [1.773,-3.320,0.279],
            [2.717,-3.263,0.413],
            [1.275,-1.812,0.097],
            [-0.708,0.063,-0.394],
            [1.042,-1.467,-0.032],
            [0.452,-0.859,-0.127],
            [-0.317,-0.298,-0.245]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 8:
        c =[[1.095,-2.370,0.172],
            [1.310,-2.417,0.210],
            [0.981,-1.841,0.066],
            [0.027,-0.869,-0.131],
            [0.526,-1.357,-0.005],
            [0.340,-1.107,-0.074],
            [-2.326,0.900,-0.817],
            [1.417,-1.891,0.078],
            [1.075,-1.537,0.002],
            [-0.855,-0.085,-0.336],
            [-0.256,-0.361,-0.282],
            [0.543,-0.884,-0.147],
            [-1.123,0.324,-0.394]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 9:
        c = [[0.990,-2.618,0.256],
             [5.819,-8.836,1.332],
             [4.575,-7.617,1.042],
             [-0.321,-0.637,-0.203],
             [1.250,-2.109,0.131],
             [7.123,-10.442,1.678],
             [4.110,-7.814,1.436],
             [6.198,-9.274,1.629],
             [2.141,-2.530,0.177],
             [4.602,-7.965,1.389],
             [4.549,-7.808,0.154],
             [0.920,-1.286,-0.089],
             [2.277,-6.227,1.082]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 15:
        c = [[-3.070,2.052,-0.598],
             [-2.473,1.335,-0.425],
             [-0.183,-0.471,-0.090],
             [1.815,-2.099,0.217],
             [0.338,-1.148,0.065],
             [-0.171,0.073,-0.369],
             [-5.150,4.328,-1.305],
             [-2.798,2.824,-0.841],
             [-1.515,-0.873,0.090],
             [0,0,0],
             [-0.919,-1.228,0.154],
             [-2.063,-0.317,-0.017],
             [-0.158,-0.062,-0.066]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 16:
        c = [[0.013,-1.140,-0.067],
             [0.910,-1.698,0.060],
             [0.707,-1.328,-0.026],
             [-0.165,-0.419,-0.204],
             [1.110,-1.597,0.065],
             [1.330,-1.789,0.105],
             [-1.357,0.865,-0.494],
             [0.641,-0.888,-0.122],
             [1.396,-1.565,0.033],
             [-0.141,-1.563,0.270],
             [-1.314,-0.674,0.074],
             [0.337,-0.517,-0.184],
             [-0.235,-1.133,0.167]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 17:
        c = [[2.515,-3.541,0.464],
             [0.638,-1.266,-0.084],
             [1.207,-1.808,0.070],
             [0.221,-0.806,-0.133],
             [1.775,-2.190,0.172],
             [1.579,-1.996,0.128],
             [1.257,-1.951,0.110],
             [1.431,-1.539,-0.013],
             [1.667,-1.727,0.028],
             [-13.057,8.470,-1.803],
             [-4.712,0.551,-0.143],
             [1.404,-1.451,-0.011],
             [-4.416,0.113,-0.039]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 34:
        c = [[4.045,-6.292,0.942],
             [2.638,-5.475,0.823],
             [2.295,-5.240,0.808],
             [0.848,-4.078,0.542],
             [-0.328,-3.313,0.422],
             [-0.848,-2.888,0.319],
             [-2.585,-1.249,-0.049],
             [-0.715,-1.453,0.168],
             [-0.522,-3.526,0.732],
             [0,0,0],
             [-0.732,-1.004,0.061],
             [0,0,0],
             [1.057,-2.011,0.272]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 35:
        c = [[-0.229,-2.906,0.323],
             [-2.349,-1.366,0.034],
             [-3.636,-0.412,-0.147],
             [-6.151,1.949,-0.617],
             [-4.687,0.506,-0.294],
             [-5.486,1.221,-0.451],
             [-2.696,-1.058,0.063],
             [-1.871,-0.999,-0.063],
             [-1.596,-1.526,0.119],
             [0.830,-2.628,0.329],
             [0,0,0],
             [-4.372,0.369,-0.139],
             [-11.419,5.490,-1.028]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)
    if L == 53:
        c = [[1.650,-2.688,0.314],
             [3.884,-6.505,0.776],
             [0.160,-0.678,-0.150],
             [0.281,-0.799,-0.133],
             [1.588,-1.940,0.157],
             [-0.827,0.252,-0.382],
             [-0.932,0.239,-0.317],
             [-2.355,1.585,-0.519],
             [-2.261,1.489,-0.501],
             [1.746,-2.462,0.269],
             [-5.211,2.884,-0.944],
             [-3.620,1.323,-0.548],
             [2.644,0.438,-0.350]]
        for i in range(len(c)):
            if i+58 == Ln:
                a = c[i][0]
                b = c[i][1]
                d = c[i][2]
                plmn = a+b*X[0]+d*X[0]**2
                return np.exp(plmn)

ovlp = [rho(ZLn,ZL[i],P[i]) for i in range(len(P))]

#a0 is the Bohr radius in meters
a0 = 5.2917721090380E-11
#r^k radial integrals in units of a0^k (Ref: Edvardsson 1998 doi: )
def rk(k):
    start4 = [i for i in range(len(content)) if content[i].find('Input Lanthanide ion') != -1][0]
    stop4 = [i for i in range(len(content)) if content[i].find('Ion dE (4f^n -> 4f^n-1 5d)') != -1][0]
    rk = content[start4:stop4]
    rk = [i.split() for i in rk]
    prov2 = []
    for i in rk:
        for j in i:
            try:
                if isinstance(float(j), float) == True:
                    prov2.append(float(j))
            except:
                pass
    rk = prov2
    if k == 2:
        return rk[0]
    if k == 4:
        return rk[1]
    if k == 6:
        return rk[2]
#Ligand distances in units of a0
RL = [P[i][0]*1E-10/a0 for i in range(len(P))]
#Spherical Harmonics for each ligand

def SH(k,q,P):
    SH = []
    for i in range(len(P)):
        sh = spe.sph_harm(q,k,P[i][1],P[i][2])
        SH.append(np.conjugate(sh))
    return SH

#Bkq computation in cm^-1:

def Bkq(k,q,P):
    if -(k+1)<q<(k+1):
        c = np.sqrt((4*np.pi/(2*k+1)))*rk(k)
        AA = []
        for i in range(len(P)):
            a = SH(k,q,P)[i]*(1/RL[i])**(k+1)*ovlp[i]*(2*(1/(1+ovlp[i])))**(k+1)*g[i]
            AA.append(a)
        return 219474.63*c*np.sum(np.array(AA))
    else:
        return 0

#The factor 219474.63 converts the value from Harthrees to cm^-1

BB = []
for i in range(1,4):
    ii = 2*i
    for j in range(0,ii+1):
        LF = Bkq(ii,j,P)
        BB.append(('k=',ii,'q=',j,LF))

#BB: Complex Bkq values

CC = []
for i in range(1,4):
    ii = 2*i
    for j in range(0,ii+1):
        LF = np.abs(Bkq(ii,j,P))
        CC.append(('k=',ii,'q=',j,LF))

#CC: Absolute Bkq values

Parm = np.array(BB)
Abs = np.array(CC)

#Crystal Field Scalar Nv for the 7F1 level:

def kn(k,P):
    kk = []
    for q in range(-k,k+1):
        sss = (4*np.pi/(2*k+1))*np.abs(Bkq(k,q,P))**2
        kk.append(sss)
    return sum(kk)

#Def: Nv_j = in which the sum truncates with max k value = j


Nv_2 = np.sqrt(kn(2,P))
Nv_4 = np.sqrt(kn(2,P)+kn(4,P))
Nv_6 = np.sqrt(kn(2,P)+kn(4,P)+kn(6,P))

#Maximum splittings according to Journal of Alloys and Compounds 228 (1995) 41-44

dE_1 = np.sqrt(3*3**2/(3*4*5*np.pi))*0.53629*Nv_2
dE_2 = np.sqrt(3*5**2/(5*6*7*np.pi))*np.sqrt(0.4319*0.3938)*Nv_4
dE_3 = np.sqrt(3*7**2/(7*8*9*np.pi))*np.cbrt(0.2265*0.1818*0.2128)*Nv_6
dE_4 = np.sqrt(3*9**2/(9*10*11*np.pi))*np.cbrt(0.1478*0.6012*0.75438)*Nv_6
dE_5 = np.sqrt(3*11**2/(11*12*13*np.pi))*np.cbrt(0.7178*0.5123*0.7212)*Nv_6
dE_6 = np.sqrt(3*13**2/(13*14*15*np.pi))*np.cbrt(1.4982*0.7080*0.2178)*Nv_6


print('Note: To calculate Bk(-q) note that Bk(-q)=(-1)^q x Bkq*')
print('Bkq values in cm-1:')
print(Parm)
print('Absolute Bkq values in cm-1:')
print(Abs)
print('Calculating Secular Determinant and Splittings for selected (2S+1)LJ levels. This may take a while...')
#Secular determinant resolution:

def Ck(k):
    if k == 2:
        return -1.366
    if k == 4:
        return 1.128
    if k ==6:
        return -1.270
    else:
        return 0


def Uk(S,L,J,k):
    def RSME(X,S,L,J,k):
        def LSME(X,S,L,k):
            primes = [2,3,5,7,11,13,17,19,23,29,31]
            f = np.array([primes[i-1]**X[i] for i in range(1,len(X))])
            return X[0]*np.sqrt(np.prod(f))
        return np.real(complex((-1)**(S+L+J+k)*(2*J+1)*wgn.wigner_6j(J, J, k, L, L, S)*LSME(X,S,L,k)))
    if ZLn == 58:
        if k == 2:
            return 1
        if k == 4:
            return 1
        if k == 6:
            return 1
    if ZLn == 59 or ZLn == 69:
        if S == 1 and L == 1: #3P
            if k == 2:
                X = [-1,-1,2,0,-1]
        if S == 1 and L == 3: #3F
            if k == 2:
                X = [-1,0,-2]
            if k == 4:
                X = [-1,0,-2]
            if k == 6:
                X = [-1,0,-2]
        if S == 1 and L == 5: #3H
            if k == 2:
                X = [1,-1,-2,0,-1,1,1]
            if k == 4:
                X = [-1,2,-2,0,-1,0,1]
            if k == 6:
                X = [-1,0,-2,1,-1,0,0,1]
        if S == 0 and L == 2: #1D
            if k == 2:
                X = [-1,-1,-1,0,-2,2]
            if k == 4:
                X = [1,2,-2,1,-2,1]
            if k == 6:
                X = [1,1,0,1,-1,-1]
        if S == 0 and L == 4: #1G
            if k == 2:
                X = [1,0,3,0,-2,-1]
            if k == 4:
                X = [-1,0,0,0,-2,-2,1,0,0,2]
            if k == 6:
                X = [1,0,5,-1,-1,-2]
        if S == 0 and L == 6: #1I
            if k == 2:
                X = [1,-1,-1,2,0,-1,1]
            if k == 4:
                X = [1,3,-2,0,0,-2,1,1]
            if k == 6:
                X = [1,0,-1,0,-1,-2,0,1,1]
        if ZLn == 59:
            return RSME(X,S,L,J,k)
        if ZLn == 69:
            return -RSME(X,S,L,J,k)
    if ZLn == 60 or ZLn == 68:
        if S == 3/2 and L == 6: #4I
            if k == 2:
                X = [1,-1,-1,0,0,-1,1]
            if k == 4:
                X = [-1,1,-2,0,0,-2,1,1]
            if k == 6:
                X = [1,0,-1,2,-1,-2,0,1,1]
        if S == 3/2 and L == 3: #4F
            if k == 2:
                X = [-1,-2]
            if k == 4:
                X = [-1,-2]
            if k == 6:
                X = [-1,-2]
        if S == 3/2 and L == 4: #4G
            if k == 2:
                X = [-1,-2,1,0,-2,-1]
            if k == 4:
                X = [1,-2,4,2,-2,-2,1]
            if k == 6:
                X = [1,-2,1,1,-1,-2,2]
        if S == 1/2 and L == 4: #2G(1)
            if k == 2:
                X = [1,2,1,2,-4,-1]
            if k == 4:
                X = [-1,2,2,0,-4,-2,1,0,2]
            if k == 6:
                X = [1,2,1,1,-1,-2]
        if S == 1/2 and L == 5: #2H
            if k == 2:
                X = [1,-1,-4,0,-1,1,1]
            if k == 4:
                X = [-1,2,-4,0,-1,0,1]
            if k == 6:
                X = [-1,0,-4,1,-1,0,0,1]
        if S == 3/2 and L == 0: #4S
            if k == 2:
                return 0
            if k == 4:
                return 0 
            if k == 6:
                return 0
        if ZLn == 60:
            return RSME(X,S,L,J,k)
        if ZLn == 68:
            return -RSME(X,S,L,J,k)
    if ZLn == 62 or ZLn == 66:
        if S == 5/2 and L == 5: #6H
            if k == 2:
                X = [-1,-1,-2,0,-1,1,1]
            if k == 4:
                X = [1,2,-2,0,-1,0,1]
            if k == 6:
                X = [1,0,-2,1,-1,0,0,1]
        if S == 3/2 and L == 4: #4G
            if k == 2:
                X = [-1,-4,1,0,-2,-1]
            if k == 4:
                X = [1,-4,4,2,-2,-2,1]
            if k == 6:
                X = [1,-4,1,1,-1,-2,2]
        if S == 3/2 and L == 9: #4M
            if k == 2:
                X = [1,-1,0,0,0,0,0,-1,1]
            if k == 4:
                X = [1,-1,-2,0,0,0,0,-1,1,1]
            if k == 6:
                X = [-1,-2,2,2,-1,0,-2,-1,1,1]
        if S == 5/2 and L == 1: #6P
            if k == 2:
                X = [1,-1,2,0,-1]
            if k == 4:
                return 0
            if k == 6:
                return 0
        if S == 3/2 and L == 6: #4I
            if k == 2:
                X = [1,-3,-1,0,0,-1,1]
            if k == 4:
                X = [-1,-1,-2,0,0,-2,1,1]
            if k == 6:
                X = [1,-2,-1,2,-1,-2,0,1,1]
        if S == 3/2 and L == 3: #4F
            if k == 2:
                X = [-1,-4]
            if k == 4:
                X = [-1,-4]
            if k == 6:
                X = [-1,-4]
        if ZLn == 62:
            return RSME(X,S,L,J,k)
        if ZLn == 66:
            return -RSME(X,S,L,J,k)
    if ZLn == 63 or  ZLn == 65:
        if S == 3 and L == 3: #7F
            if k == 2:
                X = [-1,0]
            if k == 4:
                X = [-1,0]
            if k == 6:
                X = [-1,0]
        if S == 2 and L == 2: #5D
            if k == 2:
                X = [-1,-1,-3,0,-2,0,0,2]
            if k == 4:
                X = [-1,2,-4,1,-2,1]
            if k == 6:
                return 0
        if S == 2 and L == 8: #5L
            if k == 2:
                return 0
            if k == 4:
                X = [1,-3,0,0,1,-1,-1,1,1]
            if k == 6:
                X = [1,-3,2,0,-1,-1,-2,1,1,1]
        if S == 2 and L == 4:
            if k == 2:
                X = [1,-2,-1,0,-2,-1]
            if k == 4:
                X = [-1,-2,2,2,-2,-2,1]
            if k == 6:
                return 0
        if ZLn == 63:
            return RSME(X,S,L,J,k)
        if ZLn == 65:
            return -RSME(X,S,L,J,k)
    if ZLn == 64:
        if S == 7/2 and L == 0: #8S
            if k == 2:
                return 0
            if k == 4:
                return 0
            if k == 6:
                return 0
        if S == 5/2 and L == 1: #6P
            if k == 2:
                return 0
            if k == 4:
                return 0
            if k == 6:
                return 0
        if S == 5/2 and L == 6: #6I
            if k == 2:
                return 0
            if k == 4:
                return 0
            if k == 6:
                return 0
        if S == 5/2 and L == 2: #6D
            if k == 2:
                return 0
            if k == 4:
                return 0
            if k == 6:
                return 0
    if ZLn == 67:
        if S == 2 and L == 6: #5I
            if k == 2:
                X = [-1,-1,-1,0,0,-1,1]
            if k == 4:
                X = [1,1,-2,0,0,-2,1,1]
            if k == 6:
                X = [-1,0,-1,2,-1,-2,0,1,1]
        if S == 2 and L == 3: #5F
            if k == 2:
                X = [1,-2]
            if k == 4:
                X = [1,-2]
            if k == 6:
                X = [1,-2]
        if S == 2 and L == 0: #5S
            if k == 2:
                return 0
            if k == 4:
                return 0
            if k == 6:
                return 0
        return -RSME(X,S,L,J,k)
    if ZLn == 70:
        if k == 2:
            return -1
        if k == 4:
            return -1
        if k == 6:
            return -1
        
def Vij(S,L,J,i,j):
    soma = []
    for k in range(1,4):
        if -2*k<=(i-j)<=2*k:
            vij = (-1)**(J-i)*Ck(2*k)*Uk(S,L,J,2*k)*wgn.wigner_3j(J,2*k,J,-i,(i-j),j)*Bkq(2*k,(i-j),P)
            soma.append(vij)
    return complex(np.sum(soma))


def Secdet(S,L,J):
    primsecular = []
    if ZLn%2 != 0:
        for i in range(-J,J+1):
            for j in range(-J,J+1):
                primsecular.append(Vij(S,L,J,i,j))
    if ZLn%2 == 0:
        for i in range(int(2*J+1+0.01)):
            for j in range(int(2*J+1+0.01)):
                primsecular.append(Vij(S,L,J,(i-J),(j-J)))
    pp = np.array(primsecular)
    det = pp.reshape((int(2*J+1+0.01),int(2*J+1+0.01)))
    eigenvalues,eigenvectors = np.linalg.eigh(det)
    return np.real(eigenvalues)

#Writing the output .txt file
Spec_Term = {0 : 'S',
             1 : 'P',
             2 : 'D',
             3 : 'F',
             4 : 'G',
             5 : 'H',
             6 : 'I',
             7 : 'K',
             8 : 'L',
             9 : 'M',
             10 : 'N',
             11 : 'O',}
def JNb(S,L):
    if isinstance(S,int) == True:
        return [i for i in range((np.abs(L-S)),(L+S+1))]
    if isinstance(S,int) == False:
        return [i/2 for i in range((int(np.abs(2*L-2*S))),int((2*L+2*S+1))) if i/2 != int(i/2)]

def PTerm():
    print(f'{int(2*S+0.001)+1}{Spec_Term.get(L)}:')
    if ZLn%2 == 0:
        for i in range(len(JNb(S,L))):
            print(f'\n{int(2*S+0.001)+1}{Spec_Term.get(L)}{int(2*JNb(S,L)[i]+0.01)}/2')
            for j in range(int(2*JNb(S,L)[i]+0.01+1)):
                print(f'\n{int(2*S+0.001)+1}{Spec_Term.get(L)}{int(2*JNb(S,L)[i]+0.01)}/2  E{j+1}: ')
                print(str(np.sort(Secdet(S,L,JNb(S,L)[i]))[j]))
    if ZLn%2 != 0:
        for i in range(len(JNb(S,L))):
            print(f'\n{int(2*S+0.001)+1}{Spec_Term.get(L)}{int(JNb(S,L)[i]+0.01)}')
            for j in range(int(2*JNb(S,L)[i]+0.01+1)):
                print(f'\n{int(2*S+0.001)+1}{Spec_Term.get(L)}{int(JNb(S,L)[i]+0.01)}  E{j+1}: ')
                print(str(np.sort(Secdet(S,L,JNb(S,L)[i]))[j]))
    print('Calculating Next Term Symbol, please wait...')

def WTerm():
    otpt.write(f'{int(2*S+0.001)+1}{Spec_Term.get(L)}:')
    otpt.write('\n')
    if ZLn%2 == 0:
        for i in range(len(JNb(S,L))):
            otpt.write(f'\n{int(2*S+0.001)+1}{Spec_Term.get(L)}{int(2*JNb(S,L)[i]+0.01)}/2')
            for j in range(int(2*JNb(S,L)[i]+0.01+1)):
                otpt.write(f'\nE{j+1}: ')
                otpt.write(str(np.sort(Secdet(S,L,JNb(S,L)[i]))[j]))
            otpt.write('\n')
    if ZLn%2 != 0:
        for i in range(len(JNb(S,L))):
            otpt.write(f'\n{int(2*S+0.001)+1}{Spec_Term.get(L)}{int(JNb(S,L)[i]+0.01)}')
            for j in range(int(2*JNb(S,L)[i]+0.01+1)):
                otpt.write(f'\nE{j+1}: ')
                otpt.write(str(np.sort(Secdet(S,L,JNb(S,L)[i]))[j]))
            otpt.write('\n')

otpt = open('Bkq_output.log','w')
otpt.write('SOM Ligand Field Parameters')
otpt.write('\nBy Lucca Blois')
otpt.write('\nInstitute of Chemistry, São Paulo, Brazil')
otpt.write('\nNote: To calculate Bk(-q) note that Bk(-q)=(-1)^q x Bkq*')
otpt.write('\nNote: j is the imaginary unit')
otpt.write('\n')
otpt.write('\nLanthanide ion:')
otpt.write(f'{atheadder[0]}')
otpt.write(f'\n{atheadder[0]}-L Overlap Integrals:')
otpt.write('\n')
for i in range(len(P)):
    otpt.write(f'{atheadder[0]} - {atheadder[i+1]}: ')
    otpt.write(str(ovlp[i]))
    otpt.write('\n')
otpt.write('\nL Charge Factors from JOYSpectra:')
otpt.write('\n')
for i in range(len(g)):
    otpt.write(str(g[i]))
    otpt.write('\n')
otpt.write('\n')
otpt.write('\nBkq values in cm-1:')
otpt.write('\n')
otpt.write('\nk = 2:')
otpt.write('\n')
for i in range(0,3):
    otpt.write('q=')
    otpt.write(str(i))
    otpt.write(' ')
    otpt.write(str(Bkq(2,i,P)))
    otpt.write('\n')
otpt.write('\nk = 4:')
otpt.write('\n')
for i in range(0,5):
    otpt.write('q=')
    otpt.write(str(i))
    otpt.write(' ')
    otpt.write(str(Bkq(4,i,P)))
    otpt.write('\n')
otpt.write('\nk = 6:')
otpt.write('\n')
for i in range(0,7):
    otpt.write('q=')
    otpt.write(str(i))
    otpt.write(' ')
    otpt.write(str(Bkq(6,i,P)))
    otpt.write('\n')
otpt.write('\nAbsolute Bkq values in cm-1:')
otpt.write('\n')
otpt.write('\nk = 2:')
otpt.write('\n')
for i in range(0,3):
    otpt.write('q=')
    otpt.write(str(i))
    otpt.write(' ')
    otpt.write(str(np.abs(Bkq(2,i,P))))
    otpt.write('\n')
otpt.write('\nk = 4:')
otpt.write('\n')
for i in range(0,5):
    otpt.write('q=')
    otpt.write(str(i))
    otpt.write(' ')
    otpt.write(str(np.abs(Bkq(4,i,P))))
    otpt.write('\n')
otpt.write('\nk = 6:')
otpt.write('\n')
for i in range(0,7):
    otpt.write('q=')
    otpt.write(str(i))
    otpt.write(' ')
    otpt.write(str(np.abs(Bkq(6,i,P))))
    otpt.write('\n')
otpt.write('\n')
otpt.write('Crystal Field Scalar for J<=1   Nv = ')
otpt.write(str(Nv_2))
otpt.write('\n')
otpt.write('Crystal Field Scalar for J<=2   Nv = ')
otpt.write(str(Nv_4))
otpt.write('\n')
otpt.write('Crystal Field Scalar for J>=3 Nv = ')
otpt.write(str(Nv_6))
otpt.write('\n')
otpt.write('\nSplitting of Selected SLJ levels (with respect to transition centroid):')
otpt.write('\n')
otpt.write('\n')
if ZLn == 58 or ZLn == 70:
    S = 1/2
    L = 3
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
if ZLn == 59 or ZLn == 69:
    S = 1
    L = 1
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 1
    L = 3
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 1
    L = 5
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 0
    L = 2
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 0
    L = 4
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 0
    L = 6
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
if ZLn == 60 or ZLn == 68:
    S = 3/2
    L = 6
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 3/2
    L = 4
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 3/2
    L = 3
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 3/2
    L = 0
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 1/2
    L = 5
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 1/2
    L = 4
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
if ZLn == 61 or ZLn == 67:
    S = 2
    L = 6
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 2
    L = 3
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 2
    L = 0
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
if ZLn == 62 or ZLn == 66:
    S = 5/2
    L = 5
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 5/2
    L = 1
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 3/2
    L = 9
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 3/2
    L = 6
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 3/2
    L = 4
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 3/2
    L = 3
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
if ZLn == 63 or ZLn == 65:
    S = 3
    L = 3
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 2
    L = 8
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 2
    L = 4
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 2
    L = 2
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
if ZLn == 64:
    S = 7/2
    L = 0
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 5/2
    L = 6
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 5/2
    L = 2
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
    S = 5/2
    L = 1
    WTerm()
    PTerm()
    otpt.write('\n----------------------------------')
    otpt.write('\n')
otpt.write('\nGeometry in Cartesian Coordinates (x,y,z,)=')
otpt.write ('\n[0 0 0]')
for i in range(len(G)):
    otpt.write('\n')
    otpt.write(str(G[i]))
otpt.write('\n')
otpt.write('\nGeometry in Spherical Coordinates (R,phi,theta)=')
for i in range(len(P)):
    otpt.write('\n')
    otpt.write(str(P[i]))
otpt.close()
