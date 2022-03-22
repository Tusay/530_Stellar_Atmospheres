# This code will be the function repository for useful functions from 530 Stellar Atmospheres.
# Include it in python with the import module feature (use sys to pull it from another directory if need be)

from astropy import units as u
from astropy import constants as const
from scipy.special import wofz
import numpy as np

# Calculate the Blackbody Planck Function for one frequency (nu) and one temp (T). 
# Optional coefficient A. 
# cgs units.
@u.quantity_input(T=u.K,nu=u.Hz)
def one_Planck(T,nu, A=1) -> (u.erg * u.cm**-2 * u.Hz**-1 * u.s**-1 * u.sr**-1):
    B_nu = A * (2 * const.h.cgs * nu**3 / (const.c.cgs)**2) / ((np.exp(const.h.cgs*nu/(const.k_B.cgs*T)))-1) / u.sr
    return B_nu

# Calculate the Blackbody Planck Function for an array of frequencies (nu) 
# and one temp (T) using the function one_Planck() defined above. 
# cgs units
def Planck_array_at_temp(T,nu_array,A=1) -> (u.erg * u.cm**-2 * u.Hz**-1 * u.s**-1 * u.sr**-1):
    B_nus=[]
    for nu in nu_array:
        B_nus.append(one_Planck(T,nu,A).value)
    return B_nus*(u.erg * u.cm**-2 * u.Hz**-1 * u.s**-1 * u.sr**-1)

# Calculate the Blackbody Planck Function for an array of wavenumbers (1/lambda) 
# and one temp (T) using the function one_Planck() defined above. 
# cgs units
@u.quantity_input(T=u.K,nu_array=1/u.m)
def Planck_from_wavenumbers(T,nu_array,A=1) -> (u.erg * u.cm**-2 * u.Hz**-1 * u.s**-1 * u.sr**-1):
    nu_array = (nu_array*const.c).to(u.Hz)
    B_nus = one_Planck(T,nu_array,A)
    return B_nus

# A numerical integrator for two arrays of equal length.
def my_integrator(x,y):
    dx = np.diff(x)
    ave_y = [(y[i]+y[i-1])/2 for i in range(1,len(y))]
    total = np.dot(dx,ave_y)
    return total 

# A wrapper for the integrator allowing the user to specify the function to obtain y values
def integrator_wrapper(x_min,x_max,density,x_units,function,*args):
    x_min*=x_units
    x_max*=x_units
    points = density*(x_max-x_min)
    x = np.logspace(x_min,x_max,points)
    x=x[1:]
    y = function(*args,x)
    answer = my_integrator(x,y)*const.c
    return answer

# A better integrator
def better_integrator(x,y):
    dx = np.diff(x)
    ave_y = np.diff(y)/2+y[:-1]
    total = np.dot(dx,ave_y)
    return total 

# A logspace integrator
def log_integrator(x,y):
    dlnx = np.diff(np.log(x))
    xy = x*y
    ave_xy = np.diff(xy)/2+xy[:-1]
    total = np.dot(dlnx,ave_xy)
    return total 

def better_integrator_wrapper(x_min,x_max,density,function,*args):
    points = round(density*(x_max-x_min))
    x = np.linspace(x_min,x_max,points)
    y = function(x, *args)
    return better_integrator(x,y)

def log_integrator_wrapper(x_min,x_max,density,function,*args):
    points = round(density*(x_max-x_min))
    x = 10**np.linspace(x_min,x_max,points)
    y = function(x, *args)
    return better_integrator(x,y)

def Voigt(x, x0, y0, a, sigma, gamma):
    return y0 + a * np.real(wofz((x - x0 + 1j*gamma)/sigma/np.sqrt(2))) / sigma /np.sqrt(2*np.pi)

def ave(array):
    return (array[1:]+array[:-1])/2

#Opacity Functions
# Note: valid only between 2600 to 113900 Angstroms
@u.quantity_input(T=u.K,lam=u.Angstrom)
def Kappa_Hff_2e(P_e, T, lam):
    theta = 5040*u.K/T
    lam=lam/u.Angstrom
    f0 = -2.2763 - 1.6850*np.log10(lam) + 0.76661*np.log10(lam)**2 - 0.053346*np.log10(lam)**3
    f1 = 15.2827 - 9.2846*np.log10(lam) + 1.99381*np.log10(lam)**2 - 0.142631*np.log10(lam)**3
    f2 = -197.789 + 190.266*np.log10(lam) - 67.9775*np.log10(lam)**2 + 10.6913*np.log10(lam)**3 - 0.625151*np.log10(lam)**4
    Kappa = 10**-26*P_e*10**(f0+f1*np.log10(theta)+f2*np.log10(theta)**2)
    return Kappa

# Note: valid only between 2250 to 15000 Angstroms
@u.quantity_input(T=u.K,lam=u.Angstrom)
def Kappa_Hbf_2e(P_e, T, lam):
    theta = 5040*u.K/T
    a0=1.99654
    a1=-1.18267*10**-5
    a2=2.64243*10**-6
    a3=-4.40524*10**-10
    a4=3.23992*10**-14
    a5=-1.39568*10**-18
    a6=2.78701*10**-23
    lam=lam/u.Angstrom
    alpha_bf=(a0+a1*lam+a2*lam**2+a3*lam**3+a4*lam**4+a5*lam**5+a6*lam**6)*10**-18
    Kappa = 4.158*10**-10*alpha_bf*P_e*theta**(5/2)*10**(0.754*theta)
    Kappa[Kappa<0]=0
    Kappa[np.where(Kappa==0)[0][0]:]=0
    return Kappa

@u.quantity_input(T=u.K,lam=u.Angstrom)
def Kappa_Hff(P_e, T, lam):
    theta = 5040*u.K/T
    lam=lam/u.Angstrom
    alpha_0=1.0449*10**-26
    I=13.59843401
    R=1.0968*10**-3
    X_lam = 1.2398*10**4/lam
    g_ff=1+0.3456/(lam*R)**(1/3)*(np.log10(np.exp(1))/(theta*X_lam)+0.5)
    Kappa = alpha_0*lam**3*g_ff*np.log10(np.exp(1))/2/theta/I*10**(-theta*I)
    return Kappa

@u.quantity_input(T=u.K,lam=u.Angstrom)
def Kappa_Hbf(P_e, T, lam):
    theta = 5040*u.K/T
    lam=lam/u.Angstrom
    I=13.59843401
    alpha_0=1.0449*10**-26
    Kappa=0
    R=1.0968*10**-3
    for n in range(1,100):
        g_bf = 1-0.3456/(lam*R)**(1/3)*(lam*R/n**2-0.5)
        if n==1:
            g_bf[lam>912]=0
        if n==2:
            g_bf[lam>3746]=0
        if n==3:
            g_bf[lam>8206]=0
        if n==4:
            g_bf[lam>14588]=0
        X = I*(1-1/n**2)
        Kappa += lam**3/n**3*g_bf*10**(-theta*X)
    return Kappa*alpha_0

@u.quantity_input(T=u.K,lam=u.Angstrom)
def NIST_Kappa_Total(df, Partition_Table, P_e, T, lam):
    Phi_val = (1+NIST_Phi(df,Partition_Table,'H',T)/P_e)
    theta=5040*u.K/T
    X_lam = 1.2398*10**4/(lam*(1/u.Angstrom))
    Kappa = ((Kappa_Hbf(P_e, T, lam)+Kappa_Hff(P_e, T, lam)+Kappa_Hbf_2e(P_e, T, lam))*(1-10**(-X_lam*theta))+Kappa_Hff_2e(P_e, T, lam))/Phi_val/P_e
    return Kappa

def A_sol(df,element):
    row = df[df.element==element]
    return row

@u.quantity_input(T=u.K)
def Partition(df,element, T):
    # print(element)
    if element=='H-' or element=='HII':
        return 1
    theta = 5040*u.K/T
    row = df[df.Element==element]
    thetas = list(df)[1:-1]
    U_strs = [val for val in list(row.reset_index(drop=True).iloc[0][1:-1])]
    nulls = []
    for i,U in enumerate(U_strs):
        if U == '-':
            nulls.append(i)
    for i in reversed(nulls):
        del thetas[i]
        del U_strs[i]
    Us = [float(u) for u in U_strs]
    U_r = np.interp(theta,thetas,Us)
    return 10**U_r

# Saha Phi(T) function using nist_ioniz.txt
@u.quantity_input(T=u.K)
def NIST_Phi(NIST_Table, PartitionTable, element, T):
    # print(element)
    if element=='H-':
        X = 0.754*u.eV
    elif element=='HII':
        X = 0*u.eV
    else:
        df=NIST_Table
        row = df[df.Element==element]
        X = row['Ionization_Energy'].item()*u.eV
    # k_B=const.k_B.cgs
    # h=const.h.cgs
    # m_e=const.m_e.cgs
    # pi=np.pi
    # Phi = (2*pi*m_e*k_B*T/h**2)**(3/2)*np.exp(-X/k_B/T)
    theta = 5040*u.K/T
    I = X/u.eV
    element0 = element
    if element=='H':
        element1 = element+'II'
    else:
        element1 = element+'+'
    if element == 'H-':
        element0 = 'H-'
        element1 = 'H'
    Pfrac = Partition(PartitionTable,element1,T)/Partition(PartitionTable,element0,T)
    Phi = 0.6665*5040**2.5*Pfrac*theta**(-5/2)*10**(-theta*I)
    return Phi

@u.quantity_input(T=u.K)
def NIST_Solve_Pe(Pg,elements,ATable,NIST_Table,PartitionTable,T):
    def Aj(j):
        if j == 'H-':
            j = 'H'
        return A_sol(ATable,j).A.item()
    if T.value > 30000:
        Pe = Pg/2
    else:
        Pe = np.sqrt(Pg*NIST_Phi(NIST_Table,PartitionTable,'H',T))
        # print(f'Initial Pe guess: {Pe}')
    count = 0
    Pei=0
    Pef=Pe
    while abs(Pef-Pei) > 1e-15:
        Pei=Pef
        num=0
        denom=0
        for j in elements:
            frac=NIST_Phi(NIST_Table,PartitionTable,j,T)/Pei
            num+=Aj(j)*frac/(1+frac)
            denom+=Aj(j)*(1+frac/(1+frac))
        Pef = Pg*num/denom
        if count == 100:
            print(f'Count exceeded. While loop broken.')
            print(f'Pef: {Pef}\nPei: {Pei}')
            break
        count+=1
    print(f'It took {count} iterations to calculate P_e')
    return Pef

