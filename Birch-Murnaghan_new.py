#from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
V0_exp=60.69845700   #experimental volume
V=[58.87750329,59.48448786,60.09147243,60.69845700,61.30544157,61.91242614,62.51941071]  # -3% to 3% of the volume
E_PBE=[-125.11011825,-125.12927075,-125.14059696,-125.14434100,-125.13995838,-125.08165327,-125.06523873] #calcluated energy

V0=[0];E0=[0];B0=[0];B01=[0];a=[0]
a_exp=V0_exp**(1/3.0)
def E_murnaghan(V,E0,V0,B0,B01):
        E_V=E0+9.0/16.0*B0*V0*(((V0/V)**(2.0/3.0)-1)**3*B01+((V0/V)**(2.0/3.0)-1)**2*(6-4*(V0/V)**(2.0/3.0)))
        return E_V
a_guess_PBE=[-30.0, 60.0,100.0,-20.0 ]  ##you need address parameters to fit the E_V curve
popt0, pcov0=curve_fit(E_murnaghan,V,E_PBE,p0=a_guess_PBE)
E0[0]=popt0[0];V0[0]=popt0[1];B0[0]=popt0[2];B01[0]=popt0[3]

A=len(E_PBE)
for i in range(A):
    E_PBE[i]=E_PBE[i]-E0[0]

ss='PBE'
print ('E0_:'+ss+'    '+str(E0[0]))
print ('V0_:'+ss+'    '+str(V0[0]))
print ('B0_:'+ss+'    '+str(B0[0]))
print ('B01_:'+ss+'    '+str(B01[0]))
xx=np.linspace(min(V),max(V),2000)
yy_0=E_murnaghan(xx,E0[0],V0[0],B0[0],B01[0])-E0[0]

plot1=plt.plot(V,E_PBE, '*',label='E_PBE')
plot4=plt.plot(xx,yy_0, 'r',color='black',label='PBE_curve_fit')

plt.xlabel('Volume ($\AA^3$)',fontsize=16)
plt.ylabel('Energy(eV/f.u.)',fontsize=16)
longx=max(V)-min(V)
longyy=max(E_PBE)-min(E_PBE)
plt.xlim(xmin=min(V)-0.1*longx,xmax=max(V)+0.15*longx)
plt.ylim(-0.1*longyy,1.2*longyy)
plt.legend(loc=1)
plt.title('E_V_PBE',fontsize=18)
fig=plt.gcf()
fig.subplots_adjust(top=0.88,bottom=0.12,left=0.16,right=0.9,hspace=0.2,wspace=0.2)
plt.savefig('E_V_PBE.png',dpi=300)
plt.show()



                  
