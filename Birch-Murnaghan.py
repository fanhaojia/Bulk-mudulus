#from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
V0_exp=60.69845700   #experimental volume
V=[58.87750329,59.48448786,60.09147243,60.69845700,61.30544157,61.91242614,62.51941071]  # -3% to 3% of the volume
E_PBE=[-125.11011825,-125.12927075,-125.14059696,-125.14434100,-125.13995838,-125.08165327,-125.06523873]
E_HF=[-81.20209305,-81.21078860,-81.21077325,-81.20247642,-81.18583686,-81.13013596,-81.09860297]
E_RPA=[-155.9611463156,-156.0125216025,-156.0546101171,-156.0864741805,-156.1193270035,-156.1347441874,-156.1538288153]
#print V,E_RPA
V0=[0,0,0];E0=[0,0,0];B0=[0,0,0];B01=[0,0,0];a=[0,0,0]
a_exp=V0_exp**(1/3.0)
def E_murnaghan(V,E0,V0,B0,B01):
        E_V=E0+9.0/16.0*B0*V0*(((V0/V)**(2.0/3.0)-1)**3*B01+((V0/V)**(2.0/3.0)-1)**2*(6-4*(V0/V)**(2.0/3.0)))
        return E_V
a_guess_PBE=[-30.0, 60.0,100.0,-20.0 ]
a_guess_HF=[-80.0, 60.0,100.0,-20.0 ]
a_guess_RPA=[-150.0, 60.0,100.0,-20.0 ]
#E_RPA
popt0, pcov0=curve_fit(E_murnaghan,V,E_RPA,p0=a_guess_RPA)
print popt0
popt1, pcov1=curve_fit(E_murnaghan,V,E_PBE,p0=a_guess_PBE)
popt2, pcov2=curve_fit(E_murnaghan,V,E_HF,p0=a_guess_HF)
E0[0]=popt0[0];V0[0]=popt0[1];B0[0]=popt0[2];B01[0]=popt0[3]
E0[1]=popt1[0];V0[1]=popt1[1];B0[1]=popt1[2];B01[1]=popt1[3]
E0[2]=popt2[0];V0[2]=popt2[1];B0[2]=popt2[2];B01[2]=popt2[3]
A=len(E_RPA)
for i in range(A):
    E_RPA[i]=E_RPA[i]-E0[0]
    E_PBE[i]=E_PBE[i]-E0[1]
    E_HF[i]=E_HF[i]-E0[2]
for i in range(3):
    if i == 0:
        ss='RPA'
    elif i ==1:
        ss='PBE'
    elif i == 2:
        ss='HF'
    print 'E0_:'+ss+'    '+str(E0[i])
    print 'V0_:'+ss+'    '+str(V0[i])
    print 'B0_:'+ss+'    '+str(B0[i])
    print 'B01_:'+ss+'    '+str(B01[i])
    a[i]=V0[i]**(1/3.0);
    print 'Lattice constant('+ss+'):'+'    '+str(a[i])
xx=np.linspace(min(V),max(V),2000)
yy_0=E_murnaghan(xx,E0[0],V0[0],B0[0],B01[0])-E0[0]
yy_1=E_murnaghan(xx,E0[1],V0[1],B0[1],B01[1])-E0[1]
yy_2=E_murnaghan(xx,E0[2],V0[2],B0[2],B01[2])-E0[2]
plot1=plt.plot(V,E_RPA, '*',label='E_RPA')
plot2=plt.plot(V,E_PBE, '*',label='E_PBE')
plot3=plt.plot(V,E_HF, '*',label='E_HF')
plot4=plt.plot(xx,yy_0, 'r',color='black',label='RPA_curve_fit')
plot5=plt.plot(xx,yy_1, 'r',color='blue',label='SCAN_curve_fit')
plot6=plt.plot(xx,yy_2, 'r',color='orange',label='HF_curve_fit')

plot7=plt.axvline(V0[0],color='black',linestyle='dashed',label='a_RPA: '+str('%.3f' % a[0]))#'-*r',label='fitted Volume')
plot8=plt.axvline(V0_exp, color='red',linestyle='dashed',label='a_exp: '+str('%.3f' % a_exp))#'-*r',label='experimental Volume')
plt.xlabel('Volume ($\AA^3$)',fontsize=16)
plt.ylabel('Energy(eV/f.u.)',fontsize=16)
longx=max(V)-min(V)
#E_RPA
longy=[0,0,0]
longy[0]=max(E_RPA)-min(E_RPA)
longy[1]=max(E_PBE)-min(E_PBE)
longy[2]=max(E_HF)-min(E_HF)
longyy=max(longy)
plt.xlim(xmin=min(V)-0.1*longx,xmax=max(V)+0.15*longx)
plt.ylim(-0.1*longyy,1.2*longyy)
plt.legend(loc=1)
plt.title('SCAN',fontsize=18)
fig=plt.gcf()
fig.subplots_adjust(top=0.88,bottom=0.12,left=0.16,right=0.9,hspace=0.2,wspace=0.2)
plt.savefig('SCAN.png',dpi=300)
plt.show()



                  
