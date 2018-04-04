#from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
V0_exp=64.048012001
V=[62.12657164097,62.76705176098,63.40753188099,64.04801200100,64.68849212101,65.32897224102,65.96945236103]
E_RPA=[-155.0951693564,-155.1164076206,-155.1193012099,-155.1220723272,-155.1134407981,-155.1054506802,-155.0830905321]
#E_PBE=[-125.11011825,-125.12927075,-125.14059696,-125.14434100,-125.13995838,-125.08165327,-125.06523873]
#E_HF=[-81.20209305,-81.21078860,-81.21077325,-81.20247642,-81.18583686,-81.13013596,-81.09860297]

nnn=3
V0=[0]*nnn;E0=[0]*nnn;K0=[0]*nnn;K01=[0]*nnn;a=[0]*nnn;longy=[0]*nnn
a_exp=V0_exp**(1/3.0)
def E_murnaghan(V,E0,V0,K0,K01):
    if K01 !=0 and K01 !=1.0 and V0>0:
        E_V=E0+K0*V0*((V/V0)**(1.0-K01)/(K01*(K01-1.0))+(V/V0)/K01-1.0/(K01-1.0))
        return E_V
    else:
        return 0
a_guess_RPA=[-150, 60,3,7 ]
#a_guess_PBE=[-120, 60,3,7 ]
#a_guess_HF=[-80, 60,3,7 ]

#E_RPA
popt0, pcov0=curve_fit(E_murnaghan,V,E_RPA,p0=a_guess_RPA)
#popt1, pcov1=curve_fit(E_murnaghan,V,E_PBE,p0=a_guess_PBE)
#popt2, pcov2=curve_fit(E_murnaghan,V,E_HF,p0=a_guess_HF)
E0[0]=popt0[0];V0[0]=popt0[1];K0[0]=popt0[2];K01[0]=popt0[3]
#E0[1]=popt1[0];V0[1]=popt1[1];K0[1]=popt1[2];K01[1]=popt1[3]
#E0[2]=popt2[0];V0[2]=popt2[1];K0[2]=popt2[2];K01[2]=popt2[3]
A=len(E_RPA)
for i in range(A):
    E_RPA[i]=E_RPA[i]-E0[0]
#    E_PBE[i]=E_PBE[i]-E0[1]
#    E_HF[i]=E_HF[i]-E0[2]
for i in range(1):
    if i == 0:
        ss='RPA'
#    elif i ==1:
#        ss='PBE'
#    elif i == 2:
#        ss='HF'
    print 'E0_:'+ss+'    '+str(E0[i])
    print 'V0_:'+ss+'    '+str(V0[i])
    print 'K0_:'+ss+'    '+str(K0[i])
    print 'K01_:'+ss+'    '+str(K01[i])
    a[i]=V0[i]**(1/3.0);
    print 'Lattice constant('+ss+'):'+'    '+str(a[i])
xx=np.linspace(min(V),max(V),2000)
yy_0=E_murnaghan(xx,E0[0],V0[0],K0[0],K01[0])-E0[0]
#yy_1=E_murnaghan(xx,E0[1],V0[1],K0[1],K01[1])-E0[1]
#yy_2=E_murnaghan(xx,E0[2],V0[2],K0[2],K01[2])-E0[2]
plot1=plt.plot(V,E_RPA, '*',label='E_RPA')
#plot2=plt.plot(V,E_PBE, '*',label='E_PBE')
#plot3=plt.plot(V,E_HF, '*',label='E_HF')
plot4=plt.plot(xx,yy_0, 'r',color='black',label='RPA_curve_fit')
#plot5=plt.plot(xx,yy_1, 'r',color='blue',label='SCAN_curve_fit')
#plot6=plt.plot(xx,yy_2, 'r',color='orange',label='HF_curve_fit')

plot7=plt.axvline(V0[0],color='black',linestyle='dashed',label='a_RPA: '+str('%.3f' % a[0]))#'-*r',label='fitted Volume')
plot8=plt.axvline(V0_exp, color='red',linestyle='dashed',label='a_exp: '+str('%.3f' % a_exp))#'-*r',label='experimental Volume')
plt.xlabel('Volume ($\AA^3$)',fontsize=16)
plt.ylabel('Energy(eV/f.u.)',fontsize=16)
longx=max(V)-min(V)
#E_RPA

longy[0]=max(E_RPA)-min(E_RPA)
#longy[1]=max(E_PBE)-min(E_PBE)
#longy[2]=max(E_HF)-min(E_HF)
longyy=max(longy)
plt.xlim(xmin=min(V)-0.1*longx,xmax=max(V)+0.15*longx)
plt.ylim(-0.1*longyy,1.2*longyy)
plt.legend(loc=1)
plt.title('SCAN',fontsize=18)
fig=plt.gcf()
fig.subplots_adjust(top=0.88,bottom=0.12,left=0.16,right=0.9,hspace=0.2,wspace=0.2)
plt.savefig('SCAN.png',dpi=300)
plt.show()



                  
