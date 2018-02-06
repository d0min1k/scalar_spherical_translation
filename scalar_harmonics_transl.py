# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:19:21 2018

@author: Dominic
"""

#import math
import numpy as np
import scipy.special
import matplotlib.pyplot as plt

#factorial from a table or approximation for big values
def tabfact(n):
    if n<21:
        return {0:1, 1:1, 2:2, 3:6, 4:24, 5:120, 6:720, 7:5040, 8:40320, 9:362880, 10:3628800, 11:39916800, 12:479001600, 13:6227020800, 14:87178291200, 15:1307674368000, 16:20922789888000, 17:355687428096000, 18:6402373705728000, 19:121645100408832000, 20:2432902008176640000}[n]
    else:
        return np.sqrt(2*np.pi*n)*np.exp(1/(12*n))*(n/np.e)**(n)

def wigner3ja(j1,j2,j3,m1,m2,m3):
    if 2*j1!=np.floor(2*j1) or 2*j2!=np.floor(2*j2) or 2*j3!=np.floor(2*j3) or 2*m1!=np.floor(2*m1) or 2*m2!=np.floor(2*m2) or 2*m3!=np.floor(2*m3):
        return 'non-integer input'
    if j1-m1!=np.floor(j1-m1) or j2-m2!=np.floor(j2-m2) or j3-m3!=np.floor(j3-m3):
        return 'must have same parity'
    
    if j3>j1+j2 or j3<abs(j1-j2):
        return 0;
    
    if abs(m1)>j1 or abs(m2)>j2:
        return 'out of bounds'
    
    if abs(m3)>j3:
        return 0;
        
    t1=j2-m1-j3
    t2=j1+m2-j3
    t3=j1+j2-j3
    t4=j1-m1
    t5=j2+m2
    
    tmin=max(0,max(t1,t2))
    tmax=min(t3,min(t4,t5))
    wig=0.
    
    for t in range(tmin,tmax+1):
        wig+=(-1)**t / (tabfact(t)*tabfact(t-t1)*tabfact(t-t2)*tabfact(t3-t)*tabfact(t4-t)*tabfact(t5-t))
    
    wig*= (-1.)**(j1-j2-m3)*np.sqrt(tabfact(j1+j2-j3)*tabfact(j1-j2+j3)*tabfact(-j1+j2+j3)*tabfact(j1+m1)*tabfact(j1-m1)*tabfact(j2+m2)*tabfact(j2-m2)*tabfact(j3+m3)*tabfact(j3-m3)/tabfact(j1+j2+j3+1))
    return wig

#special case for wigner-3j when all three bottom components are 0
def wigner3j0(j1,j2,j3):
#    if 2*j1!=m.floor(2*j1) or 2*j2!=m.floor(2*j2) or 2*j3!=m.floor(2*j3) or 2*m1!=m.floor(2*m1) or 2*m2!=m.floor(2*m2) or 2*m3!=m.floor(2*m3):
#        return 'non-integer input'
    if j1!=np.floor(j1) or j2!=np.floor(j2) or j3!=np.floor(j3):
        return 'must have same parity'
    
    if j3>j1+j2 or j3<abs(j1-j2):
        return 0;
            
    J=j1+j2+j3
    if J%2>0:
        return 0
    else:
        t0=J/2
        t1=t0-j1
        t2=t0-j2
        t3=t0-j3
        wig=np.sqrt(tabfact(2*t1))/tabfact(t1)
        wig*=np.sqrt(tabfact(2*t2))/tabfact(t2)
        wig*=np.sqrt(tabfact(2*t3))/tabfact(t3)        
        wig*=(-1)**(J/2)*tabfact(J/2)/np.sqrt(tabfact(J+1))        
        return wig

#spherical harmonic multiplied by bessel
def scalar_harm(m,n,x,y,z,k):
    rho=np.sqrt(x**2+y**2)
    R=np.sqrt(rho**2+z**2)
    theta=np.arctan2(rho,z)
    phi=np.arctan2(y,x)    
    bes=scipy.special.spherical_jn(n,k*R)
    leg=scipy.special.lpmv(m,n,np.cos(theta))
    eps=np.exp(1j*m*phi)
    return bes*leg*eps

def spherical_h1n(n,R):
    return scipy.special.spherical_jn(n,R)+1j*scipy.special.spherical_yn(n,R)

#spherical harmonic multiplied by hankel
def scalar_harm_h(m,n,x,y,z,k):
    rho=np.sqrt(x**2+y**2)
    R=np.sqrt(rho**2+z**2)
    theta=np.arctan2(rho,z)
    phi=np.arctan2(y,x)    
    bes=spherical_h1n(n,k*R)
    leg=scipy.special.lpmv(m,n,np.cos(theta))
    eps=np.exp(1j*m*phi)
    return bes*leg*eps
    
def fac_ratio(m,n):         # (n+m)!/(n-m)!
    if n<0 or np.abs(m)>n:
        print('lol')
        print(m,n)
#        return 0
    if m==0:
        return 1
    else:
        res=1.
        for i in np.arange(-np.abs(m)+1,np.abs(m)+1):
            res*=(n+i)
        if m<0:
            res=1./res
        return res

#coefficients for legendre addition theorem
def a_mnmunup(m,n,mu,nu,p):
    return (-1.)**(m+mu)*np.sqrt(fac_ratio(m,n)*fac_ratio(mu,nu)/fac_ratio(mu+m,p))*(2*p+1)*wigner3j0(n,nu,p)*wigner3ja(n,nu,p,m,mu,-m-mu)

def inner_coeff(m,n,mu,nu,x,y,z,k):
    A=0
    for p in np.arange(np.abs(nu-n),nu+n+1):
        if p>=np.abs(mu-m):
            harm=scalar_harm(m-mu,p,x,y,z,1)
            a=a_mnmunup(m,n,-mu,nu,p)
            A+=(-1.)**mu*1j**(nu+p-n)*(2*nu+1)*harm*a
    return A

def inner_coeff_h(m,n,mu,nu,x,y,z,k):
    A=0
    for p in np.arange(np.abs(nu-n),nu+n+1):
        if p>=np.abs(mu-m):
            harm=scalar_harm_h(m-mu,p,x,y,z,1)
            a=a_mnmunup(m,n,-mu,nu,p)
            A+=(-1.)**mu*1j**(nu+p-n)*(2*nu+1)*harm*a
    return A
    
def inner_coeff2(nn,m,n,mu,nu,x,y,z,k):
    rho=np.sqrt(x**2+y**2)
    R=np.sqrt(rho**2+z**2)
    theta=np.arctan2(rho,z)
    phi=np.arctan2(y,x)
    A=0.
    for p in np.arange(np.abs(nu-n),nu+n+1):
        if p>=np.abs(mu-m):
            bes=scipy.special.spherical_jn(n,k*R)
            harm=scipy.special.sph_harm(m,n,theta,phi)
            a=a_mnmunup(m,n,-mu,nu,p)
#            if np.isnan(a):
#                print(a,mu,nu,p)
            A+=(-1.)**mu*1j**(nu+p-n)*(2*nu+1)*bes*harm*a
#            if A==0:
#                print(bes,harm,a)
    return A
                    
def translate_scalar_harm_inner(nn,m,n,x,y,z,k,x1,x2,x3):
    b=np.zeros(np.shape(x),dtype=complex)
    for nu in np.arange(0,nn+1):
        for mu in np.arange(-nu,nu+1):
            c=inner_coeff(m,n,mu,nu,x1,x2,x3,k)
            d=scalar_harm(mu,nu,x,y,z,k)
            b+=c*d
    return b

def translate_scalar_harm_inner_h(nn,m,n,x,y,z,k,x1,x2,x3):
    b=np.zeros(np.shape(x),dtype=complex)
    for nu in np.arange(0,nn+1):
        for mu in np.arange(-nu,nu+1):
            c=inner_coeff_h(m,n,mu,nu,x1,x2,x3,k)
            d=scalar_harm(mu,nu,x,y,z,k)
            b+=c*d
    return b

def translate_scalar_harm_outer(nn,m,n,x,y,z,k,x1,x2,x3):
    b=np.zeros(np.shape(x),dtype=complex)
    for nu in np.arange(0,nn+1):
        for mu in np.arange(-nu,nu+1):
            for p in np.arange(np.abs(n-nu),nu+n+1):
                if p>=np.abs(mu-m):
                    harm=scalar_harm(m-mu,p,x,y,z,1)
                    a=a_mnmunup(m,n,-mu,nu,p)
                    A=(-1.)**mu*1j**(nu+p-n)*(2*nu+1)*harm*a
                    d=scalar_harm(mu,nu,x1,x2,x3,k)
                    b+=A*d
    return b

def translate_scalar_harm_outer_h(nn,m,n,x,y,z,k,x1,x2,x3):
    b=np.zeros(np.shape(x),dtype=complex)
    for nu in np.arange(0,nn+1):
        for mu in np.arange(-nu,nu+1):
            for p in np.arange(np.abs(n-nu),nu+n+1):
                if p>=np.abs(mu-m):
                    harm=scalar_harm_h(m-mu,p,x,y,z,1)
                    a=a_mnmunup(m,n,-mu,nu,p)
                    A=(-1.)**mu*1j**(nu+p-n)*(2*nu+1)*harm*a
                    d=scalar_harm(mu,nu,x1,x2,x3,k)
                    b+=A*d
    return b

#contour plot for real, imag or abs 
def conplt(x,y,z):
    plt.figure()
    #plt.imshow(np.abs(a))
#    plt.contour(x,y,np.real(z),30,colors='k',linewidths=.25)
#    plt.contourf(x,y,np.real(z),30)
#    plt.colorbar()
#    plt.figure()
#    plt.contour(x,y,np.imag(z),30,colors='k',linewidths=.25)
#    plt.contourf(x,y,np.imag(z),30)
#    plt.colorbar()
    plt.contour(x,y,np.abs(z),30,colors='k',linewidths=.25)
    plt.contourf(x,y,np.abs(z),30)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')

####### init
m=0
n=2
nn=7
        
xy=np.arange(-1,1,.01)*5
(x,y)=np.meshgrid(xy,xy)
z=0
k=1
RR=2

####### bessel translation
a=scalar_harm(m,n,x,y,z,k)
b=translate_scalar_harm_inner(nn,m,n,x,y,z,k,RR,0,0)
c=translate_scalar_harm_outer(nn,m,n,x,y,z,k,RR,0,0)
        
conplt(x,y,a)
conplt(x,y,c)

####### hankel translation
R=1
a=scalar_harm_h(m,n,x,y,0,1)
b=translate_scalar_harm_inner_h(nn,m,n,x,y,z,k,RR,0,0)
c=translate_scalar_harm_outer_h(nn,m,n,x,y,z,k,RR,0,0)
a[x**2+y**2<R**2]=0
b[(x+RR)**2+y**2<R**2]=0
c[(x+RR)**2+y**2<R**2]=0

b[x**2+y**2>RR**2]=0
c[x**2+y**2<RR**2]=b[x**2+y**2<RR**2]

conplt(x,y,a)
conplt(x,y,c)