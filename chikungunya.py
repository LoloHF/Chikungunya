#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 11:34:01 2022

@author: laurianepierrotdeseilligny
"""
import numpy as np
import matplotlib.pyplot as plt
import math

"""

b=5
Ke=1000
s=0.2
d=0.6
Kl=500
sl=0.3
dl=0.6
dm=0.7
r=(b/(s+d))*(s/(sl+dl))*(sl/dm)

#### ELA system is independant from the folowing parameters ####
bh=0.0000457
beth=0.75
gamh=0.1428
betm=0.5
R0=(betm*beth)/(dm*(gamh+bh))
gamE=1+((s+d)*dm*Ke)/(b*sl*Kl)
gaml=1+((sl+dl)*Kl)/(s*Ke)
f=(1-(1/r))
Xx=np.zeros((1,3))
Xx=[f*Ke/gamE,f*Kl/gaml,f*(sl/dm)*(Kl/gaml)]

SHx=bh/(beth+bh) + beth/((beth+bh)*R0)
IHx=dm*bh*(R0-1)/(betm*(beth+bh))
IMx=bh*(R0-1)/(beth+bh*R0)
####

if R0 <=1 :
	print("The system  (E,L,A);(Sh,Ih,Im) has {}; (1,0,0) as globally asymptotically stable equilibrium".format(Xx)) 
else :
	print("The system  (E,L,A);(Sh,Ih,Im) has ({},{},{}); (1,0,0) has a unique endemic equilibrium".format(SHx,IHx,IMx))

if r <=1 :
	print("The system (E,L,A) only has mosquitoe-free equilibrium (0,0,0)") 
else :
	print("The system (E,L,A) has {} as unique endemic equilibrium".format(Xx)) 
	


"""
b=6
Ke=1000
s=0.5
d=0.2
Kl=500
sl=0.5
dl=0.25
dm=0.25
r=(b/(s+d))*(s/(sl+dl))*(sl/dm)
bh=0.0000457
beth=0.75
gamh=0.1428
betm=0.5
R0=(betm*beth)/(dm*(gamh+bh))
gamE=1+((s+d)*dm*Ke)/(b*sl*Kl)
gaml=1+((sl+dl)*Kl)/(s*Ke)
f=(1-(1/r))
Xx=np.zeros((1,3))
Xx=[f*Ke/gamE,f*Kl/gaml,f*(sl/dm)*(Kl/gaml)]

SHx=bh/(beth+bh) + beth/((beth+bh)*R0)
IHx=dm*bh*(R0-1)/(betm*(beth+bh))
IMx=bh*(R0-1)/(beth+bh*R0)


if R0 <=1 :
	print("The system  (E,L,A);(Sh,Ih,Im) has {}; (1,0,0) as globally asymptotically stable equilibrium".format(Xx)) 
else :
	print("The system  (E,L,A);(Sh,Ih,Im) has ({},{},{}); (1,0,0) has a unique endemic equilibrium".format(SHx,IHx,IMx))

if r <=1 :
	print("The system (E,L,A) only has mosquitoe-free equilibrium (0,0,0)") 
else :
	print("The system (E,L,A) has {} as unique endemic equilibrium".format(Xx)) 


def mosquit(X0):
	X=np.array(X0)
	V=np.empty((1,3))
	E=X[0]
	L=X[1]
	A=X[2]
	V[0][0]=b*A*(1-E/Ke)-(s+d)*E
	V[0][1]=s*E*(1-L/Kl)-(sl+dl)*L
	V[0][2]=sl*L-dm*A
	return(V)


#CAREFUL gamh, not clear : cf article gam in the eq, gamh in the plot fig5
def human(H0):	
	Lx=Xx[1]
	Ax=Xx[2]
	H=np.array(H0)
	N=np.empty((1,3))
	Sh=H[0]
	Ih=H[1]
	Im=H[2]
	N[0][0]=-(bh+beth*Im)*Sh+bh
	N[0][1]=beth*Im*Sh-(gamh+bh)*Ih
	N[0][2]=-(sl*(Lx/Ax)+betm*Ih)*Im+betm*Ih
	return(N)
	

def norm(Y):
	Z=math.sqrt(np.sum((Y**2)))
	return(Z)

def vec(X1,V,h):
	X2=X1+h*V
	return(X2)

def ligneM(X0,n):
	Z=np.empty((n,3))
	Z[0]=X0
	V=np.empty((n,3))
	for i in range(1,n):
		V[i-1]=mosquit(Z[i-1])
		Z[i]=vec(Z[i-1],V[i-1],norm(V[i-1])**-1)
	return(Z,V)

def ligneH(H0,n):
	Z=np.empty((n,3))
	Z[0]=H0
	N=np.empty((n,3))

	for i in range(1,n):
		N[i-1]=human(Z[i-1])
		Z[i]=vec(Z[i-1],N[i-1],0.0001*norm(N[i-1])**-1)
	return(Z,N)
	
	
def ELA(X0,n):
	Z,V=ligneM(X0,n)
	#E=Z[:,0]
	#L=Z[:,1]
	#A=Z[:,2]
	return(Z)

def ShIhIm(H0,n):
	Z,V=ligneH(H0,n)
	#E=Z[:,0]
	#L=Z[:,1]
	#A=Z[:,2]
	return(Z)


def X0(k,j):
	LE_h=np.empty((k,3))
	LE_b=np.empty((k,3))
	LE_d=np.empty((k,3))
	LE_g=np.empty((k,3))
	
	AE_h=np.empty((k,3))
	AE_b=np.empty((k,3))
	AE_d=np.empty((k,3))
	AE_g=np.empty((k,3))
	
	AL_h=np.empty((k,3))
	AL_b=np.empty((k,3))
	AL_d=np.empty((k,3))
	AL_g=np.empty((k,3))
	
	for i in range(k):
		AL_h[i,:]=[0,i*200,1200] #j=0
		AL_b[i,:]=[0,i*200,1]
		AL_d[i,:]=[0,1200,i*200]
		AL_g[i,:]=[0,1,i*200]
		
		AE_h[i,:]=[i*200,0,1200]
		AE_b[i,:]=[i*200,0,1] #j=1
		AE_g[i,:]=[1,0,i*200]
		AE_d[i,:]=[1200,0,i*200]
		
		LE_h[i,:]=[i*200,1200,0]
		LE_b[i,:]=[i*200,1,0]
		LE_g[i,:]=[1,i*200,0]
		LE_d[i,:]=[1200,i*200,0] #j=2
	if j==0:
		return(AL_h,AL_b,AL_d,AL_g)
	elif j==1:
		return(AE_h,AE_b,AE_d,AE_g)
	else:
		return(LE_h,LE_b,LE_d,LE_g)

def H0(k):
	IhSh=np.empty((k,3))
	ImSh=np.empty((k,3))
	ImIh=np.empty((k,3))
	for i in range(k): #YX=0,X changes,Y constant
		ImIh[i,:]=[0,0.1+i*0.01,0.000001] #j=0
		ImSh[i,:]=[0.1+i*0.01,0,0.000001] #j=1
		IhSh[i,:]=[0.1+i*0.01,0.000001,0] #j=2
	return(ImIh,ImSh,IhSh)

def portraitM(k,n,j):
	CI=X0(k,j)
	x=np.empty((n,1))
	y=np.empty((n,1))
	for c in CI:
		for l in range(k):
			XY=np.delete(ELA(c[l],n),j,axis=1)
			x=XY[:,0]
			y=XY[:,1]
			plt.plot(x,y)
		
	plt.show()
	if j==0:
		plt.savefig('PortraitdephaseMoustique0.png',format='png',bbox_inches='tight')
		plt.xlabel('L')
		plt.ylabel('A')
	elif j==1:
		plt.savefig('PortraitdephaseMoustique1.png',format='png',bbox_inches='tight')
		plt.xlabel('E')
		plt.ylabel('A')
	else:
		plt.savefig('PortraitdephaseMoustique2.png',format='png',bbox_inches='tight')
		plt.xlabel('E')
		plt.ylabel('L')

def portraitH(k,n,j):
	#6,10000,1 
	CI=H0(k)[j]
	x=np.empty((n,1))
	y=np.empty((n,1))
	for i in CI:
		XY=np.delete(ShIhIm(i,n),j,axis=1)
		x=XY[:,0]
		y=XY[:,1]
		plt.plot(x,y)
	plt.show()	
	if j==0:
		plt.savefig('PortraitdephaseHumain0.png',format='png',bbox_inches='tight')
		plt.xlabel('Ih')
		plt.ylabel('Im')
	elif j==1:
		plt.savefig('PortraitdephaseHumain1.png',format='png',bbox_inches='tight')
		plt.xlabel('Sh')
		plt.ylabel('Im')
	else:
		plt.savefig('PortraitdephaseHumain2.png',format='png',bbox_inches='tight')
		plt.xlabel('Sh')
		plt.ylabel('Ih')
	
def TEMPS(start,n):
	H0=np.ones((1,3))*0.00001
	Sh=ShIhIm(H0,n)[start:,1]
	Ih=ShIhIm(H0,n)[start:,1]
	Im=ShIhIm(H0,n)[start:,2]
	X=np.linspace(start,n-1,n-start)
	#plt.plot(X,Sh)
	plt.plot(X,Ih)
	plt.plot(X,Im)
	plt.show()
	plt.xlabel('Temps (jour)')
	plt.ylabel('Im,Ih')
	plt.savefig('Temps.png',format='png',bbox_inches='tight')


"""Xt1=[400,1200,0]
m=10000
E1,L1=LvsE(Xt1, m)

Xt2=[600,1200,0]
m=10000
E2,L2=LvsE(Xt2, m)

Xt3=[600,1200,0]
m=10000
E2,L2=LvsE(Xt2, m)

plt.plot(L1,E1)
plt.plot(L2,E2)
plt.show()"""






	