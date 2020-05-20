# -*- coding: utf-8 -*-
"""
Created on Sun May 17 16:31:06 2020

@author: ariannagasparri
@coauthor: pa
"""

import numpy as np
import math as mt
import datetime as dt
import matplotlib
from matplotlib import pyplot as plt
#matplotlib.use("Agg")
#import matplotlib.animation as animation

# lunghezza dei link
a_1 = 1.3
a_2 = 1.6
# masse dei link approssimate
m_1 = 1500
m_2 = 1350
# accelerazione di gravità
g = 9.81

dt = 0.01
N = 700

ks = 1
kv = 10

# momenti di inerzia principali rispetto ai frame locali
I_xx_1 = (m_1*((a_1/2)**2))/3
I_zz_1 = (m_1*((a_1/2)**2))/3

I_xx_2 = (m_2*((a_2/2)**2))/3
I_zz_2 = (m_2*((a_2/2)**2))/3

# Matrice di Inerzia
def b11(q2): 
	return (m_1 * (a_1**2)/4) + (m_2 * (a_1**2 + (a_2**2)/4 + a_1*a_2*mt.cos(q2))) + I_zz_1 + I_zz_2

def b12(q2):
	return m_2 * ((a_1 * a_2/2 * mt.cos(q2)) + (a_2**2)/4) + I_zz_2

def b21(q2):
	return m_2 * ((a_1 * a_2/2 * mt.cos(q2)) + (a_2**2)/4) + I_zz_2

def b22():
	return m_2 * (a_2**2)/4 + I_zz_2
	
def B(q2):
	return np.array([[b11(q2), b12(q2)], [b21(q2), b22()]])

# Matrice di Criolis
def c11(q2, q2_dot):
	return (-0.5 * m_2 * a_1 * a_2 * mt.sin(q2)) * q2_dot

def c12(q2, q1_dot, q2_dot):
	return (-0.5 * m_2 * a_1 * a_2 * mt.sin(q2)) * (q1_dot + q2_dot)

def c21(q2, q1_dot):
	return (m_2 * a_1 * a_2 * mt.sin(q2)) * q1_dot
	
def C(q2, dq1, dq2):
	return np.array([[c11(q2, dq2), c12(q2, dq1, dq2)], [c21(q2, dq1), 0]])

# Matrice Gravitazionale
def g1(q1, q2):
	return g * (m_1*(a_1/2)*mt.cos(q1) + m_2*(a_1*mt.cos(q1) + a_2/2*mt.cos(q1 + q2)))

def g2(q1, q2):
	return g * m_2*(a_2/2)*mt.cos(q1 + q2)

def G(q1, q2):
	return np.array([[g1(q1, q2)], [g2(q1, q2)]])
	

def dynamic(q_, dq_):
	dq = np.array([dq_]).T
	q = np.array([q_]).T
	sign_vector = np.sign(dq)
	#Fv = np.matmul(np.matrix([[kv,0],[0,kv]]),dq) #dipende dalla velocità
	#Fs = np.matmul(np.matrix([[ks,0],[0,ks]]), sign_vector) #non dipende dalla velocità
	ret = np.matmul(np.linalg.inv(B(q[1])), 0 - np.matmul(C(q[1], dq[0], dq[1]), dq) - G(q[0],q[1]))



	return ret
	
q1 = mt.pi/2
q2 = mt.pi/6
dq1 = 0
dq2 = 0
q = np.zeros([2,1])
dq = np.zeros([2,1])
ddq = np.zeros([2,1])

q[0,0] = q1
q[1,0] = q2
dq[0,0] = dq1
dq[1,0] = dq2

print(B(q[1,0]))
print(C(q[1,0], dq[0,0] , dq[1,0] ))
print(G(q[0,0], q[1,0]))



for k in range(0,N):
	
	ddq_temp = np.array(dynamic(q[:,k], dq[:,k]))
	try:
		ddq = np.concatenate((ddq,ddq_temp), axis=1)
		dq_temp = np.transpose([dq[:,k]]) + np.transpose([ddq[:,k]*dt])
		dq = np.concatenate((dq,dq_temp),axis=1)

		q_temp = np.transpose([q[:,k]]) + np.transpose([dq[:,k]*dt])
		q = np.concatenate((q,q_temp),axis=1)
	except:
		break

time = np.arange(0., N*dt+dt, dt)

q0_list = q[0,:]
q1_list = q[1,:]
dq0_list = dq[0,:]
dq1_list = dq[1,:]

plt.figure()
plt.plot(time, q0_list)
plt.plot(time, q1_list)
plt.show(block=False)

plt.figure()
plt.plot(time, dq0_list)
plt.plot(time, dq1_list)
plt.show()

