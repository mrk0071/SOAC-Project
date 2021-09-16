# -*- coding: utf-8 -*-
"""
@author: Ruoke Meng

Hello Python!

"""

import numpy as np
import matplotlib.pyplot as plt

L = 2500000
dx = 25000
nx_len = 100
dt = 100
nt_len = 3000

## ----------------- ##
C = np.zeros((nt_len,nx_len))
C[0,int(1125/25):int(1375/25+1)] = 10
u = 10

for nt in range(1,nt_len):
    for nx in range(nx_len):
        C[nt,nx] = C[nt-1,nx]-u*(C[nt-1,nx]-C[nt-1,nx-1])/dx*dt
# ~~~~~~~ #
fig = plt.figure(figsize=(7,6))
ax3 = plt.axes(projection='3d')
ax3.set_xlabel("x [km]", fontsize=14, labelpad=8.5)
ax3.set_ylabel("time [s]", fontsize=14, labelpad=8.5)
ax3.set_zlabel("amplitude", fontsize=14, labelpad=8.5)
ax3.tick_params(labelsize=14)
xx = np.arange(0,L,dx)/1000
yy = np.arange(0,nt_len,1)
X, Y = np.meshgrid(xx, yy)
ax3.plot_surface(X,Y,C,rstride = 2, cstride = 2, cmap='Reds')
plt.show()
# ~~~~~~~ #
plt.figure(figsize=(7,6))
plt.xlabel("x [km]", fontsize=16)
plt.ylabel("time [s]", fontsize=16)
plt.title("Euler forward in time and upwind in space",fontsize=16)
plt.xticks(np.arange(0,101,20),np.arange(0,101,20)*25,fontsize=16)
plt.yticks(fontsize=16)
plt.imshow(C, cmap="Reds", aspect='auto' ,origin = 'lower')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=16)
plt.show()

# ~~~~~~~ #
fig = plt.figure(figsize=(13,5))
fig.suptitle("Euler forward in time and upwind in space",fontsize=18)
plt.subplots_adjust(wspace = 0.5)
ax3 = fig.add_subplot(1,2,1,projection='3d')
ax3.set_xlabel("x [km]", fontsize=16, labelpad=8.5)
ax3.set_ylabel("time [s]", fontsize=16, labelpad=8.5)
ax3.set_zlabel("amplitude", fontsize=16, labelpad=8.5)
ax3.tick_params(labelsize=14)
xx = np.arange(0,L,dx)/1000
yy = np.arange(0,nt_len,1)
X, Y = np.meshgrid(xx, yy)
ax3.plot_surface(X,Y,C,rstride=2, cstride=2, cmap='Reds')

ax = fig.add_subplot(1,2,2)
im = ax.imshow(C, cmap="Reds", aspect='auto' ,origin = 'lower')
ax.set_xlabel("x [km]", fontsize=16)
ax.set_ylabel("time [s]", fontsize=16)
plt.xticks(np.arange(0,100,20),np.arange(0,100,20)*25,fontsize=14)
ax.tick_params(labelsize=14)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax).ax.tick_params(labelsize=14)

plt.savefig("ex2_(1).png",bbox_inches="tight",dpi = 300)
plt.show()

## ----------------- ##
C = np.zeros((nt_len,nx_len))
C[0,int(1125/25):int(1375/25+1)] = 10
u = 10
# ~~~~~~~ #
for nt in range(1,nt_len):
    for nx in range(nx_len):
        C[nt,nx-1] = C[nt-1,nx-1] + (-u*((C[nt-1,nx]-C[nt-1,nx-2])/(2*dx)) + \
                                      (u**2*dt/2)*((C[nt-1,nx]-2*C[nt-1,nx-1]+C[nt-1,nx-2])/dx**2))*dt
# ~~~~~~~ #
fig = plt.figure(figsize=(7,6))
ax3 = plt.axes(projection='3d')
ax3.set_xlabel("x [km]", fontsize=14, labelpad=8.5)
ax3.set_ylabel("time [s]", fontsize=14, labelpad=8.5)
ax3.set_zlabel("amplitude", fontsize=14, labelpad=8.5)
ax3.tick_params(labelsize=14)
xx = np.arange(0,L,dx)/1000
yy = np.arange(0,nt_len,1)
X, Y = np.meshgrid(xx, yy)
ax3.plot_surface(X,Y,C,rstride = 2, cstride = 2, cmap='Reds')
plt.show()
# ~~~~~~~ #
plt.figure(figsize=(7,6))
plt.xlabel("x [km]", fontsize=16)
plt.ylabel("time [s]", fontsize=16)
plt.title("Lax-Wendroff scheme",fontsize=16)
plt.xticks(np.arange(0,101,20),np.arange(0,101,20)*25,fontsize=16)
plt.yticks(fontsize=16)
plt.imshow(C, cmap="Reds", aspect='auto' ,origin = 'lower')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=16)
plt.show()

# ~~~~~~~ #
fig = plt.figure(figsize=(13,5))
fig.suptitle("Lax-Wendroff scheme",fontsize=18)
plt.subplots_adjust(wspace = 0.5)
ax3 = fig.add_subplot(1,2,1,projection='3d')
ax3.set_xlabel("x [km]", fontsize=16, labelpad=8.5)
ax3.set_ylabel("time [s]", fontsize=16, labelpad=8.5)
ax3.set_zlabel("amplitude", fontsize=16, labelpad=8.5)
ax3.tick_params(labelsize=14)
xx = np.arange(0,L,dx)/1000
yy = np.arange(0,nt_len,1)
X, Y = np.meshgrid(xx, yy)
ax3.plot_surface(X,Y,C,rstride=2, cstride=2, cmap='Reds')

ax = fig.add_subplot(1,2,2)
im = ax.imshow(C, cmap="Reds", aspect='auto' ,origin = 'lower')
ax.set_xlabel("x [km]", fontsize=16)
ax.set_ylabel("time [s]", fontsize=16)
plt.xticks(np.arange(0,100,20),np.arange(0,100,20)*25,fontsize=14)
ax.tick_params(labelsize=14)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax).ax.tick_params(labelsize=14)

plt.savefig("ex2_(2).png",bbox_inches="tight",dpi = 300)
plt.show()

## ----------------- ##

N = L/dx
alpha_k = 0
beta_k = 0

C = np.zeros((nt_len,nx_len))
C[0,int(1125/25):int(1375/25+1)] = 10
u = 10

def loop(C):
    for k in range(int(N/2)+1):
        C_star = C[k] - ((2*np.pi*u*k/L)*1j)*C[k]*dt
        C[k] = C[k] - ((2*np.pi*u*k/L)*1j)*C_star*dt
    return C

for nt in range(1,nt_len):
    C_ft = np.fft.rfft(C[nt-1,:])
    C_ft = loop(C_ft)
    C[nt,:] = np.fft.irfft(C_ft) 
    
# ~~~~~~~ #
fig = plt.figure(figsize=(7,6))
ax3 = plt.axes(projection='3d')
ax3.set_xlabel("x [km]", fontsize=14, labelpad=8.5)
ax3.set_ylabel("time [s]", fontsize=14, labelpad=8.5)
ax3.set_zlabel("amplitude", fontsize=14, labelpad=8.5)
ax3.tick_params(labelsize=14)
xx = np.arange(0,L,dx)/1000
yy = np.arange(0,nt_len,1)
X, Y = np.meshgrid(xx, yy)
ax3.plot_surface(X,Y,C,rstride = 2, cstride = 2, cmap='Reds')
plt.show()
# ~~~~~~~ #
plt.figure(figsize=(7,6))
plt.xlabel("x [km]", fontsize=16)
plt.ylabel("time [s]", fontsize=16)
plt.title("Spectral method",fontsize=16)
plt.xticks(np.arange(0,101,20),np.arange(0,101,20)*25,fontsize=16)
plt.yticks(fontsize=16)
plt.imshow(C, cmap="Reds", aspect='auto' ,origin = 'lower')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=16)
plt.show()
# ~~~~~~~ #

fig = plt.figure(figsize=(13,5))
fig.suptitle("Spectral method",fontsize=18)
plt.subplots_adjust(wspace = 0.5)
ax3 = fig.add_subplot(1,2,1,projection='3d')
ax3.set_xlabel("x [km]", fontsize=16, labelpad=8.5)
ax3.set_ylabel("time [s]", fontsize=16, labelpad=8.5)
ax3.set_zlabel("amplitude", fontsize=16, labelpad=8.5)
ax3.tick_params(labelsize=14)
xx = np.arange(0,L,dx)/1000
yy = np.arange(0,nt_len,1)
X, Y = np.meshgrid(xx, yy)
ax3.plot_surface(X,Y,C,rstride=2, cstride=2, cmap='Reds')

ax = fig.add_subplot(1,2,2)
im = ax.imshow(C, cmap="Reds", aspect='auto' ,origin = 'lower')
ax.set_xlabel("x [km]", fontsize=16)
ax.set_ylabel("time [s]", fontsize=16)
plt.xticks(np.arange(0,100,20),np.arange(0,100,20)*25,fontsize=14)
ax.tick_params(labelsize=14)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax).ax.tick_params(labelsize=14)

plt.savefig("ex2_(3).png",bbox_inches="tight",dpi = 300)
plt.show()