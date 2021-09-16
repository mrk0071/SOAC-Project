# -*- coding: utf-8 -*-
"""
@author: Ruoke Meng

Hello Python!

"""

import numpy as np
import matplotlib.pyplot as plt


def Euler(u,v,time,dt):
    KE = [50]
    for nt in range(1,nt_len):
        time[nt] = time[nt-1] + dt
        u[nt] = u[nt-1] + (f*v[nt-1]*dt)
        v[nt] = v[nt-1] - (f*u[nt-1]*dt)
        KE.append((u[nt]**2+v[nt]**2)/2)
    return u, v, KE
    
def leapfrog(u,v,time,dt):
    KE = [50]
    u[1] = u[0] + (f*v[0]*dt)
    v[1] = v[0] - (f*u[0]*dt)
    KE.append((u[1]**2+v[1]**2)/2)
    for nt in range(2,nt_len):
        time[nt] = time[nt-1] + dt
        u[nt] = u[nt-2] + (f*v[nt-1]*2*dt)
        v[nt] = v[nt-2] - (f*u[nt-1]*2*dt)
        KE.append((u[nt]**2+v[nt]**2)/2)
    return u, v, KE

def Matsuno(u,v,time,dt):
    KE = [50]
    for nt in range(1,nt_len):
        time[nt] = time[nt-1] + dt
        u_star = u[nt-1] + f*v[nt-1]*dt
        v_star = v[nt-1] - f*u[nt-1]*dt
        u[nt] = u[nt-1] + (f*v_star*dt)
        v[nt] = v[nt-1] - (f*u_star*dt)
        KE.append((u[nt]**2+v[nt]**2)/2)
    return u, v, KE

def Heun(u,v,time,dt):
    KE = [50]
    for nt in range(1,nt_len):
        time[nt] = time[nt-1] + dt
        u_star = u[nt-1] + f*v[nt-1]*dt
        v_star = v[nt-1] - f*u[nt-1]*dt
        u[nt] = u[nt-1] + (f*v[nt-1]+f*v_star)*dt/2
        v[nt] = v[nt-1] - (f*u[nt-1]+f*u_star)*dt/2
        KE.append((u[nt]**2+v[nt]**2)/2)
    return u, v, KE

## ----- Part I ----- ##
nt_len = 720
dt = 300.0
f = 0.0001

u = np.zeros((nt_len))
v = np.zeros((nt_len))
u[0] = 10.0
u_ana = np.zeros((nt_len))
u_ana[0] = 10.0
time = np.zeros((nt_len))
for nt in range(1,nt_len):
    time[nt] = time[nt-1] + dt
    u_ana[nt] = u[0] * np.cos(f*time[nt])

time = np.zeros((nt_len))   
u = np.zeros((nt_len))
v = np.zeros((nt_len))
u[0] = 10.0 
u_Euler, v_Euler, KE_Euler = Euler(u,v,time,dt)

time = np.zeros((nt_len))   
u = np.zeros((nt_len))
v = np.zeros((nt_len))
u[0] = 10.0 
u_leapfrog, v_leapfrog, KE_leapfrog = leapfrog(u,v,time,dt)

time = np.zeros((nt_len))   
u = np.zeros((nt_len))
v = np.zeros((nt_len))
u[0] = 10.0 
u_Matsuno, v_Matsuno, KE_Matsuno = Matsuno(u,v,time,dt)

time = np.zeros((nt_len))   
u = np.zeros((nt_len))
v = np.zeros((nt_len))
u[0] = 10.0 
u_Heun, v_Heun, KE_Heun = Heun(u,v,time,dt)

plt.figure(figsize=(8,6))
plt.xlabel("time [hrs]", fontsize=15)
plt.ylabel("u [m/s]", fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title("Inertial Oscillation with different computational schemes",fontsize=15)
plt.plot(np.arange(0,dt*nt_len,dt)/3600, u_Euler, label="Euler", linewidth=2)
plt.plot(np.arange(0,dt*nt_len,dt)/3600, u_leapfrog, label="leap-frog", linewidth=2)
plt.plot(np.arange(0,dt*nt_len,dt)/3600, u_Matsuno, label="Matsuno", linewidth=2)
plt.plot(np.arange(0,dt*nt_len,dt)/3600, u_Heun, label="Heun", linewidth=2)
plt.scatter(np.arange(0,dt*nt_len,dt)/3600, u_ana, edgecolors="none",color="black", s=5, label="analytical", zorder=30, alpha=0.7)
plt.legend(scatterpoints=1)
plt.rc('legend', fontsize=14)
plt.grid(linestyle='--')
plt.savefig("ex1_results.png",bbox_inches="tight",dpi = 300)
plt.show()

plt.figure(figsize=(8,6))
plt.xlabel("time [hrs]", fontsize=15)
plt.ylabel("KE [m$^2$s$^{-2}$]", fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title("Kinetic energy with different computational schemes",fontsize=15)
# plt.plot(np.arange(0,dt*nt_len,dt)/3600, KE_Euler, label="Euler", linewidth=2)
plt.plot(np.arange(0,dt*nt_len,dt)/3600, KE_leapfrog, label="leap-frog", linewidth=2)
# plt.plot(np.arange(0,dt*nt_len,dt)/3600, KE_Matsuno, label="Matsuno", linewidth=2)
plt.plot(np.arange(0,dt*nt_len,dt)/3600, KE_Heun, label="Heun", linewidth=2)
plt.legend()
plt.rc('legend', fontsize=14)
plt.grid(linestyle='--')
plt.savefig("ex1_KE_leapfrog&Heun.png",bbox_inches="tight",dpi = 300)
plt.show()

## ----- Part II ----- ##
# dt = [10,50,150,300,600]

# error = np.zeros((4,len(dt)))

# for i in range(len(dt)):
#     nt_len = int(60*3600/dt[i])
    
#     time = np.zeros((nt_len))   
#     u = np.zeros((nt_len))
#     v = np.zeros((nt_len))
#     u[0] = 10.0 
#     u_Euler, v_Euler, KE_Euler = Euler(u,v,time,dt[i])
    
#     time = np.zeros((nt_len))   
#     u = np.zeros((nt_len))
#     v = np.zeros((nt_len))
#     u[0] = 10.0 
#     u_leapfrog, v_leapfrog, KE_leapfrog = leapfrog(u,v,time,dt[i])
    
#     time = np.zeros((nt_len))   
#     u = np.zeros((nt_len))
#     v = np.zeros((nt_len))
#     u[0] = 10.0 
#     u_Matsuno, v_Matsuno, KE_Matsuno = Matsuno(u,v,time,dt[i])
    
#     time = np.zeros((nt_len))   
#     u = np.zeros((nt_len))
#     v = np.zeros((nt_len))
#     u[0] = 10.0 
#     u_Heun, v_Heun, KE_Heun = Heun(u,v,time,dt[i])
    
#     error[0,i] = KE_Euler[-1]/(u_ana[0]**2/2)
#     error[1,i] = KE_leapfrog[-1]/(u_ana[0]**2/2)
#     error[2,i] = KE_Matsuno[-1]/(u_ana[0]**2/2)
#     error[3,i] = KE_Heun[-1]/(u_ana[0]**2/2)


# for i in range(4):
#     regress_1 = np.polyfit(dt, error[i,:], 1)
#     regress_2 = np.polyfit(dt, error[i,:], 2)

#     p_1 = np.poly1d(regress_1)
#     p_2 = np.poly1d(regress_2)
#     xp = np.linspace(1, 650, 650)
#     title = ['Euler', 'leapfrog', 'Matsuno', 'Heun']
#     plt.figure(figsize=(6,5))
#     plt.title(title[i]+' scheme', fontsize=15)
#     plt.xlabel("Time step [s]", fontsize=15)
#     plt.ylabel("Error in KE/KE(t=0) at t=60 hrs", fontsize=15)
#     plt.xticks(fontsize=15)
#     plt.yticks(fontsize=15)
#     plt.scatter(dt, error[i,:], edgecolors="none",color="black", s=30)
#     plt.plot(xp, p_1(xp), linewidth=2)
#     plt.plot(xp, p_2(xp), linewidth=2)
#     plt.legend(['$y = $'+str(format(regress_1[0],'.3f'))+'$x + $'+str(format(regress_1[1],'.3f')), '$y = $'+str(format(regress_2[0],'.3f'))+'$x^2 + $'+str(format(regress_2[1],'.3f'))+'$x + $'+str(format(regress_2[2],'.3f'))])
#     plt.grid(linestyle='--')
#     plt.savefig("ex1_error_"+title[i]+".png",bbox_inches="tight",dpi = 300)
#     plt.show()
