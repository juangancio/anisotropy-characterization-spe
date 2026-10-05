## Script for the simulation of Aliev-Panfilov 2D reaction-diffusion system, with an 
# anisotropic diffusion tensor with arbitrary orientation.
# The parameters are set for the systems to develop an chaotic regime, but can be changed 
# to produce other regimes, such as spiral waves.
#
# The script has an adaptative wave break-up process that should start the chaotic regime
# for any anisotropy. But for high anisotropy ratios it could fail to initialize or stop 
# too soon. In that case, there are two parameters that can be tweaked in order to obtain
# better results.
#  
# The simulation uses '9-point stencil' for the diffusion term, but it reduces to a 
# 5-point stencil when the diffusion axes align with the cartesian axes. This implies that
# the stencial is actually 5-point, which can cause some problems, such as an intricsic numerical
# anisotropy. This should be fixed by generalizing the 9-point stencil for an arbitrary orientation.

import xarray as xr
import numpy as np
import tqdm, string, random
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline


def AlievPanfilov_ODE(t,u,d1,d2,phi,N):

    # Conditions implemented back in Göttingen
    k = 8
    a = 0.02
    m1 = 0.09
    m2 = 0.3
    e = 0.01
    L = N # size of domain
    #N = 200  # no. of grid points

    '''#Conditions as Ulrich paper (Entropy, 2015)
    k = 8
    a = 0.05
    m1 = 0.2
    m2 = 0.3
    e = 0.01
    L = N # size of domain'''




    h=L/N
    hsq=h**2
    #hsq6=6*hsq
    uu = u[0,:,:]
    vv = u[1,:,:]

    # no-flux boundary condition
    uuu=np.zeros((N+2,N+2))
    
    uuu[0,0]=uu[1,1]
    uuu[0,1:-1]= uu[1,:]
    uuu[0,-1]=uu[1,-2]

    uuu[1:-1,0]=uu[:,1]
    uuu[1:-1,1:-1]=uu
    uuu[1:-1,-1]=uu[:,-2]
    
    uuu[-1,0]=uu[-2,1]
    uuu[-1,1:-1]=uu[-2,:]
    uuu[-1,-1]=uu[-2,-2]

    '''
    #9 point stencil 
    diff_term = d1 * ( 
                4 * (uuu[1:-1,:-2] + uuu[1:-1,2:] + 
                        uuu[2:,1:-1] + uuu[:-2,1:-1] ) + 
                        uuu[2:,2:] + uuu[2:,:-2] + uuu[:-2,2:] + 
                        uuu[:-2,:-2] - 20 *uu  )  / hsq6  '''
    
    '''#5 point stencil isotropic
    diff_term=d*(uuu[1:-1,:-2] + uuu[1:-1,2:] + 
                        uuu[2:,1:-1] + uuu[:-2,1:-1]  - 4 *uu  )  / hsq
    
    #5 point stencil anisotropic (axis aligned)
    diff_term=(d1 * uuu[1:-1,:-2] + d1 * uuu[1:-1,2:] + 
                       d2 * uuu[2:,1:-1] + d2 * uuu[:-2,1:-1]  - d1 * 2 *uu - d2 * 2 * uu  )  / hsq
    '''

    #5 point stecil arbitrary direction
    a11 = d1*np.cos(phi)**2+d2*np.sin(phi)**2
    a22 = d2*np.cos(phi)**2+d1*np.sin(phi)**2

    a12 = (d1-d2)*np.cos(phi)*np.sin(phi)

    diff_term = (2 * a11 * uuu[1:-1,:-2] + 2 * a11 * uuu[1:-1,2:] + 
                        2 * a22 * uuu[2:,1:-1] + 2 * a22 * uuu[:-2,1:-1] 
                        + a12 * uuu[:-2,:-2]  + a12 * uuu[2:,2:]
                        - a12 * uuu[2:,:-2]  - a12 * uuu[:-2,2:]
                          - a11 * 4 *uu - a22 * 4 * uu  )  / ( 2  * hsq )


    # finite difference formula
    du1 =  k*uu*(uu-a)*(1-uu)- vv*uu + diff_term  
    du2=  (e + (m1*vv)/(m2+uu))*(-vv - k*uu*(uu-a-1))

    return np.array([du1,du2])

def get_angle(N):
    x = np.arange(0, N+60, 60)    
    y = np.arange(0, N+60, 60)  
    z = np.random.uniform(-.1,.1,(x.size,y.size))
    func = RectBivariateSpline(x, y, z, s=0)

    xnew = np.linspace(0, N, N)
    ynew = np.linspace(0, N, N)
    znew = func(xnew, ynew)

    znew = 90*(znew-min(znew.flatten()))/(max(znew.flatten())-min(znew.flatten()))
    return np.radians(znew)


def id_generator(size=4, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

#################################################################################
## Set parameters

#np.random.seed(9)# r=2
#np.random.seed(10)# r=3
#np.random.seed(11)# r=4
np.random.seed(12)# r=5

N = 200 #sqrt of the number of grid points 
angles = np.linspace(0,90,19) #angles to test (in degrees)
path = '/Volumes/T7/RD_data/AP/' #output directory
step=.1
simulation_time = 6500 # Simulation time, simulation_time/step simulation steps will be 
# performed after the wave break-up.
transitory_time = 1500 # Removed transitory time after wave break-up, before saving data.
save_every = 500 #Iteration steps between saved snapshots


## Wave breaking parameters:

# Starting a chaotic regime that chaotic regime that persists long enough to be analyzed
# can be challenging. The following two parameters control the wave breakup that initializes 
# this behaviour, and can be slighted tweaked if the chaotic regime is not being reached or
# is not self-sustaining for long enough.

planar_wave_direction = 'vertical' # 'horizontal' or 'vertical'. 
# In high anisotropy ratio cases the initial propagation direction of the initial planar
# wave can make a difference. Try chage direction if the chaotic regime is not starting.

delay_breakup_param = .65 #between 0 and 1
# This parameter determines when the wave break-up is started. It sets the position of a 
# trigger in the 'tail' of the planar wave, that its activated when the wave drops below .5
# This trigger is only activated if the front of the planar wave has passed through the 
# middle of the grid.
# 0 corresponds to the border, and 1 to the middle of the grid.
# Try between .7 and .85 (values outside this interval are untested).
# Values between 1 and 2 correspond to the front of the wave, and probably would make the 
# wave breakup triggers as soon as it has reached the middle fo the grid (code would need
# to be modified if the wave break-up wants to be started later).
# Values outside (0,2) should raise an error since correspond to points outside the grid.

#################################################################################
for it in range(1):
    id = id_generator()
    for phi_deg in np.linspace(0,90,19):

        if phi_deg != 'spatially_dependent':
            phi =   np.radians(phi_deg)
        else: phi = get_angle(N)

        d1 = 4.8
        d2 = .8

        #Initial State
        us = 0.0  #  u steady state
        vs = 0.0   # v steady state
        uu = np.array([us] * N**2).reshape(N,N)
        vv = np.array([vs] * N**2).reshape(N,N)

        #uu[:, 1:2*N//3] = .9
        #vv[:2*N//3 - N//12, :] = .9
        if planar_wave_direction == 'horizontal':
            uu[:,:int(N/20)] = .9
        elif planar_wave_direction == 'vertical':
            uu[:int(N/20),:] = .9   
        else:
            raise Exception("Specified planar wave propagation direction is not valid, try 'vertical' or 'horizontal'.")

        u0 = np.array([uu,vv])
        t=1
        u = u0
        u_s = []

        first_wave = 0
        trigger = True
        #for i in tqdm.tqdm(range(int(110/step))):
        while trigger:
            u_dot = AlievPanfilov_ODE(t,u,d1,d2,phi,N)
            u = np.array([u[0]+u_dot[0]*step, u[1]+u_dot[1]*step])
            if first_wave == 1:
                if u[0,int(delay_breakup_param*N/2),int(delay_breakup_param*N/2)] < .5:
                    trigger = False
            elif u[0,int(N/2),int(N/2)] > .5:
                first_wave = 1
            

            #if i%save_every==0:
            #        u_s.append(u[0].flatten())

        plt.imshow(u[0])
        plt.colorbar()
        plt.show()

        uu = u[0]
        vv = u[1]

        if planar_wave_direction == 'horizontal':
            uu[:int(N/2),:] = 1
        elif planar_wave_direction == 'vertical':
            uu[:,:int(N/2)] = 1
        
        u0 = np.array([uu,vv])
        u = u0

        for i in tqdm.tqdm(range(int(simulation_time/step))):


            u_dot = AlievPanfilov_ODE(t,u,d1,d2,phi,N)
            u = np.array([u[0]+u_dot[0]*step, u[1]+u_dot[1]*step])
            if i > transitory_time/step:
                if i % save_every == 0:
                    #u_s.append(u[0])
                    u_s.append(u[0].flatten())
                    #u_toplot = u[0]

        plt.imshow(u[0], origin='lower')
        plt.colorbar(); plt.show()
        #plt.imshow(phi,origin='lower')
        #plt.colorbar(); plt.show()

        data = np.stack([u_s[i] for i in range(len(u_s))])
        #np.savetxt(path+'AP.csv',data, delimiter=",")
        #np.savetxt(path+'AP_d1_'+str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi_deg))+'_ulrich.csv',data, delimiter=",")
        #np.savetxt(path+'AP_d1_'+str(d1)+'_d2_'+str(d2)+'_spatially_dependent_angle_ulrich.csv',data, delimiter=",")

        d1_mat = d1*np.ones(phi.shape)
        d2_mat = d2*np.ones(phi.shape)

        dataset = xr.Dataset(
            {
                "data": (("t","x"),data),
            }
        )
        filename = path+'AP_d1_'+str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi_deg))+'_ulrich.h5'
        dataset.to_netcdf(filename,engine="h5netcdf")

        '''dataset = xr.Dataset(
            {
                "data": (("t","x","y"),data),
                "angles": (("x","y"),phi),
                "d1": (("x","y"),d1_mat),
                "d2": (("x","y"),d2_mat),
            }
        )
        filename = path+'AP_d1_'+str(d1)+'_d2_'+str(d2)+'_'+id+'.h5'
        dataset.to_netcdf(filename,engine="h5netcdf")'''


