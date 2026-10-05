## Script for the simulation of 2D Complex Ginzburg-Landau Equation (GLE), with an 
# anisotropic diffusion tensor with arbitrary orientation.
# The parameters are set for the systems to develop an chaotic regime, but can be changed 
# to produce other regimes, such as spiral waves.
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

def GLE_ODE(t,uu,d1,d2,phi,N):
    a = 1.1
    b = -1
    L = N#90 # size of domain
    #N = 200  # no. of grid points
    boundary = 'periodic' # 'periodic'/'noflux'
    h=L/N
    hsq=h**2
    #hsq6=6*hsq

    uuu=np.zeros((N+2,N+2),dtype = 'complex_')

    if boundary == 'no-flux':

        # no-flux boundary condition
        
        uuu[0,0]=uu[1,1]
        uuu[0,1:-1]= uu[1,:]
        uuu[0,-1]=uu[1,-2]

        uuu[1:-1,0]=uu[:,1]
        uuu[1:-1,1:-1]=uu
        uuu[1:-1,-1]=uu[:,-2]
        
        uuu[-1,0]=uu[-2,1]
        uuu[-1,1:-1]=uu[-2,:]
        uuu[-1,-1]=uu[-2,-2]

    elif boundary == 'periodic':
        # periodic boundary conditions
        uuu[0,0]=uu[-1,-1]
        uuu[0,1:-1]= uu[-1,:]
        uuu[0,-1]=uu[-1,0]

        uuu[1:-1,0]=uu[:,-1]
        uuu[1:-1,1:-1]=uu
        uuu[1:-1,-1]=uu[:,0]

        uuu[-1,0]=uu[0,-1]
        uuu[-1,1:-1]=uu[0,:]
        uuu[-1,-1]=uu[0,0]

    else:
        raise Exception("Boundary mode not specified or incorrect, mode has to be  'periodic' or 'noflux'")

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
    ## Spatially dependent 

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
    du1 = (1+1j*a)*diff_term + uu - (1+1j*b)*uu*abs(uu)**2

    return np.array(du1)

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

#np.random.seed(1)# r=2
np.random.seed(2)# r=3
#np.random.seed(3)# r=4
#np.random.seed(4)# r=5

N = 200 #sqrt of the number of grid points 
angles = np.linspace(0,90,19) #angles to test (in degrees)
path = '/Volumes/T7/RD_data/GLE/' #output directory
step=.01
simulation_time = 660 # Simulation time, simulation_time/step simulation steps will be 
# performed after the wave break-up.
transitory_time = 160 # Removed transitory time after wave break-up, before saving data.
save_every = 500 #Iteration steps between saved snapshots

#################################################################################
#'spatially_dependent'
for it in range(1):
    id = id_generator()
    for phi_deg in [0]:#angles:

        if phi_deg != 'spatially_dependent':
            phi =  np.radians(phi_deg)
        else: phi = get_angle(N)

        d1 = 2.4
        d2 = .4

        #Initial State
        u = np.random.uniform(-1, 1, (N,N)) + 1.j * np.random.uniform(-1, 1, (N,N))

        u_s = []
        plt.imshow(np.angle(u)); plt.title('Initial conditions')
        plt.set_cmap('hsv'); plt.colorbar()
        plt.show()

        for i in tqdm.tqdm(range(int(simulation_time/step))):


            u_dot = GLE_ODE(1,u,d1,d2,phi,N)
            u = np.array(u+u_dot*step)
            if i > transitory_time/step:
                if i % save_every == 0:
                    #u_s.append(np.angle(u))
                    u_s.append(np.angle(u.flatten()))
                    u_toplot = u
                    #plt.imshow(np.angle(u))
                    #plt.colorbar()
                    #plt.show()


        #plt.imshow(np.angle(u),origin='lower')
        plt.imshow(np.angle(u_toplot),origin='lower')
        plt.colorbar(); plt.show()
        #plt.imshow(phi,origin='lower')
        #plt.colorbar(); plt.show()

        data = np.stack([u_s[i] for i in range(len(u_s))])
        #np.savetxt(path+'AP.csv',data, delimiter=",")
        #np.savetxt(path+'GLE_d1_'+str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi_deg))+'_chaos.csv',data, delimiter=",")
        #np.savetxt(path+'GLE_d1_'+str(d1)+'_d2_'+str(d2)+'_'+id+'.csv',data, delimiter=",")

        '''dataset = xr.Dataset(
            {
                "data": (("t","x"),data),
            }
        )
        filename = path+'GLE_d1_'+str(d1)+'_d2_'+str(d2)+'_phi_'+str(int(phi_deg))+'_chaos.h5'
        dataset.to_netcdf(filename,engine="h5netcdf")'''

        '''d1_mat = d1*np.ones(phi.shape)
        d2_mat = d2*np.ones(phi.shape)

        dataset = xr.Dataset(
            {
                "data": (("t","x","y"),data),
                "angles": (("x","y"),phi),
                "d1": (("x","y"),d1_mat),
                "d2": (("x","y"),d2_mat),
            }
        )
        filename = path+'GLE_d1_'+str(d1)+'_d2_'+str(d2)+'_'+id+'.h5'
        dataset.to_netcdf(filename,engine="h5netcdf")'''



