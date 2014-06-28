import numpy as np

OCEAN = dict(salt=[0., 38.], temp=[0., 31.], ubar=[-0.6, 0.6], vbar=[-0.6, 0.6], 
	         zeta=[-1., 1.], u=[-0.6, 0.6], v=[-0.6, 0.6])

TIDE = dict(tide_Cangle=[0, 2*np.pi], tide_Cmax=[0, 2], tide_Cmin=[0, 2],
	        tide_Cphase=[0, 2*np.pi], tide_Eamp=[0, 5], tide_Ephase=[0, 2*np.pi])

METEO = dict(apratesfc=[0, 0.05], dlwrfsfc=[0, 500], dswrfsfc=[0, 1000], 
	         mslp=[99200, 103500], rh2m=[0, 100], tmp2m=[-30, 40], tmpsfc=[243, 313], 
	         ugrd10m=[-15, 15], vgrd10m=[-15, 15])
