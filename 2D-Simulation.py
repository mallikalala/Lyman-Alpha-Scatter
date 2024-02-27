

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from random import randint, uniform

# graph
frames = 20
graph = plt.figure(figsize=(4,4), dpi=80)
ax = plt.gca()

def simulate (images):
	#PART 1: Simulating the hydrogen particles
	# Simulation parameters
	v0H = 2.0 # velocity
	angle = 6.28 # (in radians)
	B= 5 # size of box
	R= 1 # interaction radius
	dt= 0.1 # time step, in seconds
	t= 20 # number of time steps
	num_H_Atoms= 100 # number of H Atoms
	graphRealTime = True
	images = images

	#PART 2: Simulating the Lyman Alpha particles
	# Simulation parameters
	v0L = 2.0 # velocity
	angle = 6.28 # random angle (in radians)
	B= 5 # size of box
	R= 1 # interaction radius
	num_Ly_Atoms= 50 # number of H Atoms
	graphRealTime = True



	#Hydrogen atom positions
	xH = np.random.rand(num_H_Atoms,1)*B 
	yH = np.random.rand(num_H_Atoms,1)*B

	# Initialize
	np.random.seed(17) # set the random number generator seed



	# atom velocities
	theta = 2 * np.pi * np.random.rand(num_H_Atoms,1)
	vxH = v0H * np.cos(theta)
	vyH = v0H * np.sin(theta)
		

	# Simulation Main Loop
	for i in range(t):

		# move
		xH += vxH*dt
		yH += vyH*dt
		# apply periodic BCs
		xH = xH % B
		yH = yH % B
			
		# find mean angle of neighbors within R
		mean_theta = theta
		for b in range(num_H_Atoms):
			neighbors = (xH-xH[b])**2+(yH-yH[b])**2 < R**2
			sxH = np.sum(np.cos(theta[neighbors]))
			syH = np.sum(np.sin(theta[neighbors]))
			mean_theta[b] = np.arctan2(syH, sxH)
				
			# add random perturbations
			theta = mean_theta + angle*(np.random.rand(num_H_Atoms,1)-0.5)
			
			# update velocities
			vxH = v0H * np.cos(theta)
			vyH = v0H * np.sin(theta)

			#update positions
			xH = xH + vxH*images*dt
			yH = yH + vyH*images*dt
		

	# photons positions	
	xL = np.random.rand(num_Ly_Atoms,1)*B
	yL = np.random.rand(num_Ly_Atoms,1)*B

	# photons velocities
	theta = 2 * np.pi * np.random.rand(num_Ly_Atoms,1)
	vxL = v0L * np.cos(theta)
	vyL = v0L * np.sin(theta)
		
	# Simulation Main Loop LYMAN
	for i in range(t):

		# move
		xL += vxL*dt
		yL += vyL*dt
			
		# apply periodic BCs
		xL = xL % B
		yL = yL % B
			
		# find mean angle of neighbors within R
		mean_theta = theta
		for b in range(num_Ly_Atoms):
			neighbors = (xL-xL[b])**2+(yL-yL[b])**2 < R**2
			sxL = np.sum(np.cos(theta[neighbors]))
			syL = np.sum(np.sin(theta[neighbors]))
			mean_theta[b] = np.arctan2(syL, sxL)
				
		# add direction
		theta = mean_theta + angle*(np.random.rand(num_Ly_Atoms,1)-0.5)
			
		# update velocities
		vxL = v0L * np.cos(theta)
		vyL = v0L * np.sin(theta)

	for i in xH:
		for j in yH:
			for k in xL:
				for l in yH:
					if (i == k and j == l):
						vxL = v0L * np.cos(theta)
						vyL = v0L * np.sin(theta)
					#update positions
	xL = xL + vxL*dt*images
	yL = yL + vyL*dt*images
	#input variables

	rad_object = 5
	x_coord= 0 
	y_coord = 0 
	mass = 10000 #kg (not accurate)
	photons_sec = 10*45 #the sun
	Temp_Object = 30000 # in Kelvin


    #calculated variables
	wavelength = 2.897771955*10**(-3)/Temp_Object #wein's constant divided by temperature yields wavelength
	R_x = rad_object + np.absolute(x_coord)
	R_y = rad_object + np.absolute(y_coord)
	Radius = np.sqrt(R_x**2 + R_y**2)
	watts = mass*Temp_Object*299792458 #speed of light
	intensity = watts/((4*np.pi*rad_object**2)*R**2) #watts/meters^2
	photons_dt = photons_sec*dt #photons per meter^2 per time unit
	#print(photons_dt)
	#print (intensity)

    #Simulating the background photons
    # Simulation parameters
	v0B = np.sqrt(watts/0.5)
	angle = np.pi/6 # angle of dispersion in radians)
	Bx_init= x_coord + rad_object # size of box
	By_init= y_coord + rad_object # size of box
	fig = plt.figure()
	num_photons = 500 #scaled down for run time.
	N = num_photons
    # 0 mean, 0.2 std

	xB = np.random.normal(x_coord,Bx_init,N)
	yB = np.random.normal(y_coord,By_init,N)
	theta = 2 * np.pi * np.random.rand(num_photons,1)
	vxB = v0B * np.cos(theta)
	vyB = v0B * np.sin(theta)
	print(vyB)

    # calculate the distance to (0, 0).
	C = np.sqrt((xB-x_coord)**2 + (yB-y_coord)**2)

	fig = plt.figure()
	# graph in real time
	if graphRealTime:
		plt.cla()
		plt.quiver(xH,yH,vxH,vyH, color = "brown")
		plt.quiver(xL,yL,vxL,vyL, color = "red")
		plt.quiver(xB , yB, vxB,vyB, C,cmap='autumn', alpha = 0.3)
		ax.set(xlim=(0, B), ylim=(0, B))
		ax.set_aspect('equal')	
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)
		plt.pause(0.001)
		print("woah")
	
	return fig

i = 1
while i<frames:
	name = "2dSimulation9" + str(i) + ".png"
	print(simulate(i))
	plt.savefig(name,dpi=240)
	plt.show()
	i+=1






plt.show()