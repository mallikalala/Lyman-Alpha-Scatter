import simpy
import random
import math
import time
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
frames = 3
import time

start = time.time()



# graph
graph = plt.figure(figsize=(4,4), dpi=80)
fig = plt.figure()
ax = plt.figure().add_subplot(projection='3d')
import statistics

def simulate(images) :
	#PART 1: Simulating the hydrogen particles
	# Simulation parameters
	Temperature = 270 #Kelvin
	sqrt_T = math.sqrt(Temperature)
	v0H = 111.6769 * sqrt_T #meters per second
	angle = 6.28 # random angle (in radians)
	B= 10 # size of box
	R= 1 # interaction radius
	dt= 0.2 # time step
	t= 30 # number of time steps uj
	num_H_Atoms= 100 # number of H Atoms
	graphRealTime = True

	#PART 2: Simulating the Lyman Alpha particles
	# Simulation parameters
	v0L = 2.0 # velocity
	angle = 6.28 # random angle (in radians)
	B= 10 # size of box
	R= 1 # interaction radius
	num_Ly_Atoms= 50 # number of L Atoms
	graphRealTime = True


	
	""" Finite Volume simulation """
	
	# Initialize
	np.random.seed(num_Ly_Atoms) # set the random number generator 

	# atom positions
	xH = np.random.rand(num_H_Atoms,1)*B
	yH = np.random.rand(num_H_Atoms,1)*B
	zH = np.random.rand(num_H_Atoms,1)*B

	# atom velocities
	theta = 2 * np.pi * np.random.rand(num_H_Atoms,1)
	phi = 2 * np.pi * np.random.rand(num_H_Atoms,1)
	vxH = v0H * np.cos(phi) * np.sin(theta)
	vyH = v0H * np.sin(phi) *np.sin(theta)
	vzH = v0H * np.cos(theta)

	# Simulation Main Loop
	for i in range(t):

		# move
		xH += vxH*dt
		yH += vyH*dt
		zH += vzH*dt 
		# apply periodic BCs
		xH = xH % B
		yH = yH % B
		zH = zH % B 
			
		# find mean direction (angle from 0) of particles within R
		mean_theta = theta
		for b in range(num_H_Atoms):
			particles = (xH-xH[b])**2+(yH-yH[b])**2 < R**2 # Translates directly to the pythagorean theorem a^2 + b^2 = c^2. R is the radius of interaction that the function looks at, which I set to fill the entire box because I want all of the atoms to be able to interact with each other. 
			sxH = np.sum(np.cos(theta[particles]))
			syH = np.sum(np.sin(theta[particles]))
			mean_theta[b] = np.arctan2(syH, sxH) 
			# essentially, this for loop is taking the positions of all of the particles relative to x, the value of all of them and uses a little trigonometry to find the angle of all of them. It then constructs an array of all of the angles in radians. 
		mean_phi = phi
		for b in range(num_H_Atoms):
			particles = (xH-xH[b])**2+(zH-zH[b])**2 < R**2 
			sxH = np.sum(np.cos(phi[particles]))
			szH = np.sum(np.sin(phi[particles]))
			mean_phi[b] = np.arctan2(sxH, szH) 


		# add random directions
		theta = mean_theta + angle*(np.random.rand(num_H_Atoms,1)-0.5)
		phi = mean_phi + angle*(np.random.rand(num_H_Atoms,1)-0.5)

		# update velocities
		vxH = v0H * np.cos(phi) * np.sin(theta)
		vyH = v0H * np.sin(phi) *np.sin(theta)
		vzH = v0H * np.cos(phi)

		#update positions
		xH = xH + vxH*images*dt
		yH = yH + vyH*images*dt
		zH = zH + vzH*images*dt

	# Initialize
	np.random.seed(17) # set the random number generator 

	# photon positions
	xL = np.random.rand(num_H_Atoms,1)*B
	yL = np.random.rand(num_H_Atoms,1)*B
	zL = np.random.rand(num_H_Atoms,1)*B

	# photon velocities
	theta = 2 * np.pi * np.random.rand(num_H_Atoms,1)
	phi = 2 * np.pi * np.random.rand(num_H_Atoms,1)
	vxL = v0L * np.cos(phi) * np.sin(theta)
	vyL = v0L * np.sin(phi) *np.sin(theta)
	vzL = v0L * np.cos(theta)

	# Simulation Main Loop
	for i in range(t):

		# move
		xL += vxL*dt
		yL += vyL*dt
		zL += vzL*dt 
		# apply periodic BCs
		xL = xL % B
		yL = yL % B
		zL = zL % B 
		
		# find mean direction (angle from 0) of particles within R
		mean_theta = theta
		for b in range(num_H_Atoms):
			particles = (xL-xL[b])**2+(yL-yL[b])**2 < R**2 # Translates directly to the pythagorean theorem a^2 + b^2 = c^2. R is the radius of interaction that the function looks at, which I set to fill the entire box because I want all of the atoms to be able to interact with each other. 
			sxL = np.sum(np.cos(theta[particles]))
			syL = np.sum(np.sin(theta[particles]))
			mean_theta[b] = np.arctan2(syL, sxL) 
			# essentially, this for loop is taking the positions of all of the particles relative to x, the value of all of them and uses a little trigonometry to find the angle of all of them. It then constructs an array of all of the angles in radians. 
		mean_phi = phi
		for b in range(num_H_Atoms):
			particles = (xL-xL[b])**2+(zL-zL[b])**2 < R**2 
			sxL = np.sum(np.cos(phi[particles]))
			szL = np.sum(np.sin(phi[particles]))
			mean_phi[b] = np.arctan2(sxL, szL) 


		# add random directions
		theta = mean_theta + angle*(np.random.rand(num_H_Atoms,1)-0.5)
		phi = mean_phi + angle*(np.random.rand(num_H_Atoms,1)-0.5)

		# update velocities
		vxL = v0L * np.cos(phi) * np.sin(theta)
		vyL = v0L * np.sin(phi) *np.sin(theta)
		vzL = v0L * np.cos(phi)

		#update position
		xL = xL + vxL*dt*images
		yL = yL + vyL*dt*images
		zL = zL +vzL*dt*images
	coordinates = []
	while i <num_H_Atoms:
		new_coordinates = [(xH[i],yH[i],zH[i])]
		coordinates.append(new_coordinates)
		i+=1

	"""def binary_search(x):
		temp_coordinates = coordinates
		for i in x:
			length_coord = len(temp_coordinates)
			if i < length_coord:
				length_coord = length_coord/2
				temp_coordinates = temp_coordinates[0:length_coord:1]
			if i > coordinates.length:
				temp_coordinates = temp_coordinates[length_coord/2:length_coord:1]
		return temp_coordinates

	test=[1,2,3,4,5]"""

	


	"""for i in temp_xH:
		print("i'm a friend of coal")
		for j in temp_yH:
			print("no")
			for k in temp_zH:
				print("bye")
				for l in temp_xL:
					print("hu")
					for m in temp_yL:
						print(i)
						for n in temp_zL:
							if (i == l and j == m and k == n):
								vxL = v0L * np.cos(phi) * np.sin(theta)
								vyL = v0L * np.sin(phi) *np.sin(theta)
								vzL = v0L * np.cos(phi)
								print("x")"""
							
	"""rad_object = 5
	x_coord= 2 
	y_coord = -4 
	z_coord = 2
	mass = 10000 #kg (not accurate)
	photons_sec = 10*45 #the sun
	Temp_Object = 30000 # in Kelvin


    #calculated variables
	wavelength = 2.897771955*10**(-3)/Temp_Object #wein's constant divided by temperature yields wavelength
	R_x = rad_object + np.absolute(x_coord)
	R_y = rad_object + np.absolute(y_coord)
	R_z = rad_object + np.absolute(z_coord)
	Radius = np.sqrt(R_x**2 + R_y**2 + R_z**2)
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
	Bz_init= z_coord + rad_object # size of box

	fig = plt.figure()
	num_photons = 500 #scaled down for run time.
	N = num_photons
    # 0 mean, 0.2 std

	xB = np.random.normal(x_coord,Bx_init,N)
	yB = np.random.normal(y_coord,By_init,N)
	zB = np.random.normal(z_coord,Bz_init,N)

	vxB = v0B * np.cos(phi) * np.sin(theta)
	vyB = v0B * np.sin(phi) *np.sin(theta)
	vzB = v0B * np.cos(phi)"""

	# graph in real time
	#calc distance
	#C = np.sqrt((xB-x_coord)**2 + (yB-y_coord)**2 + (zB-z_coord)**2)
	if graphRealTime or (i == t-1):
		plt.cla()
	
		ax.set_xlim3d(0, B)
		ax.set_ylim3d(0, B)
		ax.set_zlim3d(0, B)
		ax.quiver(xH, yH, zH, vxH, vyH, vzH, color='brown',length = 0.5, arrow_length_ratio=0.5, normalize = True)
		ax.quiver(xL, yL, zL, vxL, vyL, vzL, color='red',length = 0.5, arrow_length_ratio=0.5, normalize = True)
		#ax.quiver(xB, yB, zB, vxB, vyB, vzB, cmap='autumn', length = 0.8, arrow_length_ratio=0.8, normalize = True)
		ax.set(xlim=(0, B), ylim=(0, B))
		ax.set_aspect('equal')	
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)
		ax.get_zaxis().set_visible(False)
		plt.pause(0.001)
					
		# Save figure			
		
		
	return fig
i = 1
while i<frames:
	name = "3dsimulationframe" + str(i) + ".png"
	print(simulate(i))
	plt.savefig(name,dpi=240)
	plt.show()
	i+=1

end = time.time()
print(end - start)