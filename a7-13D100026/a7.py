import matplotlib.pyplot as plt
import numpy as np
import random
from copy import deepcopy

time = 3.0
sides = 30

def plot_vortices(z_vort,A):
	x_nex_pos,y_nex_pos,x_nex_neg,y_nex_neg = [],[],[],[]	
	for a in range(0,len(z_vort)):	
		if (A[a] > 0.0):
			x_nex_pos.append(z_vort[a].real)
			y_nex_pos.append(z_vort[a].imag)
		elif (A[a] < 0.0):
			x_nex_neg.append(z_vort[a].real)
			y_nex_neg.append(z_vort[a].imag)
	plt.figure()
	xc = np.linspace(-1,1,100)
	yc = np.sqrt(1.0 - xc*xc)
	plt.plot(xc,yc,color= 'b')
	yc = -np.sqrt(1.0 - xc*xc)
	plt.plot(xc,yc,color = 'b')
	plt.scatter(x_nex_pos,y_nex_pos,color='r')
	plt.scatter(x_nex_neg,y_nex_neg,color='b')
	plt.axis([-2, 2, -2, 2])
	axis_function()
	plt.title("Vorticity Distribution around Cylinder after 3 seconds")
	plt.savefig('three.png')

def axis_function():
	ax = plt.gca()  
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_position(('data',0))
	ax.yaxis.set_ticks_position('left')
	ax.spines['left'].set_position(('data',0))

def polygon(sides, radius):
	theta = np.pi * 2.0 / sides
	x = [np.cos(theta*i) * radius for i in range(sides+1)]
	y = [np.sin(theta*i) * radius for i in range(sides+1)]
	normal = [(y[i]-y[i+1]) - 1j*(x[i]-x[i+1])  for i in range(sides)]
	normal = [z/np.absolute(z) for z in normal]
	control_points = [1.0/2.0*(x[i+1]+x[i]) + 1.0j/2.0*(y[i+1]+y[i] ) for i in range(0,sides)]
	return x,y,normal,control_points

def vortex_velocity_blob(zpos,z_vort, gamma0):
	delta = 0.1
	vel = 0.0
	for i in range(0,len(z_vort)):
		vel += gamma0[i]*((z_vort[i].imag - zpos.imag) - 1j*(z_vort[i].real- zpos.real))/(2*np.pi*((abs(zpos - z_vort[i]))**2+ delta**2)) 
	return vel

def vel_panel(gamma1,gamma2 , z , z0, z1):
	theta = np.angle(z0-z1)
	panel_length = abs(z0-z1)
	z_modified = (z-z1)*np.exp(-1j*theta)
	vel = (-1j/(2*np.pi)*((gamma1-gamma2)*(z_modified/panel_length*np.log((z_modified - panel_length)/z_modified) + 1.0) - gamma1*np.log((z_modified - panel_length)/z_modified))).conjugate()
	return vel*np.exp(1j*theta)

def computeA():
	matrix = [[0.0 for w in range(sides)] for h in range(sides+1)]
	for j in range(sides):
		matrix[sides][j] =1.0
	for i in range(0,sides):
		for j in range(0,sides):
			vel = vel_panel(1.0, 0.0, control_points[i] , x[j+1] + 1j*y[j+1] ,x[j] + 1j*y[j])
			vel2 = vel_panel(0.0, 1.0, control_points[i] , x[j+1] + 1j*y[j+1] ,x[j] + 1j*y[j])
			matrix[i][j] += np.dot( [vel.real,vel.imag] , [normal[i].real , normal[i].imag] )
			if ( j == sides-1) :			
				matrix[i][0] += np.dot( [vel2.real,vel2.imag] , [normal[i].real , normal[i].imag] )
			else:
				matrix[i][j+1] += np.dot( [vel2.real,vel2.imag] , [normal[i].real , normal[i].imag] )
	return matrix

def panel_method(z_vort, A,matrix):
	b = [0.0+1.0j*0.0]*(sides+1)
	velocity = [-vel_fs]*sides
	for i in range(0,sides):
		velocity[i] -= vortex_velocity_blob(control_points[i],z_vort, A)
		b[i] = np.dot( [velocity[i].real , velocity[i].imag] , [normal[i].real , normal[i].imag] )
	b[sides] = 0.0
	gamma = np.linalg.lstsq(matrix, b)[0]
	return gamma

def velocity_any_point(x,y,points,gamma,z_vort,gamma0):
	velocity = [vel_fs]*len(points)
	for i in range(0,len(points)):
		velocity[i] += vortex_velocity_blob(points[i],z_vort, gamma0)
		for j in range(0,sides):
			if( j == sides-1 ):
				velocity[i] += vel_panel(gamma[j], gamma[0] , points[i] , x[j+1] + 1j*y[j+1] ,x[j] + 1j*y[j])
			else:
				velocity[i] += vel_panel(gamma[j], gamma[j+1] , points[i] , x[j+1] + 1j*y[j+1] ,x[j] + 1j*y[j])
	return velocity

def plot_streamlines(gamma,z_vort,gamma0):
	xa = np.linspace(0, 2, 20)
	ya = np.linspace(0, 2, 20)
	xv, yv = np.meshgrid(xa, ya)
	z_point = [xi + 1j*xj for xi in xv for xj in yv]
	veloc = velocity_any_point(x,y,z_point,gamma,z_vort,gamma0)
	u = [veloc[i].real for i in range(len(z_point))]
	v = [veloc[i].imag for i in range(len(z_point))]
	plt.figure()
	xc = np.linspace(-1,1,100)
	yc = np.sqrt(1.0 - xc*xc)
	plt.plot(xc,yc,color= 'b')
	yc = -np.sqrt(1.0 - xc*xc)
	plt.plot(xc,yc,color = 'b')
	plt.quiver(xv, yv, u, v )
	plt.axis([0, 2, 0, 2])
	plt.title("Velocity vector plot around cylinder for a Viscous flow")
	plt.savefig('streamlines_three.png')

def diffuse(z_vor):
	z_vortices = deepcopy(z_vor)
	z_vortices += np.random.normal(mean, sigma, 1)[0] + 1j*np.random.normal(mean, sigma, 1)[0]		
	return z_vortices

def rvm_new_blobs(cp,num,gamma_0,num_of_vort):
	if (num_of_vort == 0):
		return [],[]
	z_new = [cp]*num_of_vort
	gamma_new = [gamma_0/num_of_vort*panel_length]*num_of_vort
	for k in range(0,len(z_new)):		
		z_new[k] = diffuse(z_new[k])
		z_new[k] = reflect2(z_new[k],num)
	return z_new,gamma_new

def reflect2(z,rr):
	def line_eq(xx,yy,r):
		f = y[r] - yy + ((y[r+1] - y[r])/(x[r+1]-x[r]))*(xx-x[r])
		return f	
	if (line_eq(0.0,0.0,rr)*line_eq(z.real,z.imag,rr) > 0):
		theta = np.angle((x[rr+1]-x[rr])+1.0j*(y[rr+1]-y[rr]))
		z_reflect = (z - control_points[rr])*np.exp(-1.0j*theta)
		z_reflected = (z_reflect.conjugate())*np.exp(1.0j*theta)+control_points[rr]
	else:
		z_reflected = z
	return z_reflected

def reflect(z):
	theta = np.arctan2(z.imag,z.real)
	z += 2*(1.0 - abs(z))*(np.cos(theta) + 1j*np.sin(theta))
	return z

def rk2_advect(z_vor,matrix,gamma,A):
	z_vore = [0.0]*len(z_vor)
	for k in range(0,len(z_vor)):
		z_pos_mid = deepcopy(z_vor)
		u = velocity_any_point(x,y,[z_vor[k]],gamma,z_vor,A)		
		z_pos_mid[k] = z_vor[k] + 0.5*delta_t * u[0]
		if(radius - abs(z_pos_mid[k]) > 1e-3):
			z_pos_mid[k] = reflect(z_pos_mid[k])
		gamma_new = panel_method(z_pos_mid,A,matrix)
		u_mid = velocity_any_point(x,y,[z_pos_mid[k]],gamma_new,z_pos_mid,A)		
		z_vore[k] = z_vor[k] + delta_t * u_mid[0]
		if(radius - abs(z_vor[k]) > 1e-3):	
			z_vore[k] = reflect(z_vore[k])
	return z_vore

def vortex_momentum(z_vort,A):
	F= 0
	for d in range(0,len(z_vort)):
		F += (z_vort[d].imag -1j*z_vort[d].real)*A[d]
	return F

def viscous_flow_across_cylinder():
	matrix = computeA()
	t = 0.0
	gamma_max =0.1
	z_vort = []
	gamma0 = []
	I = []
	while(t < time):
		print(len(z_vort),t)
		t+= delta_t
		gamma = panel_method(z_vort,gamma0,matrix)
		z_vort = rk2_advect(z_vort,matrix,gamma,gamma0)
		for rv in range(0,len(gamma)):
			if (rv == len(gamma)-1):
				gamma_rvm = (gamma[rv] + gamma[0])/2.0
			else:
				gamma_rvm = (gamma[rv] + gamma[rv+1])/2.0
			num_of_vort = int(round(abs(gamma_rvm/gamma_max*panel_length)))
			z_new,gamma_new = rvm_new_blobs(control_points[rv],rv,gamma_rvm,num_of_vort)
			z_vort += z_new	
			gamma0 += gamma_new
		I.append(vortex_momentum(z_vort,gamma0)) 
	gamma = panel_method(z_vort,gamma0,matrix)
	plot_vortices(z_vort,gamma0)
	plot_streamlines(gamma,z_vort,gamma0)
	cd = [(I[h+1].real - I[h].real)/delta_t for h in range(0,len(I)-1)]
	cd_moving = [0.0]*(len(cd)-2)
	for c in range(0,len(cd)-2):
		cd_moving[c] = cd[c] + cd[c+1] + cd[c+2]
	tt = [i*delta_t for i in range(1,int(round(time/delta_t))-1)]
	plt.figure()
	plt.title("Coefficient of drag Vs. Time")
	plt.xlabel('Time')
	plt.ylabel('C_d')
	plt.plot(tt,cd_moving)
	plt.savefig('cd_V_time.png')

if __name__ == '__main__':
	vel_fs = 1.0
	radius = 1.0
	x,y,normal,control_points = polygon(sides, radius)
	panel_length = 2*radius*np.sin(np.pi/sides)
	delta_t = 0.1
	Re = 1000
	mu = 2*radius*vel_fs/Re #viscocity
	mean = 0.0	
	sigma = np.sqrt(2*mu*delta_t)	
	viscous_flow_across_cylinder()
