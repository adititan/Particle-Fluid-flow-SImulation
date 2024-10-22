import matplotlib.pyplot as plt
import numpy as np

sides = 100
radius = 1.0
speed  = 0.0
def axis_function():
	ax = plt.gca()  
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_position(('data',0))
	ax.yaxis.set_ticks_position('left')
	ax.spines['left'].set_position(('data',0))

def polygon(sides, radius, speed, time):
	theta = np.pi * 2.0 / sides
	x = [np.cos(theta*i) * radius + speed.real* time for i in range(sides+1)]
	y = [np.sin(theta*i) * radius +speed.imag* time for i in range(sides+1)]
	normal = [(y[i]-y[i+1]) - 1j*(x[i]-x[i+1])  for i in range(sides)]
	normal = [z/np.absolute(z) for z in normal]
	control_points = [3.0/4.0*(x[i+1] -x[i])+x[i] + 1j*(3.0/4.0*(y[i+1] -y[i])+y[i] ) for i in range(0,sides)]
	return x,y,normal,control_points

def vel_cylinder_exact(z):
	theta = np.arctan2(z.imag,z.real)
	vel = ( vel_fs * ( 1.0 - (radius/abs(z))**2) *np.cos(theta) + 1j *vel_fs * (1 + (radius/abs(z))**2)*np.sin(theta))*np.exp(-1j*theta)
	return vel.conjugate()

def test_vel_cylinder_exact():
	z = 1e-16 +1j*1.0
	expected = 2*vel_fs
	assert (vel_cylinder_exact(z) - expected) < 1e-2

def rk2_exact(delta_t,z_pos):
	for j in range(0,1000):
		u = vel_cylinder_exact(z_pos)
		z_pos_mid = z_pos + 0.5*delta_t * u
		u_mid = vel_cylinder_exact(z_pos_mid)
		z_next = z_pos + delta_t * u_mid
        	plt.plot([z_pos.real,z_next.real],[z_pos.imag,z_next.imag])
		z_pos = z_next
	return z_pos

def vortex_velocity_blob(zpos,z_vort, gamma0):
	delta = 0.0001
	vel = 0.0
	for i in range(0,len(z_vort)):
		vel += gamma0[i]*((z_vort[i].imag - zpos.imag) - 1j*(z_vort[i].real- zpos.real))/(2*np.pi*((abs(zpos - z_vort[i]))**2+ delta**2)) 
	return vel

def test_vortex_blob():
	z = 0.0
	z_vort = [1.0+1j*0.0]
	gamma = [2*np.pi]
	expected = 1j*-1.0
	assert (abs(vortex_velocity_blob(z,z_vort,gamma) - expected)) < 1e-3

def source_velocity(z_pos,z, Q):
	vel = 0.0
	for i in range(0,len(Q)):
		if(np.any(z[i] != z_pos)):
			vel += Q[i]/(2*np.pi*(z_pos - z[i])).conjugate()
	return vel

def test_source_velocity():
	Q = [1,1]
	z = [-1+1j*0,1+1j*0]
	z_pos = 0+1j*0
	expect = 0.0- 1j*0.0
	result = source_velocity(z_pos,z, Q)	
	assert abs(result - expect) < 1e-14

def exact_solution(points,z_vort,A,z_source,strength):
	velocity = [0.0]*len(points)
	for i in range(len(points)):
		velocity[i] += vel_cylinder_exact(points[i])
	z_images = [radius**2/(z.conjugate()) for z in z_vort]
	A_images = [-A[i] for i in range(len(z_vort))]
	z_centre = [0.0 +1j*0.0]*len(z_vort)
	A_centre = [A[i] for i in range(len(z_vort))]
	z_new = z_vort + z_images + z_centre 
	A_new = A + A_images + A_centre
	for i in range(len(points)):
		velocity[i] += vortex_velocity_blob(points[i] ,z_new, A_new)
		velocity[i] += source_velocity(points[i],z_source, strength)
	return velocity

def exact_solution3(delta_t,z_pos,A):		#question3--- path of vortex arounf cylinder, no free stream
	for j in range(0,1000):
		z_images = [radius**2/(z_pos.conjugate())]
		A_images = [-A]
		z_centre = [0.0 +1j*0.0]
		A_centre = [A]
		z_new = [z_pos] + z_images + z_centre 
		A_new = [A] + A_images + A_centre
		u = vortex_velocity_blob(z_pos ,z_new, A_new)
		z_pos_mid = z_pos + 0.5*delta_t * u
		z_new1 = [z_pos_mid] + z_images + z_centre 
		u_mid = vortex_velocity_blob(z_pos_mid ,z_new1, A_new)
		z_next = z_pos + delta_t * u_mid
        	plt.plot([z_pos.real,z_next.real],[z_pos.imag,z_next.imag])
		z_pos = z_next
	return z_pos

def vel_panel(gamma , z , z0, z1):
	theta = np.arctan2((z0.imag - z1.imag),(z0.real-z1.real))
	z_modified = (z-z1)*np.exp(-1j*theta)
	vel1 = -1j*gamma/(2*np.pi)* np.log((z_modified - abs(z0-z1))/z_modified)*np.exp(-1j*theta)
	return vel1.conjugate()

def test_vel_panel():
	gamma =1.0
	z = 0.0+1j*0.0
	z0 , z1 = -1.0 , 1.0
	expected = -0.5
	assert(abs(vel_panel(gamma , z , z0, z1) - expected)) < 1e-3

def flow_across_cylinder(z_vort, A):
	b = [0.0]*(sides+1)
	c = [0.0]*(sides+1)
	velocity = [-vel_fs + speed]*len(control_points)
	for i in range(0,len(control_points)):
		velocity[i] -= source_velocity(control_points[i],z_source, strength)
		velocity[i] -= vortex_velocity_blob(control_points[i],z_vort, A)
		b[i] = np.dot( [velocity[i].real , velocity[i].imag] , [normal[i].real , normal[i].imag] )
	matrix = [[0.0 for w in range(sides)] for h in range(sides+1)]
	for j in range(sides):
		matrix[sides][j] =1.0
	for i in range(0,sides):
		for j in range(0,sides):
			vel = vel_panel(1.0, control_points[i] , x[j+1] + 1j*y[j+1] ,x[j] + 1j*y[j])
			matrix[i][j] = np.dot( [vel.real,vel.imag] , [normal[i].real , normal[i].imag] )
	gamma = np.linalg.lstsq(matrix, b)[0]
	return gamma
	
def question1n2(num_points,distance):
	theta = np.pi * 2 / num_points 
	points = [np.cos(theta*i)*distance + 1j *np.sin(theta*i) * distance for i in range(num_points)]
	velocity = velocity_any_point(x,y,points,gamma,z_vort,A)
	x1 = [points[i].real for i in range(len(points))]
	y1 = [points[i].imag for i in range(len(points))]
	u = [velocity[i].real for i in range(len(velocity))]
	v = [velocity[i].imag for i in range(len(velocity))]
	plt.figure()
	plt.quiver(x1,y1,u,v)		
	axis_function()
	plt.title("Velocity vectors around a cylinder from the approximate method")
	plt.savefig('approx_q1.png')
	velocity2 = exact_solution(points,z_vort,A,z_source,strength)
	error = [0.0]*len(points)
	for i in range(len(points)):	
		error[i] = abs(velocity[i] - velocity2[i])/abs(velocity2[i])
	return (sum(error)/len(error))

def velocity_any_point(x,y,points,gamma,z_vort,gamma0):
	velocity = [vel_fs - speed]*len(points)
	for i in range(0,len(points)):
		velocity[i] += vortex_velocity_blob(points[i],z_vort, gamma0)
		velocity[i] += source_velocity(points[i],z_source, strength)
		for j in range(0,sides):
			velocity[i] += vel_panel(gamma[j] , points[i] , x[j+1] + 1j*y[j+1] ,x[j] + 1j*y[j])
	return velocity

def rk2_panel(delta_t,z_pos,A):
	for j in range(0,1000):
		gamma = flow_across_cylinder([z_pos],[A])
		u = velocity_any_point(x,y,[z_pos],gamma,[z_pos],[A])
		z_pos_mid = z_pos + 0.5*delta_t * u[0]
		u_mid = velocity_any_point(x,y,[z_pos_mid],gamma,[z_pos_mid],[A])
		z_next = z_pos + delta_t * u_mid[0]
        	plt.plot([z_pos.real,z_next.real],[z_pos.imag,z_next.imag])
		z_pos = z_next
	return z_pos

def question3():	
	vel_fs = 0.0
	plt.figure()		
	plt.clf()
	z_pos_nex = rk2_panel(0.005,-1.5+1j*0.0,2*np.pi)
	plt.axis([-2, 2, -2, 2])
	axis_function()
	plt.title("Question3: Using Rk2 to capture movement of vortex around cylinder using Panel Method")
	plt.savefig('q3_panel.png')
	plt.figure()		#q3
	plt.clf()
	z_pos_ex = exact_solution3(0.005,-1.5+1j*0.0,2*np.pi)
	axis_function()
	plt.title("Question3: Using Rk2 to capture movement of vortex around cylinder using Exact Solution")
	plt.savefig('q3_exact.png')
	error = abs(z_pos_nex - z_pos_ex)
	print(error)
	return error

def plot_panel():
	x = np.linspace(-2.5,2.5, 20)  # plot velocity distribution around vortex sheet 
	y = np.linspace(-2.5, 2.5,20)
	X, Y = np.meshgrid(x, y)
	points = [x+1j*y for x in X for y in Y]
	velocity = [vel_panel(1.0 , points[i] , -1.0 + 1j*2.0 , 1.0) for i  in range(len(points))]
	u = [velocity[i].real for i in range(len(velocity))]
	v = [velocity[i].imag for i in range(len(velocity))]
	plt.figure()
	plt.axis([-3, 3, -3, 3])
	plt.plot([-1,1],[2,0],linewidth = 3)
	plt.quiver(X,Y,u,v)
	axis_function()
	plt.title("Velocity vectors around a vortex sheet panel extending from -1,2 to 1,0")
	plt.savefig('panel.png')


if __name__ == '__main__':
	vel_fs = 0.0
	z_vort,A,z_source,strength = [],[],[],[]
	test_vel_cylinder_exact()	#test_functions
	test_source_velocity()
	test_vortex_blob()
	test_vel_panel()
	x,y,normal,control_points = polygon(sides,radius,speed,0.0)	
	question3() 
	plot_panel()
	#ques1&2
	vel_fs = 1.0
	plt.figure()
	gamma  = flow_across_cylinder(z_vort,A)
	distance = [1.5,2.0,2.5,5.0,10.0]		
	errorr = [question1n2(200,distance[i]) for i in range(len(distance))]
	plt.figure()		#plot
	plt.plot(distance,errorr)
	plt.xlabel('Distance')
	plt.ylabel('Error')
	plt.title("Distance VS. Error with 200 panels on the cylinder")
	plt.savefig('result.jpg') 
	plt.figure()		#q1
	plt.clf()
	plt.axis('equal')
	plt.axis([-2, 2, -2, 2])
	rk2_exact(0.01,-1.5+1j*0.00001)
	axis_function()
	plt.title("Question3: Using Rk2 to Capture streamlines around cylinder using Exact Solution")
	plt.savefig('streamline_exact.png')
	plt.figure()
	num_points ,num_points1 ,num_points2 = 20 , 10 , 30		#velocity plot
	distance, distance1, distance2 = 1.5, 1.0, 2.0
	theta, theta1, theta2 = np.pi * 2 / num_points , np.pi * 2 / num_points1 , np.pi * 2 / num_points2
	points1 = [np.cos(theta1*i)*distance1 + 1j *np.sin(theta1*i) * distance1 for i in range(num_points1)]
	points3 = [np.cos(theta*i)*distance + 1j *np.sin(theta*i) * distance for i in range(num_points)]
	points2 = [np.cos(theta2*i)*distance2 + 1j *np.sin(theta2*i) * distance2 for i in range(num_points2)]
	points = points1 + points2 + points3
	z_vort,A,z_source,strength = [],[],[],[]
	velocity = exact_solution(points,z_vort,A,z_source,strength)
	x1 = [points[i].real for i in range(len(points))]
	y1 = [points[i].imag for i in range(len(points))]
	u = [velocity[i].real for i in range(len(velocity))]
	v = [velocity[i].imag for i in range(len(velocity))]
	plt.axis([-3, 3, -3, 3])
	plt.quiver(x1,y1,u,v)
	axis_function()
	plt.title("Question3: Velocity vectors around cylinder using Exact Solution")
	plt.savefig('velocity_plot_exact.png') 000):
		gamma = flow_across_cylinder([z_pos],[A])
		u = velocity_any_point(x,y,[z_pos],gamma,[z_pos],[A])
		z_pos_mid = z_pos + 0.5*delta_t * u[0]
		u_mid = velocity_any_point(x,y,[z_pos_mid],gamma,[z_pos_mid],[A])
		z_next = z_pos + delta_t * u_mid[0]
        	plt.plot([z_pos.real,z_next.real],[z_pos.imag,z_next.imag])
		z_pos = z_next
	return z_pos

def question3():	
	vel_fs = 0.0
	plt.figure()		
	plt.clf()
	z_pos_nex = rk2_panel(0.005,-1.5+1j*0.0,2*np.pi)
	plt.axis([-2, 2, -2, 2])
	axis_function()
	plt.title("Question3: Using Rk2 to capture movement of vortex around cylinder using Panel Method")
	plt.savefig('q3_panel.png')
	plt.figure()		#q3
	plt.clf()
	z_pos_ex = exact_solution3(0.005,-1.5+1j*0.0,2*np.pi)
	axis_function()
	plt.title("Question3: Using Rk2 to capture movement of vortex around cylinder using Exact Solution")
	plt.savefig('q3_exact.png')
	error = abs(z_pos_nex - z_pos_ex)
	print(error)
	return error

def plot_panel():
	x = np.linspace(-2.5,2.5, 20)  # plot velocity distribution around vortex sheet 
	y = np.linspace(-2.5, 2.5,20)
	X, Y = np.meshgrid(x, y)
	points = [x+1j*y for x in X for y in Y]
	velocity = [vel_panel(1.0 , points[i] , -1.0 + 1j*2.0 , 1.0) for i  in range(len(points))]
	u = [velocity[i].real for i in range(len(velocity))]
	v = [velocity[i].imag for i in range(len(velocity))]
	plt.figure()
	plt.axis([-3, 3, -3, 3])
	plt.plot([-1,1],[2,0],linewidth = 3)
	plt.quiver(X,Y,u,v)
	axis_function()
	plt.title("Velocity vectors around a vortex sheet panel extending from -1,2 to 1,0")
	plt.savefig('panel.png')


if __name__ == '__main__':
	vel_fs = 0.0
	z_vort,A,z_source,strength = [],[],[],[]
	test_vel_cylinder_exact()	#test_functions
	test_source_velocity()
	test_vortex_blob()
	test_vel_panel()
	x,y,normal,control_points = polygon(sides,radius,speed,0.0)	
	question3() 
	plot_panel()
	#ques1&2
	vel_fs = 1.0
	plt.figure()
	gamma  = flow_across_cylinder(z_vort,A)
	distance = [1.5,2.0,2.5,5.0,10.0]		
	errorr = [question1n2(200,distance[i]) for i in range(len(distance))]
	plt.figure()		#plot
	plt.plot(distance,errorr)
	plt.xlabel('Distance')
	plt.ylabel('Error')
	plt.title("Distance VS. Error with 200 panels on the cylinder")
	plt.savefig('result.jpg') 
	plt.figure()		#q1
	plt.clf()
	plt.axis('equal')
	plt.axis([-2, 2, -2, 2])
	rk2_exact(0.01,-1.5+1j*0.00001)
	axis_function()
	plt.title("Question3: Using Rk2 to Capture streamlines around cylinder using Exact Solution")
	plt.savefig('streamline_exact.png')
	plt.figure()
	num_points ,num_points1 ,num_points2 = 20 , 10 , 30		#velocity plot
	distance, distance1, distance2 = 1.5, 1.0, 2.0
	theta, theta1, theta2 = np.pi * 2 / num_points , np.pi * 2 / num_points1 , np.pi * 2 / num_points2
	points1 = [np.cos(theta1*i)*distance1 + 1j *np.sin(theta1*i) * distance1 for i in range(num_points1)]
	points3 = [np.cos(theta*i)*distance + 1j *np.sin(theta*i) * distance for i in range(num_points)]
	points2 = [np.cos(theta2*i)*distance2 + 1j *np.sin(theta2*i) * distance2 for i in range(num_points2)]
	points = points1 + points2 + points3
	z_vort,A,z_source,strength = [],[],[],[]
	velocity = exact_solution(points,z_vort,A,z_source,strength)
	x1 = [points[i].real for i in range(len(points))]
	y1 = [points[i].imag for i in range(len(points))]
	u = [velocity[i].real for i in range(len(velocity))]
	v = [velocity[i].imag for i in range(len(velocity))]
	plt.axis([-3, 3, -3, 3])
	plt.quiver(x1,y1,u,v)
	axis_function()
	plt.title("Question3: Velocity vectors around cylinder using Exact Solution")
	plt.savefig('velocity_plot_exact.png')
