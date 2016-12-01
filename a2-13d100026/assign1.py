import webbrowser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def phi(z_pos,U,Q,z):
	phi = U*z_pos
	for i in range(0,len(Q)):
		if(np.any(z[i] != z_pos)):
			phi += 1/(2 * np.pi)*Q[i]*np.log(z_pos - z[i])
	return phi

def test_phi():
	Q = [1,-1]
	z = [-1+1j*0,1+1j*0]
	z_pos = 0+1j*0
	expect = 0 - 1j*0.5
	result = phi(z_pos,U,Q,z)
	assert abs(result.imag - expect.imag) < 1e-14
	assert abs(result.real - expect.real) < 1e-14

def axis_function():
	ax = plt.gca()  
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_position(('data',0))
	ax.yaxis.set_ticks_position('left')
	ax.spines['left'].set_position(('data',0))


def problem1():
	plt.figure(1)
	plt.clf()
	plt.axes([0.025, 0.025, 0.95, 0.95])
	plt.contourf(X, Y,  phi(X+ 1j*Y,U,Q,z).real , 40, alpha=.75, cmap=plt.cm.hot)
	C = plt.contour(X, Y, phi(X+ 1j*Y,U,Q,z).real, 40, colors='black', linewidth=.5)
	plt.contourf(X, Y, phi(X+ 1j*Y,U,Q,z).imag, 40, alpha=.75, cmap=plt.cm.hot)
	S = plt.contour(X, Y,  phi(X+ 1j*Y,U,Q,z).imag, 40, colors='blue', linewidth=.5)
	axis_function()
	plt.xticks([-1,0, 1],[r'$-1$', r'$0$', r'$1$'])
	plt.yticks([-2, 2],[r'$-2$', r'$2$'])
	plt.savefig('one.png') 

####Ques2
def source_velocity(z_pos,z, Q,U):
	vel = U
	for i in range(0,len(Q)):
		if(np.any(z[i] != z_pos)):
			vel += Q[i]/(2*np.pi*(z_pos - z[i])).conjugate()
	return vel.real,vel.imag

def test_source_velocity():
	Q = [1,-1]
	z = [-1+1j*0,1+1j*0]
	z_pos = 0+1j*0
	expect = 1.3183098861837905 - 1j*0
	[result_real,result_imag] = source_velocity(z_pos,z, Q, U)	
	assert abs(result_imag - expect.imag) < 1e-14
	assert abs(result_real - expect.real) < 1e-14

def euler_ques2(delta_t):
	ques2 = np.linspace(-2,2,10)
	for j in range(0,10):
		x_0 , y_0 = -2 , ques2[j]
		for i in range(0,500):
			[u,v] = source_velocity(x_0+ 1j*y_0,z,Q,U)
			x_1 = x_0 + delta_t * u
			y_1 = y_0 + delta_t * v
			plt.plot([x_0,x_1],[y_0,y_1])
			x_0 , y_0 = x_1 , y_1
	return x_0,y_0

def runge_ques2(delta_t):
	for j in range(0,10):
		ques2 = np.linspace(-2,2,10)
		x_0 , y_0 = -2 , ques2[j]
		for i in range(0,500):
			[u,v] = source_velocity(x_0+ 1j*y_0,z,Q,U)
			x_mid = x_0 + delta_t * u
			y_mid = y_0 + delta_t * v
			[u1,v1] = source_velocity(x_mid + 1j*y_mid,z,Q,U)
			x_1 = x_0 + 0.5*delta_t * u + 0.5*delta_t * u1
			y_1 = y_0 + 0.5*delta_t * v + 0.5*delta_t * v1
			plt.plot([x_0,x_1],[y_0,y_1])
			x_0 , y_0 = x_1 , y_1
	return x_0,y_0

def problem2():
	plt.figure()
	plt.clf()
	euler_ques2(0.01)
	axis_function()
	plt.title("Question2: Using Euler Integration for t=0.01")
	plt.savefig('two.png')	
	plt.figure()
	plt.clf()
	runge_ques2(0.01)
	axis_function()
	plt.title("Question2: Using Runge-Kutta Integration for t=0.01")
	plt.savefig('three.png')

####Ques3
def vortex_velocity(zpos,z_vort, A):
	vel = 0
	for i in range(0,len(A)):
		if(np.any(z_vort[i] != zpos)):
			vel += 1j*A[i]/(2*np.pi*(zpos - z_vort[i])).conjugate()
	return vel

def test_vortex_velocity():
	A=[2*np.pi,2*np.pi]
	z_vor=[0+1j*0 , 1+1j*0] 
	z_pos = 0+1j*0
	expect = 0.0-1j*1.0
	result = vortex_velocity(z_pos,z_vor,A)	
	assert abs(result.imag - expect.imag) < 1e-14
	assert abs(result.real - expect.real) < 1e-14


def euler_ques3(delta_t,z_vor,A):
	result = z_vor[:]
	for j in range(0,1000):
		for i in range(0,len(A)):
			v = vortex_velocity(result[i],result,A)
			z_next = result[i] + delta_t * v
	        	plt.plot([result[i].real,z_next.real],[result[i].imag,z_next.imag]) 
			result[i] = z_next
	return result

def runge_ques3(delta_t,z_vor,A):
	result2,dummy = z_vor[:],z_vor[:]
	for j in range(0,1000):
		for i in range(0,len(A)):
			u = vortex_velocity(result2[i],result2,A)
			dummy[i] = result2[i]+ 0.5*delta_t * u
			u_mid = vortex_velocity(dummy[i],dummy,A)
			z_next = result2[i] + delta_t * u_mid
        		plt.plot([result2[i].real,z_next.real],[result2[i].imag,z_next.imag]) 
			result2[i],dummy[i] = z_next,z_next
	return result2

def problem3():
	plt.figure()
	plt.clf()
	euler_ques3(0.005,z_vor,A)
	axis_function()
	plt.title("Question 3: Using Euler Integration for t=0.01")
	plt.savefig('four.png')
	plt.figure()
	plt.clf()
	runge_ques3(0.01,z_vor,A)
	axis_function()
	plt.title("Question3: Using Runge-Kutta Integration for t =0.01")
	plt.savefig('five.png')

def actual_position_ques3(time,z_vor,A):
	result = z_vor[:]
	r = abs(z_vor[1]-z_vor[0])/2.0
	for i in range(0,2):
		u = (vortex_velocity(z_vor[i],z_vor,A))
		theta = abs(u)*time/r + np.arctan2(u.real,u.imag)
		result[i] = (z_vor[0] +z_vor[1])/2.0 + r*np.cos(theta) + 1j*r*np.sin(theta)
	return result

def compare_euler_runge(del_t,z_vor,A):
	log_e =[]	 #euler error
	log_r =[]	 #runge-kutta error
	for delta_t in del_t:
		time = 1.0/delta_t * 1000
		actual_result = actual_position_ques3(time,z_vor,A) #########
		euler_result = euler_ques3(1.0/delta_t,z_vor,A)
		runge_result = runge_ques3(1.0/delta_t,z_vor,A)
		error_euler,error_runge = 0.0,0.0
		for i in range(0,2):
			error_euler += abs(euler_result[i] - actual_result[i])
			error_runge += abs(runge_result[i] - actual_result[i]) 
		log_e.append(np.log(error_euler/2.0))
		log_r.append(np.log(error_runge/2.0))
	plt.figure()
	plt.clf()
	plt.plot(del_t,log_e, color="blue", linewidth=1.0, linestyle="-", label="euler")
	plt.plot(del_t,log_r, color="green", linewidth=1.0, linestyle="-", label="runge")
	plt.xlabel('1/delta_t')
	plt.ylabel('log(error)')
	plt.legend(loc='upper right')
	plt.title('Results from Question 3.')
	plt.savefig('six.png')

def write_html():
	message = """<html>
	<h1>Course: AE 625 - Particle Methods for Fluid Flow Simulation</h1>
	<h2>Assignment1</h2>
	<head><title>A1-13D100026</title></head>
	<body><p>Q1. Plot the streamlines and potential lines using the complex potential generated by this.</p></body>
	<img src='one.png'>
	<body><p>Q2. Now consider a set of tracer points starting at x=-2 (consider a line with say 10 points between y=-2 to 2).  Find the trajectory of these tracer points by integrating them given the velocity of the points.  Use both an Euler integrator and a Runge-Kutta second order to study the results.</p></body>
	<img src='two.png'></img>
	<img src='three.png'></img>
	<body><p>Q3. Consider the motion of two point vortices in isolation.  Consider two vortices of the same sign and strength = 2*pi.  Place them a unit distance apart.  Integrate the motion of these vortices in time and use the exact solution to test your implementation of the integrator.  Test the implementation of Euler and RK2 with this.  Study the convergence of the two integrators as you reduce the timestep used for the integrator. </p></body>
	<img src='four.png'></img>
	<img src='five.png'></img>
	<body><p>Convergence Plot</p></body>
	<img src='six.png'></img>
	</html>"""
	f.write(message)

if __name__ == '__main__':
	f = open('A1-13D100026.html','w')
	U = 1  		#freestream velocity      	
	x = np.linspace(-2.5,2.5, 250)
	y = np.linspace(-2.5, 2.5,250)
	X, Y = np.meshgrid(x, y)
	Q = [1,-1, -1]
	z = [-1+1j*0,1+1j*0, 0+1j*1]
	Q = [1,-1] 	#strength of source/sink
	z = [-1+1j*0,1+1j*0] #position of sources
	A=[2*np.pi,2*np.pi]	#strength of vortices 
	z_vor=[0+1j*0,1+1j*0]	#position of vortices
	del_t = [1000,500,200,100,50,20]
	test_phi()
	problem1()
	test_source_velocity()
	problem2()
	test_vortex_velocity()
	problem3()
	compare_euler_runge(del_t,z_vor,A)
	write_html()
	f.close()
webbrowser.open_new_tab('A1-13D100026.html')
