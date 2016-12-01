import webbrowser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def axis_function():
	ax = plt.gca()  
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_position(('data',0))
	ax.yaxis.set_ticks_position('left')
	ax.spines['left'].set_position(('data',0))

def pos_vort_uniform(num_vort,b):
	y = np.linspace(-b/2.0,b/2.0,num_vort+1)
	dx = [0.0]*num_vort
	z_vor = [0+1j*0]*num_vort
	for i in range(0,num_vort):
		z_vor[i] = complex(0,(y[i+1]+y[i])/2)
		dx[i] = abs(y[i] - y[i+1])
	return z_vor,dx

def pos_vort_sine(num_vort,b):
	x = np.linspace(-np.pi/2.0,np.pi/2.0,num_vort+1)
	y = b*np.sin(x)/2.0
	dx = [0.0]*num_vort
	z_vor = [0+1j*0]*num_vort
	for i in range(0,num_vort):
		z_vor[i] = complex(0,(y[i+1]+y[i])/2)
		dx[i] = abs(y[i] - y[i+1])
	return z_vor,dx

def strength_of_vort(z_vort,b,dx):
	gamma0 =1
	A = [0.0]*len(z_vort)
	for i in range(0,len(z_vort)):
		A[i] = dx[i] * gamma0 * 4 * z_vort[i].imag/b**2 /np.sqrt(1-(2*z_vort[i].imag/b)**2)
	return A

def vortex_velocity_blob(zpos,z_vort, A,delta):
	vel = 0
	for i in range(0,len(A)):
		vel += A[i]*((z_vort[i].imag - zpos.imag) - 1j*(z_vort[i].real- zpos.real))/(2*np.pi*((abs(zpos - z_vort[i]))**2+ delta**2)) 
	return vel

def vortex_velocity_point(zpos,z_vort, A,delta):
	vel = 0
	for i in range(0,len(A)):
		if(np.any(z_vort[i] != zpos)):
			vel += 1j*A[i]/(2*np.pi*(zpos - z_vort[i])).conjugate()
	return vel

def runge_ques3(n,z_vor,A,vortex_velocity,delta):
	result2,dummy = z_vor[:],z_vor[:]
	delta_t = 0.01
	for j in range(0,n):
		for i in range(0,len(A)):
			u = vortex_velocity(result2[i],result2,A,delta)
			dummy[i] = result2[i]+ 0.5*delta_t * u
			u_mid = vortex_velocity(dummy[i],dummy,A,delta)
			z_next = result2[i] + delta_t * u_mid
        		result2[i],dummy[i] = z_next,z_next
	return result2

def vortex_sheet_rollup():
	result = [0.0 +1j*0.0]*num_vort
	plt.figure()		#1
	plt.clf()
	result = runge_ques3(1000,z_vor,A,vortex_velocity_point,0.001)
	for i in range(0,num_vort-1):
		plt.scatter([result[i].real,result[i+1].real],[result[i].imag,result[i+1].imag])
	axis_function()
	plt.title("Uniform distribution of 100 point vortices on vortex sheet")
	plt.savefig('one.png')	
	plt.figure()		#2
	plt.clf()
	result = runge_ques3(1000,z_vor2,A2,vortex_velocity_point,0.001)
	for i in range(0,num_vort-1):
		plt.scatter([result[i].real,result[i+1].real],[result[i].imag,result[i+1].imag])
	axis_function()
	plt.title("Sine distribution of 100 point vortices on vortex sheet")
	plt.savefig('two.png')	
	plt.figure()		#3
	plt.clf()
	result = runge_ques3(1000,z_vor,A,vortex_velocity_blob,0.001)
	for i in range(0,num_vort-1):
		plt.scatter([result[i].real,result[i+1].real],[result[i].imag,result[i+1].imag])
	axis_function()
	plt.title("Uniform distribution of 100 vortices blob on vortex sheet")
	plt.savefig('three.png')	
	plt.figure() 		#4
	plt.clf()
	result = runge_ques3(1000,z_vor2,A2,vortex_velocity_blob,0.001)
	for i in range(0,num_vort-1):
		plt.scatter([result[i].real,result[i+1].real],[result[i].imag,result[i+1].imag])
	axis_function()
	plt.title("Sine distribution of 100 point vortices blob on vortex sheet(delta = 0.001)")
	plt.savefig('four.png')	
	plt.figure() 		#7
	plt.clf()
	result = runge_ques3(0,z_vor2,A2,vortex_velocity_blob,0.001)
	for i in range(0,num_vort-1):
		plt.scatter([result[i].real,result[i+1].real],[result[i].imag,result[i+1].imag])
	axis_function()
	plt.savefig('seven.png')	
	plt.figure() 		#8
	plt.clf()
	result = runge_ques3(10,z_vor2,A2,vortex_velocity_blob,0.001)
	for i in range(0,num_vort-1):
		plt.scatter([result[i].real,result[i+1].real],[result[i].imag,result[i+1].imag])
	axis_function()
	plt.savefig('eight.png')
	plt.figure() 		#9
	plt.clf()
	result = runge_ques3(100,z_vor2,A2,vortex_velocity_blob,0.001)
	for i in range(0,num_vort-1):
		plt.scatter([result[i].real,result[i+1].real],[result[i].imag,result[i+1].imag])
	axis_function()
	plt.savefig('nine.png')	

def write_html():
	message = """<html>
	<h1>Course: AE 625 - Particle Methods for Fluid Flow Simulation</h1>
	<h2>Assignment 3</h2>
	<h3>Vortex sheet rolling using Point vortices and Vortex blobs</h3>
	<head><title>A3-13D100026</title></head>
	<body><p> The vortex sheet when discretized using point vortices, in uniform distribution or sine distrbution is as shown in figure.</p></body>
	<img src='one.png'></img>
	<img src='two.png'></img>
	<body><p> The vortex sheet when discretized using vortex blob methods, in uniform distribution or sine distrbution is as shown in figure.</p></body>
	<img src='three.png'></img>
	<img src='four.png'></img>
	<body><p> The evolution of the roll-up of vortex sheet discretized using 100 vortex blobs (Kransy blob method)occurs as follows:</p></body>
	<body><p> At t = 0 sec:</p></body>
	<img src='seven.png'></img>
	<body><p> At t = 0.1 sec:</p></body>
	<img src='eight.png'></img>
	<body><p> At t = 1 sec:</p></body>
	<img src='nine.png'></img>
	</html>"""
	f.write(message)
	
if __name__ == '__main__':
	f = open('A3-13D100026.html','w')
	num_vort = 100	
	z_vor,dx =pos_vort_uniform(num_vort,1.0)	
	z_vor2,dx2 =pos_vort_sine(num_vort,1.0)
	A = strength_of_vort(z_vor,1.0,dx)
	A2 = strength_of_vort(z_vor2,1.0,dx2)	 
	z_vor50,dx50 = pos_vort_sine(50,1.0)
	A50 = strength_of_vort(z_vor50,1.0,dx50)
	z_vor20,dx20 = pos_vort_sine(20,1.0)
	A20 = strength_of_vort(z_vor20,1.0,dx20)
	vortex_sheet_rollup()
	write_html()
	f.close()
webbrowser.open_new_tab('A3-13D100026.html')
