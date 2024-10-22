import random
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

num_of_vort =10000

def axis_function():
	ax = plt.gca()  
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_position(('data',0))
	ax.yaxis.set_ticks_position('left')
	ax.spines['left'].set_position(('data',0))

def exact_sol(z,t,h):
	return 1/(4*np.pi*mu*t)*np.exp(-abs(z)**2/(4*mu*t)) * h**2

def remesh(x,y,gamma_list,x_lim,y_lim,h):
	nx, ny = 2*x_lim/h +1, 2*y_lim/h +1
	xa = np.linspace(-x_lim, x_lim, nx)
	ya = np.linspace(-y_lim, y_lim, ny)
	xv, yv = np.meshgrid(xa, ya)
	x_new,y_new,gamma_new = [],[], []
	for m in range(len(x)):
		pos_x,pos_y,gamma0 = x[m],y[m],gamma_list[m] 
		x0 = pos_x - int(round(pos_x/h - 0.5))*h
		y0 = pos_y - int(round(pos_y/h - 0.5))*h
		if (h - pos_y % h  < 1e-5) and (h - pos_x % h < 1e-5):
			x_new += [pos_x-x0]
			y_new += [pos_y-y0]
			gamma_new += [gamma0]
		elif (h - pos_x % h < 1e-5):
			x_new += [pos_x-x0, pos_x-x0]
			y_new += [pos_y-y0, pos_y-y0+h]
			gamma_new += [gamma0*(1-y0) , gamma0*y0]
		elif (h - pos_y % h < 1e-5):
			x_new += [pos_x-x0, pos_x-x0+h]
			y_new += [pos_y-y0, pos_y-y0]
			gamma_new += [gamma0*(1-x0) , gamma0*x0]
		else:
			x_new += [pos_x-x0, pos_x-x0+h, pos_x-x0, pos_x-x0+h]
			y_new += [pos_y-y0, pos_y-y0, pos_y-y0+h, pos_y-y0+h]
			gamma_new += [gamma0*(1-y0)*(1-x0) , gamma0*x0*(1-y0) , gamma0*(1-x0)*y0,gamma0*x0*y0]
	plt.scatter(xv,yv,color='r')
	return x_new,y_new,gamma_new

def error(x_lim,y_lim,x_remesh,y_remesh,gamma_remesh):
	nx, ny = 2*x_lim/h +1, 2*y_lim/h +1
	xa = np.linspace(-x_lim, x_lim, nx)
	ya = np.linspace(-y_lim, y_lim, ny)
	xv, yv = np.meshgrid(xa, ya)
	error = [[0.0 for i in range(len(xa))] for j in range(len(ya))]
	approx_gamma = [[0.0 for i in range(len(xa))] for j in range(len(ya))]
	for k in range(len(x_remesh)):	
		i = int((x_remesh[k] + x_lim)*(nx-1)/2.0/x_lim)
		j = int((y_remesh[k] + y_lim)*(ny-1)/2.0/y_lim)
		approx_gamma[j][i] += gamma_remesh[k]
	for i,xv_pos in enumerate(xa):
		for j,yv_pos in enumerate(ya):
			z_pos = xv_pos +1j*yv_pos
			exact_gamma = exact_sol(z_pos,t,h)
			error[j][i]  = (approx_gamma[j][i] - exact_gamma )/exact_gamma
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot_surface(xv, yv, error , rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	plt.title("Relative Error Profile with grid size = 8x8")
	plt.savefig('error_10000.png')
	rms_error =sum((error[j][i])**2 for i in range(len(xa)) for j in range(len(xa)))
	return (np.sqrt((rms_error)/nx/ny))
			
def point_vortex_diffusion(gamma):
	x0,y0 = [0.0],[0.0]
	gamma_list = [1.0/num_of_vort]*num_of_vort
	x_next,y_next = x0*num_of_vort,y0*num_of_vort
	for j in range(0,num_of_vort):
		x_next[j] +=  np.random.normal(mean, sigma, 1)[0]
		y_next[j] +=  np.random.normal(mean, sigma, 1)[0] 
	plt.figure()
	plt.scatter(x_next,y_next)
	axis_function()
	plt.title("Vorticity distribution after 1 sec with num of vortices = 100000.(Using delta_t = 0.1 sec)")
	plt.savefig('vort_10000.png')
	plt.figure()
	x,y,gamma_list = remesh(x_next,y_next,gamma_list,x_lim,y_lim,h)
	plt.scatter(x,y)
	axis_function()
	plt.title(" Vorticity Distribution just after remeshing ( num of vortices = 100000)")
	plt.savefig('remesh_10000.png')
	return x,y,gamma_list

if __name__ == '__main__':
	gamma = 1.0
	mu = 0.1		#viscosity
	t =1.0				
	mean = 0.0
	delta_t = 0.1
	h = 0.05
	sigma = np.sqrt(2*mu*delta_t)
	x_lim,y_lim = int(round(3*sigma))+1 , int(round(3*sigma))+1
	x_remesh,y_remesh,gamma_remesh = point_vortex_diffusion(gamma)
	error(x_lim,y_lim,x_remesh,y_remesh,gamma_remesh)
	num_vort = [100,500,1000,3000,10000,100000]
	errorr = [0.0]*len(num_vort)
	for a in range(len(num_vort)):
		num_of_vort = num_vort[a]
		x_remesh,y_remesh,gamma_remesh = point_vortex_diffusion(gamma)	
		errorr[a] = error(x_lim,y_lim,x_remesh,y_remesh,gamma_remesh) 
	plt.figure()
	plt.plot(np.log(num_vort),errorr)
	plt.title("Variation of Error with Number of Vortices")
	plt.xlabel('log(#Vortices)')
	plt.ylabel('RMS Error')
	plt.savefig('error_plot.png')		
