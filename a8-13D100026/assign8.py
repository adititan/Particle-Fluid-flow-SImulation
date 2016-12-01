import matplotlib.pyplot as plt
import numpy as np

def cubic_kernel(x,h):
	q = abs(x)/h
	if (q <= 1.0):
		value = 2.0/(3.0*h)*(1.0-(3.0/2.0*q*q*(1-q/2.0)))
	elif(q <= 2.0):
		value = 1.0/(6.0*h)*(2.0-q)**3
	else:
		value = 0.0
	return value

def cubic_der_kernel(x,h):
	q = abs(x)/h
	if (q == 0.0):
		value = 0.0
	elif (q <= 1.0):
		value = -q*(2.0-1.5*q)/h*x/abs(x)/h
	elif(q <= 2.0):
		value = -1*(1.0/(2.0*h)*(2.0-q)**2.0)*x/abs(x)/h
	else:
		value = 0.0
	return value

def gaussian_kernel(x,h):
	q = abs(x)/h
	if (q <= 3.0):
		value = 1/(np.sqrt(np.pi)*h)*np.exp(-q*q)
	elif (q > 3.0):
		value = 0.0
	return value

def gaussian_der_kernel(x,h):
	q = abs(x)/h
	if (q == 0.0):
		value = 0.0
	elif (q <= 3.0):
		value = -2.0*q/(np.sqrt(np.pi)*h)*np.exp(-q*q)*x/abs(x)/h
	elif (q > 3.0):
		value = 0.0
	return value

def sine_func(x):
	return -1*np.sin(np.pi*x)

def cos_func(x):
	return -1*np.pi*np.cos(np.pi*x)

def sph_approximate_1d(xi,func,dx,hdx_g,hdx_c):
	h_g = hdx_g * dx
	h_c = hdx_c * dx
	fi = [func(xi[i]) for i in range(0,int(2/dx))]
	f_cubic = [0.0]*int(2/dx)
	f_gauss = [0.0]*int(2/dx)
	error_g = [0.0]*int(2/dx)
	error_c = [0.0]*int(2/dx)
	for i in range(0,int(2/dx)):
		for j in range(0,int(2/dx)):
			f_gauss[i] += fi[j]*gaussian_kernel(xi[i] - xi[j],h_g)*dx
			f_cubic[i] += fi[j]*cubic_kernel(xi[i] - xi[j],h_c)*dx
		error_g[i] = (f_gauss[i] - fi[i])**2
		error_c[i] = (f_cubic[i] - fi[i])**2
	error_g_l2 = np.sqrt(sum(error_g)/len(error_g))
	error_c_l2 = np.sqrt(sum(error_c)/len(error_g))
	plt.figure()
	plt.plot(xi,fi,label='exact')
	plt.plot(xi,f_gauss,label='gaussian')
	plt.plot(xi,f_cubic,label='cubic')
	plt.legend()
	plt.title('Approximation of Sine function with SPH')
	plt.savefig('one.png')
	return error_g_l2,error_c_l2

def sph_der_approximate_1d(xi,func,dx,hdx_g,hdx_c):
	h_g = hdx_g * dx
	h_c = hdx_c * dx
	fi = [func(xi[i]) for i in range(0,int(2/dx))]
	fi_cos = [cos_func(xi[i]) for i in range(0,int(2/dx))]
	f_der_cubic = [0.0]*int(2/dx)
	f_der_gauss = [0.0]*int(2/dx)
	error_g = [0.0]*int(2/dx)
	error_c = [0.0]*int(2/dx)
	for i in range(0,int(2/dx)):
		for j in range(0,len(xi)):
			f_der_gauss[i] += fi[j]*gaussian_der_kernel(xi[i] - xi[j],h_g)*dx
			f_der_cubic[i] += fi[j]*cubic_der_kernel(xi[i] - xi[j],h_c)*dx
		error_g[i] = (f_der_gauss[i] - fi_cos[i])**2
		error_c[i] = (f_der_cubic[i] - fi_cos[i])**2
	error_g_l2 = np.sqrt(sum(error_g)/len(error_g))
	error_c_l2 = np.sqrt(sum(error_c)/len(error_g))
	plt.figure()
	plt.plot(xi,fi_cos,label='exact')
	plt.plot(xi,f_der_gauss,label='gaussian')
	plt.plot(xi,f_der_cubic,label='cubic')
	plt.legend()
	plt.title('Derivative Approximation of Sine function with SPH')
	plt.savefig('two.png')
	return error_g_l2,error_c_l2

def noise_addition(xi):
	x_new = [0.0]*len(xi)
	for i in range(0,len(xi)):
		x_new[i] = xi[i] + np.random.uniform(-0.01, 0.01,1)
	return x_new
	
if __name__ == '__main__':
	dx = 0.1
	hdx_g = 0.7
	hdx_c = 0.7
	hdx_gg = 1.0
	hdx_cc = 1.0
	dx_list = [0.01,0.02,0.05,0.1,0.2]		#plotting for sine approx
	errorr_g = [0.0]*len(dx_list)
	errorr_c = [0.0]*len(dx_list)
	for k in range(0,len(dx_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx_list[k]))
		errorr_g[k],errorr_c[k] = sph_approximate_1d(xi,sine_func,dx_list[k],hdx_g,hdx_c)
	plt.figure()
	plt.plot(dx_list,errorr_c,label='cubic')
	plt.plot(dx_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('dx')
	plt.ylabel('L2 Error')
	plt.title('Error Vs. dx in Derivative Approximation of Sine function')
	plt.savefig('error_one.png')
	h_list = [1.0,2.0,3.0,4.0,5.0]
	errorr_g = [0.0]*len(dx_list)
	errorr_c = [0.0]*len(dx_list)
	for k in range(0,len(dx_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx))
		errorr_g[k],errorr_c[k] = sph_approximate_1d(xi,sine_func,dx,h_list[k],h_list[k])
	plt.figure()
	plt.plot(h_list,errorr_c,label='cubic')
	plt.plot(h_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('hdx')
	plt.ylabel('L2 Error')
	plt.title('Error Vs. hdx in Approximation of Sine function')
	plt.savefig('error_two.png')
	dx_list = [0.01,0.02,0.05,0.1,0.2]		#plotting for approx of derivative of sine
	errorr_g = [0.0]*len(dx_list)
	errorr_c = [0.0]*len(dx_list)
	for k in range(0,len(dx_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx_list[k]))
		errorr_g[k],errorr_c[k] = sph_der_approximate_1d(xi,sine_func,dx_list[k],hdx_gg,hdx_cc)
	plt.figure()
	plt.plot(dx_list,errorr_c,label='cubic')
	plt.plot(dx_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('dx')
	plt.ylabel('L2 Error')
	plt.title('Error Vs. dx in Derivative Approximation of Sine function')
	plt.savefig('error_three.png')
	h_list = [1.0,2.0,3.0,4.0,5.0]
	errorr_g = [0.0]*len(dx_list)
	errorr_c = [0.0]*len(dx_list)
	for k in range(0,len(dx_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx))
		errorr_g[k],errorr_c[k] = sph_der_approximate_1d(xi,sine_func,dx,h_list[k],h_list[k])
	plt.figure()
	plt.plot(h_list,errorr_c,label='cubic')
	plt.plot(h_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('hdx')
	plt.ylabel('L2 Error')
	plt.title('Error Vs. hdx in Derivative Approximation of Sine function')
	plt.savefig('error_four.png')
	dx_list = [0.01,0.02,0.05,0.1,0.2]			#plotting of sine approx with noise
	errorr_g = [0.0]*len(dx_list)
	errorr_c = [0.0]*len(dx_list)
	for k in range(0,len(dx_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx_list[k]))
		xi = noise_addition(xi)
		errorr_g[k],errorr_c[k] = sph_approximate_1d(xi,sine_func,dx_list[k],hdx_g,hdx_c)
	plt.figure()
	plt.plot(dx_list,errorr_c,label='cubic')
	plt.plot(dx_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('dx')
	plt.ylabel('L2 Error')
	plt.title('Error Vs. dx in Approximation of Sine function with Noise')
	plt.savefig('error_five.png')
	h_list = [1.0,2.0,3.0,4.0,5.0]
	errorr_g = [0.0]*len(h_list)
	errorr_c = [0.0]*len(h_list)
	for k in range(0,len(h_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx))
		xi = noise_addition(xi)
		errorr_g[k],errorr_c[k] = sph_approximate_1d(xi,sine_func,dx,h_list[k],h_list[k])
	plt.figure()
	plt.plot(h_list,errorr_c,label='cubic')
	plt.plot(h_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('hdx')
	plt.title('Error Vs. hdx in Approximation of Sine function with Noise')
	plt.ylabel('L2 Error')
	plt.savefig('error_six.png')
	hdx_gg = 1.0
	hdx_cc = 1.0
	dx_list = [0.01,0.02,0.05,0.1,0.2]		#plotting of sine derivative approx with noise
	errorr_g = [0.0]*len(dx_list)
	errorr_c = [0.0]*len(dx_list)
	for k in range(0,len(dx_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx_list[k]))
		xi = noise_addition(xi)
		errorr_g[k],errorr_c[k] = sph_der_approximate_1d(xi,sine_func,dx_list[k],hdx_gg,hdx_cc)
	plt.figure()
	plt.plot(dx_list,errorr_c,label='cubic')
	plt.plot(dx_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('dx')
	plt.ylabel('L2 Error')
	plt.title('Error Vs. dx in Derivative Approximation of Sine function with Noise')
	plt.savefig('error_seven.png')
	h_list = [1.0,2.0,3.0,4.0,5.0]
	errorr_g = [0.0]*len(h_list)
	errorr_c = [0.0]*len(h_list)
	for k in range(0,len(h_list)):
		xi = np.linspace(-1.0,1.0,int(2/dx))
		xi = noise_addition(xi)
		errorr_g[k],errorr_c[k] = sph_der_approximate_1d(xi,sine_func,dx,h_list[k],h_list[k])
	plt.figure()
	plt.plot(h_list,errorr_c,label='cubic')
	plt.plot(h_list,errorr_g,label='gaussian')
	plt.legend()
	plt.xlabel('hdx')
	plt.ylabel('L2 Error')
	plt.title('Error Vs. hdx in Derivative Approximation of Sine function with Noise')
	plt.savefig('error_eight.png')
	xi = np.linspace(-1.0,1.0,int(2/dx))
	sph_approximate_1d(xi,sine_func,dx,hdx_g,hdx_c)	
	sph_der_approximate_1d(xi,sine_func,dx,hdx_gg,hdx_cc)
