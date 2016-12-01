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
		value = q*(-2.0+1.5*q)/h*x/abs(x)/h
	elif(q <= 2.0):
		value = (-1.0/(2.0*h)*(2.0-q)**2.0)*x/abs(x)/h
	else:
		value = 0.0
	return value

def euler_integrate(delta_t,z_pos,velocity):
	z_next = z_pos + delta_t*velocity
        return z_next

def shock_tube():
	m = [rho_l*dx_l]*360
	v = [0.0]*360
	x_l = np.linspace(-0.5,0.0,320,endpoint=False)
	x_r = np.linspace(0.0,0.5,40)
	r = np.concatenate([x_l, x_r])
	p = [0.0]*360
	rho = [0.0]*360
	for i in range(0,320):
		rho[i] = rho_l
		p[i] = p_l
	for i in range(320,360):
		rho[i] = rho_r
		p[i] = p_r
	e = [0.0]*360
	for i in range(0,360):
		e[i] = p[i]/((gamma-1.0)*rho[i])
	time = 0.0
	epsilon= (0.01*h*h)
	while (time < final_t):
		print(time)
		time += dt
		dv_dt = [0.0]*360
		de_dt = [0.0]*360
		dx_dt = [0.0]*360
		rho = [0.0]*360 
		for i in range(0,360):
			for j in range(0,360):
				rho[i] += m[j]*cubic_kernel(r[i]-r[j], h)
		for i in range(0,360):
			for j in range(0,360):
				rij = r[i]-r[j]
				if (abs(rij) < 2*h):
					vij = v[i] - v[j]
					rhoij = (rho[j] + rho[i])/2.0
					if (rho[i] <=0 or rho[j] <= 0):
						print("found error1",rho[i],rho[j])
						c_ij = 0
					if (p[i] <=0 or p[j] <= 0):
						print("found error2",p[i],p[j],e[i],e[j])
						c_ij = 0
					else:
						c_ij = 0.5*(np.sqrt(gamma*p[i]/rho[i]) + np.sqrt(gamma*p[j]/rho[j]))
					if(vij*rij < 0.0 ):
						mu = h*vij*rij/(rij**2 + epsilon**2)
						pi = (-alpha*c_ij*mu + beta*mu*mu )*rhoij
					else:
						pi = 0.0
					dv_dt[i] += -m[j]*(p[j]/(rho[j]*rho[j]) + p[i]/(rho[i]*rho[i]) + pi)*cubic_der_kernel(rij,h)
					de_dt[i] += 0.5*m[j]*(p[j]/(rho[j]*rho[j]) + p[i]/(rho[i]*rho[i]) + pi)*vij*cubic_der_kernel(rij,h)	#change energy equation from pdf
					dx_dt[i] += 0.5*m[j]/rhoij*cubic_kernel(rij,h)*(-vij)
		for i in range(0,360):
			v[i] += dt*dv_dt[i]
			e[i] += dt*de_dt[i]
			p[i] = rho[i]*e[i]*(gamma-1.0)
			r[i] += dt*(v[i]+dx_dt[i])
	return p,rho,v,r,e

def exact():
	p_e,rho_e,u_e,e_e = [0.0]*8, [0.0]*8, [0.0]*8, [0.0]*8
	p_e[0],p_e[1],p_e[6],p_e[7] = p_l,p_l,p_r,p_r
	rho_e[0],rho_e[1],rho_e[6] ,rho_e[7] = rho_l,rho_l,rho_r,rho_r
	k = (gamma-1.0)/(gamma+1.0)
	l = (gamma-1.0)/(gamma+1.0)
	p2_guess = np.linspace(p_l, p_r, 1000)
	delta = np.zeros(len(p2_guess))
	for i,p2 in enumerate(p2_guess):
		u4 = (p2 - p_e[7])*np.sqrt((1.0-k)/(rho_e[7]*(p2 + k*p_e[7])))
		u2 = (p_e[0]**l - p2**l)* np.sqrt((1.0-k*k)*p_e[0]**(1.0/gamma)/(rho_l*k*k))
		delta[i] = abs(u4 - u2)
	p_e[2] = p2_guess[np.argmin(delta)]
	rho_e[2] = rho_e[0]*((p_e[2]/p_e[0])**(1.0/gamma))
	rho_e[4] = rho_e[6]*(p_e[2] + k*p_e[6])/(p_e[6] + k*p_e[2])
	u_e[2] = u_e[7] + (p_e[2] - p_e[7])/np.sqrt(0.5*rho_e[6]*((gamma+1.0)*p_e[2] + (gamma-1.0)*p_e[6]))
	u_e[3] ,u_e[4] ,u_e[5] =  u_e[2],u_e[2],u_e[2]
	rho_e[3],rho_e[5] = rho_e[2],rho_e[4]  
	p_e[3],p_e[4], p_e[5] = p_e[2],p_e[2],p_e[2]
	e_e = [p_e[i]/((gamma-1.0)*rho_e[i]) for i in range(0,8)]
	x_e = [-0.5, -0.25, 0.0, 0.2, 0.2, 0.4, 0.4, 0.5]
	return p_e ,rho_e ,u_e ,e_e ,x_e

def shock_tube_plots():	
	p_exact,rho_exact,v_exact,e_exact,x_exact = exact()
	p,rho,v,r,e = shock_tube()
	plt.figure()
	plt.plot(r, v, label='SPH')
    	plt.plot(x_exact, v_exact, label='Exact')
    	plt.xlim(-0.5, 0.5)
	plt.title('Velocity vs x')
	plt.legend()
	plt.savefig('vel.png')
	plt.figure()
	plt.plot(r, rho, label='SPH')
    	plt.plot(x_exact, rho_exact, label='Exact')
    	plt.xlim(-0.5, 0.5)
	plt.title('Density vs x')
	plt.legend()
	plt.savefig('density.png')
	plt.figure()
	plt.plot(r, p, label='SPH')
	plt.plot(x_exact, p_exact, label='Exact')
	plt.title('Pressure vs x')
	plt.xlim(-0.5, 0.5)
	plt.legend()
	plt.savefig('pressure.png')
	plt.figure()
	plt.plot(r, e, label='SPH')
	plt.plot(x_exact, e_exact, label='Exact')
	plt.xlim(-0.5, 0.5)
	plt.title('Energy vs x')
	plt.legend()
	plt.savefig('energy.png')

if __name__ == '__main__':
	rho_l = 1.0
	rho_r = 0.125
	p_l = 1.0
	p_r = 0.1
	u_l = 0.0
	u_r = 0.0
	dx_l = 0.0015625
	dx_r = 0.0125
	alpha =1.0
	beta = 1.0
	gamma =1.4
	h = 2.0*dx_r
	dt = 1e-4
	final_t = 0.2
	shock_tube_plots()
