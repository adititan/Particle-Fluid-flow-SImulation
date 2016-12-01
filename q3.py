import numpy as np
import matplotlib.pyplot as plt
def vortex_velocity(z, z_vor, gamma):
    return (-1j*gamma/(2*np.pi*(z - z_vor))).conjugate()

def test_vortex_velocity():
    z = complex(2.0, 1)
    z_vor = complex(1., 1.0)
    expect = 1j

    result = vortex_velocity(z, z_vor, np.pi*2)
    assert abs(result.imag - expect.imag) < 1e-14
    assert abs(result.real - expect.real) < 1e-14

def vortex_array_velocity(pos, gamma):
    vel = np.zeros_like(pos)
    for i, z_i in enumerate(pos):
        for j, z_j in enumerate(pos):
            if i != j:
                vel[i] += vortex_velocity(z_i, z_j, gamma[j])
    return vel

def euler_integrate(pos, gamma, dt, tf):
    result = [pos]
    t = 0.0
    while t < tf:
        vel = vortex_array_velocity(pos, gamma)
        pos += vel*dt
        result.append(pos.copy())
        t += dt
    return np.asarray(result)

def problem3():
    vor_z = np.asarray([complex(0.5, 0), complex(-0.5, 0)])
    gamma = np.asarray([2*np.pi, 2*np.pi])
    result = euler_integrate(vor_z, gamma, dt=0.01, tf=5.0)
    plt.plot(result)
    return result

if __name__ == '__main__':
    test_vortex_velocity()
    problem3()
	
