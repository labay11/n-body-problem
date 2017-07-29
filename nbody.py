from __future__ import division
import sys
import numpy as np 
import matplotlib.pyplot as plt

COLORS = ['b', 'g', 'r', 'y', 'm', 'k', 'c']

G = 6.67259e-11 # gravitational constant

class Particle(object):

	def __init__(self, name, mass, r0, v0):
		self.name = name
		self.mass = mass
		self.r0 = r0
		self.v0 = v0

def read_particles_data(file_name):
	global NAMES, MASSES, POS_0, VEL_O, PLANETS
	def _format_vector(vec):
		"""
			Utility function that formats an input string
			in the format: '%f;%f' into a tuple (x, y)

			::vec:: input string
			::return:: vector as tuple
		"""
		xy = vec.split(';')
		return (float(xy[0].strip()), float(xy[1].strip()))

	data = []

	with open(file_name, 'r') as file:
		content = file.read()
		particles = content.split('\n\n')

		for part in particles:
			p = part.split('\n')
			name = p[0].strip()
			mass = float(p[1].strip())
			r0 = _format_vector(p[2].strip())
			v0 = _format_vector(p[3].strip())

			data.append(Particle(name, mass, r0, v0))

	return data

def draw(data, STEPS, T_MAX):
	PARTICLES = len(data)
	h = T_MAX / STEPS # step size

	m = np.zeros(PARTICLES) # masses
	#			  time	 	planet     x,y
	r = np.zeros((STEPS, 	PARTICLES, 	2)) # positions
	#				planet   v_x,v_y
	v_s = np.zeros((PARTICLES,	2)) # velocities at step s

	# apply initial conditions
	for p in xrange(PARTICLES):
		m[p] = data[p].mass
		for i in xrange(2):
			r[0][p][i] = data[p].r0[i]
			v_s[p][i] = data[p].v0[i]

	def mod3(v):
		""" returns the lenght of a 2d-vector to the third power """
		return np.power(v[0]*v[0] + v[1]*v[1], 1.5)

	def g(p, rs):
		""" 
			return the intensity of the gravitational force in particle p due to the others
				given by Newton's law of universal gravitation:
					g_p = - G * sum_k [m_k * (r_p - r_k) / |r_p - r_k|^2 where k != p]
		"""
		v = [0, 0]
		for k in xrange(PARTICLES):
			if k == p:
				continue

			v += m[k] * (rs[p] - rs[k]) / mod3(rs[p] - rs[k])

		return - G * v

	alpha = [0, 0, 0.5, 0.5, 1] # first index for simplicity
	rk4 = np.zeros((5, PARTICLES, 2))

	for s in xrange(1, STEPS):
		for j in xrange(1, 5):
			for p in xrange(PARTICLES):
				rk4[j][p] = g(p, r[s - 1] + alpha[j] * h * rk4[j - 1])

		v_s += (h / 6) * (rk4[1] + 2 * (rk4[2] + rk4[3]) + rk4[4])
		r[s] = r[s - 1] + h * v_s

	# plot the trajectory
	for p in xrange(PARTICLES):
		plt.plot(r[:,p,0], r[:,p,1], label=data[p].name, lw=2, color=COLORS[p % len(COLORS)])

	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
	           ncol=3, mode="expand", borderaxespad=0.) # legend location at top
	# set x, y plot limits for convenience
	#plt.xlim(-1.6e11, 1.6e11)
	#plt.ylim(-1.6e11, 1.6e11)
	plt.savefig("nbody.pdf")
	plt.show()

if __name__ == '__main__':
	file_name = None
	if len(sys.argv) > 1:
		file_name = sys.argv[1]
	else:
		file_name = 'ss_data.txt'

	data = read_particles_data(file_name)
	draw(data, 4000, 365 * 24 * 3600.0)