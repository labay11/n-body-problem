from __future__ import division
import sys
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation

COLORS = ['b', 'g', 'r', 'y', 'm', 'k', 'c']

G = 6.67259e-11 # gravitational constant

"""
NAMES  = ['Sun',		'Earth', 		'Moon',			'Mercury', 		'Venus',		'Mars', 		'Jupiter', 		'Saturn',		'Uranus',		'Neptune',		'Pluto']
MASSES = [1.98855e30, 	5.97237e24, 	7.3477e20,		0.33e24,		4.8675e24, 		0.6417e24,		1898.2e24,		584.34e24,		86.813e24,		102.413e24,		13.03e21]
POS_0  = [(0, 0),		(147.09e9, 0),	(147.44e9, 0),	(46e9, 0),		(107.48e9, 0),	(0, 249.23e9),	(0, 816.62e9),	(1.352e12, 0),	(0, 3.003e12),	(4.444e12, 0),	(7.375e12, 0)] # initial positions
VEL_O  = [(0, 0),		(0, 30.29e3),	(0, 31.395e3),	(0, 58.98e3),	(0, 35.26e3),	(21.97e3, 0),	(12.44e3, 0),	(0, 10.18e3),	(6.49e3, 0),	(0, 5.5e3),		(0, 3.71e3)] # initial velocities

PLANETS = len(MASSES)
"""

NAMES = []
MASSES = []
POS_0 = []
VEL_O = []
PLANETS = None # number of particles

def read_particles_data(file_name):
	global NAMES, MASSES, POS_0, VEL_O, PLANETS
	def _format_vector(vec):
		"""
			Utility function that formats an input string
			in the format: '%f;%f' into a tuple (x, y)

			::vec:: input string
			::return:: vector as tuple
		"""
		data = vec.split(';')
		return (float(data[0].strip()), float(data[1].strip()))

	with open(file_name, 'r') as file:
		content = file.read()
		particles = content.split('\n\n')

		PLANETS = len(particles)

		for part in particles:
			print part
			data = part.split('\n')
			name = data[0].strip()
			mass = float(data[1].strip())
			r0 = _format_vector(data[2].strip())
			v0 = _format_vector(data[3].strip())

			NAMES.append(name)
			MASSES.append(mass)
			POS_0.append(r0)
			VEL_O.append(v0)


STEPS = 5000
T_MAX, T_MIN = 365 * 24 * 3600.0, 0.
h = (T_MAX - T_MIN) / STEPS

def draw():
	#			  time	 	planet     x,y
	r = np.zeros((STEPS, 	PLANETS, 	2)) # positions
	#				planet   v_x,v_y
	v_s = np.zeros((PLANETS,	2)) # velocities at step s

	# apply initial conditions
	for p in xrange(PLANETS):
		for i in xrange(2):
			r[0][p][i] = POS_0[p][i]
			v_s[p][i] = VEL_O[p][i]

	fig, ax = plt.subplots()
	lines = []

	def init():
		if len(lines) > 0:
			# prevents from executing a second time
			return lines

		ax.set_xlim(-5.6e12, 5.6e12)
		ax.set_ylim(-5.6e12, 5.6e12)

		for p in xrange(PLANETS):
			line = Line2D([], [], label=NAMES[p], lw=2, color=COLORS[p % len(COLORS)])
			ax.add_line(line)
			lines.append(line)

		ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
	           ncol=3, mode="expand", borderaxespad=0.)

		return lines

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
		for k in xrange(PLANETS):
			if k == p:
				continue

			v += MASSES[k] * (rs[p] - rs[k]) / mod3(rs[p] - rs[k])

		return - G * v

	alpha = [0, 0, 0.5, 0.5, 1] # first index for simplicity
	rk4 = np.zeros((5, PLANETS, 2))

	s = 0
	def update(t):
		global s, r, v_s
		if s > 0 and s < 5000:
			for j in xrange(1, 5):
				for p in xrange(PLANETS):
					rk4[j][p] = g(p, r[s - 1] + alpha[j] * h * rk4[j - 1])

			v_s += (h / 6) * (rk4[1] + 2 * (rk4[2] + rk4[3]) + rk4[4])
			r[s] = r[s - 1] + h * v_s

		for p in xrange(PLANETS):
			lines[p].set_data(r[:s,p,0], r[:s,p,1])

		s += 1

		return lines

	ani = FuncAnimation(fig, update, frames=np.linspace(0, T_MAX, STEPS),
	                    init_func=init, blit=True, interval=0.5)
	plt.show()

"""for s in xrange(1, STEPS):
	for j in xrange(1, 5):
		for p in xrange(PLANETS):
			rk4[j][p] = g(p, r[s - 1] + alpha[j] * h * rk4[j - 1])

	v_s += (h / 6) * (rk4[1] + 2 * (rk4[2] + rk4[3]) + rk4[4])
	r[s] = r[s - 1] + h * v_s

for p in xrange(PLANETS):
	plt.plot(r[:,p,0], r[:,p,1], label=NAMES[p], lw=2, color=COLORS[p % len(COLORS)])

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)
#plt.xlim(-1.6e11, 1.6e11)
#plt.ylim(-1.6e11, 1.6e11)
plt.savefig("solar_system.pdf")
plt.show()"""

if __name__ == '__main__':
	file_name = None
	if len(sys.argv) > 1:
		file_name = sys.argv[1]
	else:
		file_name = 'ss_data.txt'

	read_particles_data(file_name)
	draw()
