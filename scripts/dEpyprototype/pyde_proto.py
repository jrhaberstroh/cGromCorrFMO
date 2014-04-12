import numpy as np
import math

filename = "de_cgrom.csv"

dat = []
with open(filename) as f:
	print f.readline()
	for line in f:
		line_float = np.array([float(x) for x  in line.split()[1:]])
		print line_float
		dat.append(line_float)

energy_tot = 0
energy_tot2 = 0
for x1 in dat:
	for x2 in dat:
		if x1 is not x2:
			dist_nm = 0.;
			for i in xrange(3):
				dist_nm += (x1[i] - x2[i])**2
			dist_nm = math.sqrt(dist_nm)
			if dist_nm < .2:
				print dist_nm
			assert dist_nm != 0
			esConst_kCalnm_e2 = 33.
			q0_1 = x1[3];
			qE_1 = x1[4];
			dq_1 = x1[4] - x1[3];
			q0_2 = x2[3];
			qE_2 = x2[4];
			dq_2 = x2[4] - x2[3];
			energy_tot += .5 * esConst_kCalnm_e2 / dist_nm * (qE_1 * qE_2 - q0_1 * q0_2)
			energy_tot2 += esConst_kCalnm_e2 / dist_nm * (dq_1 * q0_2 + .5 * dq_1 * dq_2)


print "Total dE using qE: ",energy_tot, "kCal/mol"
print "Total dE using dq: ",energy_tot2, "kCal/mol"
print "Program done"
