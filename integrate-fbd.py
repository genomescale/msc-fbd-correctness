import csv
import numpy
import scipy.integrate

birth_rate = 2.0
death_rate = 1.0
sampling_rate = 0.5
removal_rate = 0.9

y1 = 2.0
y2 = 1.0
y3 = 0.0

rate_diff = birth_rate - death_rate - sampling_rate

c1 = numpy.sqrt(rate_diff * rate_diff + 4 * birth_rate * sampling_rate)
c2 = rate_diff / c1

def q(x):
	return 4.0 / (2.0 * (1.0 - c2 * c2) +
		numpy.exp(-c1 * x) * (1 - c2) * (1 - c2) +
		numpy.exp(c1 * x) * (1 + c2) * (1 + c2))

def p0(x):
	return (birth_rate + death_rate + sampling_rate + c1 * (
		(numpy.exp(-c1 * x) * (1 - c2) - (1 + c2)) /
		(numpy.exp(-c1 * x) * (1 - c2) + (1 + c2)))) / (2.0 * birth_rate)

y_partial = 1.0

for y_i in [y1, y2, y3]:
	y_partial *= sampling_rate * (removal_rate + (1 - removal_rate) * p0(y_i)) / q(y_i)

def fbd(t_or, x1, x2):
	return (q(t_or) *
		2.0 * birth_rate * q(x1) *
		2.0 * birth_rate * q(x2) *
		y_partial / (6.0 * (1.0 - p0(t_or))))

global_upper = 10.0

marginal_likelihood = scipy.integrate.tplquad(fbd, y2, global_upper, lambda x2: max(x2, y1), lambda x2: global_upper, lambda x2, x1: x1, lambda x2, x1: global_upper)[0]

output_filename = "true-marginals.csv"
output_file = open(output_filename, "w")
output_writer = csv.writer(output_file)

output_header = ["variable", "lower_bound", "upper_bound", "marginal_probability"]
output_writer.writerow(output_header)

# Bifurcation and origin time probability masses
for i in range(100):
	lower_bound = i / 10.0
	upper_bound = (i + 1) / 10.0

	if i < 20:
		t_or_mass = 0.0
		x1_mass = 0.0
	else:
		t_or_mass = scipy.integrate.tplquad(fbd, y2, upper_bound, lambda x2: max(x2, y1), lambda x2: upper_bound, lambda x2, x1: max(x1, lower_bound), lambda x2, x1: upper_bound)[0]
		x1_mass = scipy.integrate.tplquad(fbd, y2, upper_bound, lambda x2: max(x2, lower_bound), lambda x2: upper_bound, lambda x2, x1: x1, lambda x2, x1: global_upper)[0]

	if i < 10:
		x2_mass = 0.0
	else:
		x2_mass = scipy.integrate.tplquad(fbd, lower_bound, upper_bound, lambda x2: max(x2, y1), lambda x2: global_upper, lambda x2, x1: x1, lambda x2, x1: global_upper)[0]

	output_writer.writerow(["t_or", lower_bound, upper_bound, t_or_mass / marginal_likelihood])
	output_writer.writerow(["x1", lower_bound, upper_bound, x1_mass / marginal_likelihood])
	output_writer.writerow(["x2", lower_bound, upper_bound, x2_mass / marginal_likelihood])

	output_file.flush()

output_file.close()
