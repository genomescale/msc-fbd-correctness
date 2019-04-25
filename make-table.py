import csv
import scipy.stats

true_probs = [
	0.778327,
	0.078642,
	0.038657,
	0.043189,
	0.043189,
	0.004135,
	0.006930,
	0.006930,
]

topologies_filename = "topologies.csv"
topologies_file = open(topologies_filename)
topologies_reader = csv.reader(topologies_file)

inside_counts = {}

for row in topologies_reader:
	config = row[0]
	ess = int(round(float(row[2])))

	if config not in inside_counts:
		inside_counts[config] = [0] * 8

	for i in range(8):
		probability_interval = scipy.stats.binom.ppf([0.025, 0.975], n = ess, p = true_probs[i]) / ess
		observed_probability = float(row[i + 3])

		if probability_interval[0] < observed_probability < probability_interval[1]:
			inside_counts[config][i] += 1

print(inside_counts)
