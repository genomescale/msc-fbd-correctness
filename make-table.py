import csv
import numpy
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

n_raw = 20000.0

topologies_filename = "topologies.csv"
topologies_file = open(topologies_filename)
topologies_reader = csv.reader(topologies_file)

inside_counts = {}

topology_output_header = ["config", "rep_i", "topology_code", "topology_count", "autocorrelation_time"]

for row_i, row in enumerate(topologies_reader):
	if row_i >= 1:
		config = row[0]
		topology_code = int(row[2])
		topology_count = float(row[3])
		autocorrelation_time = float(row[4])
		ess = int(round(n_raw / max(1.0, autocorrelation_time)))
		# ess = 1000.0

		if config not in inside_counts:
			inside_counts[config] = [0] * 8

		probability_interval = scipy.stats.binom.ppf([0.025, 0.975], n = ess, p = true_probs[topology_code]) / ess

		if probability_interval[0] < (topology_count / n_raw) < probability_interval[1]:
			inside_counts[config][topology_code] += 1

output_filename = "table.csv"
output_file = open(output_filename, "w")
output_writer = csv.writer(output_file)
for config in sorted(inside_counts):
	output_writer.writerow([config] + inside_counts[config])

output_file.close()
