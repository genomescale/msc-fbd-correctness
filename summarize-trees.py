import pymcmcstat.chain.ChainStatistics
import dendropy
import os
import csv
import numpy

n_topologies = 8
n_reps = 100
burnin = 5001
n_raw = 20000

configurations = ["full", "nodereheight", "coordinated", "msc", "sa", "updownless"]

marginals_filename = "estimated-marginals.csv"
marginals_file = open(marginals_filename, "w")
marginals_writer = csv.writer(marginals_file)

marginals_header = ["config", "rep_i", "topology_code", "variable", "lower_bound", "upper_bound", "marginal_probability"]
marginals_writer.writerow(marginals_header)

topology_output_filename = "topologies.csv"
topology_output_file = open(topology_output_filename, "w")
topology_output_writer = csv.writer(topology_output_file)

topology_header = ["config", "rep_i", "topology_code", "topology_count", "autocorrelation_time"]
topology_output_writer.writerow(topology_header)

for config in configurations:
	for rep_i in range(n_reps):
		config_filename = "%s.%03d.csv" % (config, rep_i)
		config_path = os.path.join(config, config_filename)
		config_file = open(config_path)
		config_reader = csv.reader(config_file)

		t_or_bins = numpy.zeros((n_topologies, 100))
		x1_bins = numpy.zeros((n_topologies, 100))
		x2_bins = numpy.zeros((n_topologies, 100))

		rep_topologies = numpy.zeros((n_raw, n_topologies), dtype = numpy.int32)

		for row_i, row in enumerate(config_reader):
			if row_i >= burnin:
				topology_code = int(row[0])

				rep_topologies[row_i - burnin][topology_code] = 1

				t_or = float(row[1])
				x1 = float(row[2])
				x2 = float(row[3])

				if t_or > 10.0:
					print([config, rep_i] + row)

				t_or_bin = min(int(t_or * 10.0), 99)
				x1_bin = int(x1 * 10.0)
				x2_bin = int(x2 * 10.0)

				t_or_bins[topology_code][t_or_bin] += 1
				x1_bins[topology_code][x1_bin] += 1
				x2_bins[topology_code][x2_bin] += 1

		topology_counts = numpy.sum(rep_topologies, axis = 0)
		autocorrelation_times = pymcmcstat.chain.ChainStatistics.integrated_autocorrelation_time(rep_topologies)[0]

		print(topology_counts)
		print(autocorrelation_times)

		for topology_code in range(n_topologies):
			topology_count = topology_counts[topology_code]
			autocorrelation_time = autocorrelation_times[topology_code]

			topology_output_row = [config, rep_i, topology_code, topology_count, autocorrelation_time]
			topology_output_writer.writerow(topology_output_row)

			for variable_bin in range(100):
				if topology_count > 0:
					t_or_count = t_or_bins[topology_code][variable_bin] / topology_count
					x1_count = x1_bins[topology_code][variable_bin] / topology_count
					x2_count = x2_bins[topology_code][variable_bin] / topology_count
				else:
					t_or_count = 0.0
					x1_count = 0.0
					x2_count = 0.0

				lower_bound = variable_bin / 10.0
				upper_bound = (variable_bin + 1) / 10.0

				marginals_writer.writerow([config, rep_i, topology_code, "t_or", lower_bound, upper_bound, t_or_count])
				marginals_writer.writerow([config, rep_i, topology_code, "x1", lower_bound, upper_bound, x1_count])
				marginals_writer.writerow([config, rep_i, topology_code, "x2", lower_bound, upper_bound, x2_count])
				marginals_file.flush()

marginals_file.close()
topology_output_file.close()
