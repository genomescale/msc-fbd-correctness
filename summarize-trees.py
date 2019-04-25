import dendropy
import tensorflow
import tensorflow_probability
import os
import csv
import numpy

session = tensorflow.InteractiveSession()

n_reps = 100
burnin = 250

configurations = ["full", "nodereheight", "coordinated", "msc", "sa"]
topology_counts = {config: [] for config in configurations}
rep_ess_tensors = {config: [] for config in configurations}

marginals_filename = "estimated-marginals.csv"
marginals_file = open(marginals_filename, "w")
marginals_writer = csv.writer(marginals_file)

marginals_header = ["config", "rep_i", "topology_code", "variable", "lower_bound", "upper_bound", "marginal_probability"]
marginals_writer.writerow(marginals_header)

for config in configurations:
	for rep_i in range(n_reps):
		config_filename = "%s.%03d.csv" % (config, rep_i)
		config_path = os.path.join(config, config_filename)
		config_file = open(config_path)
		config_reader = csv.reader(config_file)

		rep_topologies = []
		rep_counts = numpy.zeros(8)
		marginal_bins = numpy.zeros((3, 8, 100))
		for row_i, row in enumerate(config_reader):
			if row_i > burnin:
				topology_code = int(row[0])
				rep_counts[topology_code] += 1
				rep_topologies.append(topology_code)

				t_or_bin = int(float(row[1]) * 10.0)
				x1_bin = int(float(row[2]) * 10.0)
				x2_bin = int(float(row[3]) * 10.0)

				if t_or_bin >= 100:
					print(row)
					print(config, rep_i, row_i)
				assert t_or_bin < 100

				marginal_bins[0][topology_code][t_or_bin] += 1
				marginal_bins[1][topology_code][x1_bin] += 1
				marginal_bins[2][topology_code][x2_bin] += 1

		for topology_code in range(8):
			for variable_bin in range(100):
				if rep_counts[topology_code] > 0:
					t_or_count = marginal_bins[0][topology_code][variable_bin] / rep_counts[topology_code]
					x1_count = marginal_bins[1][topology_code][variable_bin] / rep_counts[topology_code]
					x2_count = marginal_bins[2][topology_code][variable_bin] / rep_counts[topology_code]
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

		topologies_tensor = tensorflow.Variable(rep_topologies, dtype = tensorflow.float64)
		ess_tensor = tensorflow_probability.mcmc.effective_sample_size(topologies_tensor)

		topology_counts[config].append(rep_counts)
		rep_ess_tensors[config].append(ess_tensor)

marginals_file.close()

session.run(tensorflow.global_variables_initializer())

topology_output_filename = "topologies.csv"
topology_output_file = open(topology_output_filename, "w")
topology_output_writer = csv.writer(topology_output_file)

for config in configurations:
	for rep_i in range(n_reps):
		ess = rep_ess_tensors[config][rep_i].eval()
		output_row = [config, rep_i, ess]
		for topology_code in range(8):
			topology_prob = topology_counts[config][rep_i][topology_code] / 1000.0
			output_row.append(topology_prob)

		topology_output_writer.writerow(output_row)
		topology_output_file.flush()

topology_output_file.close()
