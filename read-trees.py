import csv
import dendropy
import sys
import subprocess
import os

def calculate_sa_code(tree):
	tree.update_bipartitions()

	sa_code = 0
	for leaf in tree.leaf_nodes():
		if leaf.edge.length == 0.0:
			sa_code += leaf.leafset_bitmask

	return sa_code

configurations = ["updownless", "sa", "msc", "coordinated", "nodereheight", "full"]
n_reps = 100

# avoid dendropy memory leak by calling self
if len(sys.argv) == 1:
	for config in configurations:
		for rep_i in range(n_reps):
			rep_base_filename = "%s.%03d" % (config, rep_i)
			rep_base_path = os.path.join(config, rep_base_filename)
			print(rep_base_path)
			subprocess.check_call(["python3", "read-trees.py", rep_base_path])
else:
	trees_path = sys.argv[1] + ".trees"
	log_path = sys.argv[1] + ".log"
	output_path = sys.argv[1] + ".csv"

	topologies_path = "sa-topologies.newick"
	topologies = dendropy.TreeList.get(path = topologies_path, schema = "newick", rooting = "default-rooted")

	topology_sa_codes = []
	for topology in topologies:
		sa_code = calculate_sa_code(topology)
		topology_sa_codes.append(sa_code)

	trees = dendropy.TreeList.get(path = trees_path, schema = "nexus", rooting = "default-rooted", taxon_namespace = topologies.taxon_namespace)

	log_file = open(log_path)
	l = log_file.readline()
	while l[0] == "#":
		l = log_file.readline()

	origin_index = l.strip().split("\t").index("originFBD.t:Species")

	log_reader = csv.reader(log_file, dialect = csv.excel_tab)

	output_file = open(output_path, "w")
	output_writer = csv.writer(output_file)

	topology_counts = [0] * len(topologies)
	for row_i, row in enumerate(log_reader):
		output_row = []

		tree = trees[row_i]
		sa_code = calculate_sa_code(tree)
		for topology_i, topology in enumerate(topologies):
			if sa_code == topology_sa_codes[topology_i]:
				rf_distance = dendropy.calculate.treecompare.symmetric_difference(tree, topology, is_bipartitions_updated = True)
				if rf_distance == 0:
					output_row.append(topology_i)

		origin = float(row[origin_index])
		output_row.append(origin)

		root_height = max(tree.calc_node_root_distances(return_leaf_distances_only = False))
		for node in tree.preorder_internal_node_iter():
			output_row.append(root_height - node.root_distance)

		output_writer.writerow(output_row)

	output_file.close()
