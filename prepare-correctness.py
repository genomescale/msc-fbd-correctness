import os
import numpy
import xml.etree.ElementTree as ET

configurations = ["updownless", "sa", "msc", "coordinated", "nodereheight", "full"]

for configuration in configurations:
	template_fn = configuration + ".xml"
	template_path = os.path.join(configuration, template_fn)
	xml_tree = ET.parse(template_path)
	xml_root = xml_tree.getroot()

	n_replicates = 100

	numpy.random.seed()

	master_fn = "kickoff.sh"
	master_path = os.path.join(configuration, master_fn)
	master_file = open(master_path, "w")

	beast_cmd = "~/beast/bin/beast"

	for rep_i in range(n_replicates):
		rep_name = "%s.%03d" % (configuration, rep_i)
		for logger_element in xml_root.findall(".//logger"):
			if "fileName" in logger_element.attrib:
				original_fn = logger_element.attrib["fileName"]
				original_ext = original_fn[original_fn.rfind("."):]
				new_fn = rep_name + original_ext
				logger_element.attrib["fileName"] = new_fn

		output_fn = rep_name + ".xml"
		output_path = os.path.join(configuration, output_fn)
		xml_tree.write(output_path)

		random_seed = numpy.random.randint(1, 2**32 - 1)

		beast_args = ["-overwrite", "-seed", str(random_seed), "-statefile", rep_name + ".state", "-threads", "1", output_fn]
		master_file.write(beast_cmd + " " + " ".join(beast_args) + "\n")

	master_file.close()
