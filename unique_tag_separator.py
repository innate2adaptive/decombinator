f = open('/home/tom/complex/MP1/bashtest/dcr_beta_VT8_N712_S502_R2_forward.n12', 'r')

v_tags = {}
j_tags= {}
for line in f:
	dcr = line.rstrip().split(", ")

	if dcr[0] != "n/a" and dcr[1] == "n/a":
		if dcr[0] in v_tags:
			v_tags[dcr[0]].append(line)
		else:
			v_tags[dcr[0]] = []
			v_tags[dcr[0]].append(line)

	elif dcr[0] == "n/a" and dcr[1] != "n/a":
		if dcr[1] in j_tags:
			j_tags[dcr[1]].append(line)
		else:
			j_tags[dcr[1]] = []
			j_tags[dcr[1]].append(line)
f.close()

import pprint
pp = pprint.PrettyPrinter(indent=2)
print "v_tags"
pp.pprint(v_tags)
print "j_tags:"
pp.pprint(j_tags)
