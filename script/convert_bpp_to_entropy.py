import math
import sys

file = sys.argv[1]
stem = {}
entro = {}
with open(file) as f:
	for line in f.readlines():
		if line == "" or line[0] == "#":	continue
		contents = line.rstrip('\n').split('\t')[1:]
		contents[2] = float(contents[2])
		stem[contents[0]] = stem.get(contents[0], 0.0)+contents[2]
		stem[contents[1]] = stem.get(contents[1], 0.0)+contents[2]
		if contents[2] == 0.0:
			ent = 0.0
		else:
			ent = contents[2]*math.log2(contents[2])
		entro[contents[0]] = entro.get(contents[0], 0.0)-ent
		entro[contents[1]] = entro.get(contents[1], 0.0)-ent

for key_num in sorted(list(map(int, entro.keys()))):
	key = str(key_num)
	if stem[key] > 0.0:
		entro[key] -= (1.0-stem[key])*math.log2(1.0-stem[key])
	print("\t".join(["*", key, str(entro[key])]))
