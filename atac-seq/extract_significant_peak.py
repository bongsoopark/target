import sys

# peak file
f1 = sys.argv[1]

# merged significant peak file
f2 = sys.argv[2]

f = open(f2, "r")
x = {}
for line in f:
	line = line.strip()
	data = line.split("\t")
	the_key = data[0] + ":" + data[1] + "-" + data[2]
	x.update({the_key:1})
f.close()

f = open(f1, "r")
for line in f:
	line = line.strip()
	data = line.split("\t")
	the_key = data[0] + ":" + data[1] + "-" + data[2]
	if x.has_key(the_key):
		print line+"\tReproducible IDR p-value < 0.05"
f.close()
