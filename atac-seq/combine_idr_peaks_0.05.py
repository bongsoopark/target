# Peak calling reproducibility test
# Input files: Three narrowPeak from IDR
# Output files:
# Johns Hopkins Environmental Health and Engineering
# Bongsoo Park, Ph.D
# Last update: 2018-05-05

# Import libraries
import sys
import operator

# Retrieve three peak files (narrowPeak format from IDR)
# file #1
f1 = sys.argv[1]

# file #2
f2 = sys.argv[2]

# file #3
f3 = sys.argv[3]

# option (merging option)
f4 = "mean"

print "==========================================="
print "combine IDR tested narrowPeaks"
print "peak1 vs 2", f1
print "peak2 vs 3", f2
print "peak3 vs 1", f3
# narrowPeak1,2,3: save peak information
# narrowPeakSize1,2,3: save peak position
narrowPeak = {}
narrowPeak1 = {}
narrowPeakSize1 = {}
narrowPeakEveryPosition = {}
f = open(f1, "r")
for line in f:
	line = line.strip()
	data = line.split("\t")
	# 0.05 have a score of int(-125log2(0.05)) = 540
	# IDR score > 540: cutoff 0.05
	if int(data[4]) > 540:
		the_key = data[0]+":"+data[1]+"-"+data[2]
		the_key2 = data[0]+":"+data[16]+"-"+data[17]
		narrowPeak1.update({the_key:int(data[9])})
		narrowPeak.update({the_key:int(data[9])})
		narrowPeak.update({the_key2:int(data[19])})
		for i in range(int(data[1]), int(data[2])+1):
			# Save one nucleotide based peak position
			the_position = data[0]+":"+str(i)
			narrowPeakSize1.update({the_position:1})	
			narrowPeakEveryPosition.update({the_position:1})
f.close()
print "==========================================="
print "Total number of peak1:", len(narrowPeak1)
print "Total size of peak1:", len(narrowPeakSize1)
print "Total number of peaks:", len(narrowPeak)
print "==========================================="

narrowPeak2 = {}
narrowPeakSize2 = {}
f = open(f2, "r")
for line in f:
	line = line.strip()
	data = line.split("\t")
	# 0.05 have a score of int(-125log2(0.05)) = 540
	# IDR score > 540: cutoff 0.05
	if int(data[4]) > 540:
		the_key = data[0]+":"+data[1]+"-"+data[2]
		the_key2 = data[0]+":"+data[16]+"-"+data[17]
		narrowPeak2.update({the_key:int(data[9])})
		narrowPeak.update({the_key:int(data[9])})
		narrowPeak.update({the_key2:int(data[19])})
		for i in range(int(data[1]), int(data[2])+1):
			# Save one nucleotide based peak position
			the_position = data[0]+":"+str(i)
			narrowPeakSize2.update({the_position:1})	
			if narrowPeakEveryPosition.has_key(the_key):
				narrowPeakEveryPosition[the_key] + 1
			else:
				narrowPeakEveryPosition.update({the_position:1})
			
f.close()
print "Total number of peak2:", len(narrowPeak2)
print "Total size of peak2:", len(narrowPeakSize2)
print "Total number of peaks:", len(narrowPeak)
print "==========================================="

narrowPeak3 = {}
narrowPeakSize3 = {}
f = open(f3, "r")
for line in f:
	line = line.strip()
	data = line.split("\t")
	# 0.05 have a score of int(-125log2(0.05)) = 540
	# IDR score > 540: cutoff 0.05
	if int(data[4]) > 540:
		the_key = data[0]+":"+data[1]+"-"+data[2]
		the_key2 = data[0]+":"+data[16]+"-"+data[17]
		narrowPeak3.update({the_key:int(data[9])})
		narrowPeak.update({the_key:int(data[9])})
		narrowPeak.update({the_key2:int(data[19])})
		for i in range(int(data[1]), int(data[2])+1):
			# Save one nucleotide based peak position
			the_position = data[0]+":"+str(i)
			narrowPeakSize3.update({the_position:1})	
			if narrowPeakEveryPosition.has_key(the_key):
				narrowPeakEveryPosition[the_key] + 1
			else:
				narrowPeakEveryPosition.update({the_position:1})
f.close()
print "Total number of peak3:", len(narrowPeak3)
print "Total size of peak3:", len(narrowPeakSize3)
print "Total number of peaks:", len(narrowPeak)
print "==========================================="
print "narrowPeak", len(narrowPeak)
print "narrowPeak EveryPosition", len(narrowPeakEveryPosition)

f = open("merged.narrowPeak","w")
sortedPeak = sorted(narrowPeak.items(), key=operator.itemgetter(1), reverse=True)
for the_key, the_value in sortedPeak:
	coor = the_key.split(":")
	chrom = coor[0]
	start_pos = coor[1].split("-")[0]
	end_pos = coor[1].split("-")[1]
	f.write(chrom+"\t"+start_pos+"\t"+end_pos+"\t.\t"+str(the_value)+"\n")
f.close()
	
