# Peak calling reproducibility test
# Input files: Two narrowPeak from MACS
# Output files:
# Johns Hopkins Environmental Health and Engineering
# Bongsoo Park, Ph.D
# Last update: 2018-05-05

# Import libraries
import sys

# Retrieve two peak files (narrowPeak format from MACS)
# file #1
f1 = sys.argv[1]

# file #2
f2 = sys.argv[2]

print "==========================================="
print "peak reproducibility test, developed by Bongsoo Park. Johns Hopkins"
print "peak1", f1
print "peak2", f2
# narrowPeak1,2: save peak information
# narrowPeakSize1,2: save peak position
narrowPeak1 = {}
narrowPeakSize1 = {}
peak1_score = 0
f = open(f1, "r")
for line in f:
	line = line.strip()
	data = line.split("\t")
	the_key = data[0]+":"+data[1]+"-"+data[2]
	narrowPeak1.update({the_key:int(data[4])})
	peak1_score += int(data[4])
	for i in range(int(data[1]), int(data[2])+1):
		# Save one nucleotide based peak position
		the_position = data[0]+":"+str(i)
		narrowPeakSize1.update({the_position:1})	
f.close()
print "==========================================="
print "Total number of peak1:", len(narrowPeak1)
print "total score of peak1:", peak1_score
print "Total size of peak1:", len(narrowPeakSize1)
print "==========================================="

narrowPeak2 = {}
narrowPeakSize2 = {}
peak2_score = 0
overlap_peak_20 = 0
overlap_peak_50 = 0
overlap_peak_20_hscore = 0
overlap_peak_50_hscore = 0
overlap_pos = 0
f = open(f2, "r")
for line in f:
	line = line.strip()
	data = line.split("\t")
	the_key = data[0]+":"+data[1]+"-"+data[2]
	narrowPeak2.update({the_key:int(data[4])})
	peak2_score += int(data[4])
	flag = 0
	for i in range(int(data[1]), int(data[2])+1):
		# Save one nucleotide based peak position
		the_position = data[0]+":"+str(i)
		narrowPeakSize2.update({the_position:1})	
		if narrowPeakSize1.has_key(the_position):
			overlap_pos += 1
			flag += 1
	# if overlapped position is more than 20%, we consider these are overlapped peak
	# print flag, int(data[2]) - int(data[1])
	if flag > (int(data[2]) - int(data[1]))/5:
			overlap_peak_20 += 1 
			if int(data[4]) > 100:
				overlap_peak_20_hscore += 1
	if flag > (int(data[2]) - int(data[1]))/2:
			overlap_peak_50 += 1 
			if int(data[4]) > 100:
				overlap_peak_50_hscore += 1
		
f.close()
print "Total number of peak2:", len(narrowPeak2)
print "total score of peak2:", peak2_score
print "Total size of peak2:", len(narrowPeakSize2)
print "==========================================="
print "Overlap peak(>20%): reproducible", overlap_peak_20
print "Overlap peak(>50%): highly reproducible", overlap_peak_50
print "Overlap peak(>20%): reproducible hscore (enrichment > 100)", overlap_peak_20_hscore
print "Overlap peak(>50%): highly reproducible hscore (enrichment > 100)", overlap_peak_50_hscore
print "Overlapped position between peak1 and peak2:", overlap_pos
