# Generate reproducible sub peaks from the MACS result
# Bongsoo Park, Johns Hopkins
# Import libraries
import sys
import os

# Import three biological replicates
f1 = sys.argv[1]
f2 = sys.argv[2]
f3 = sys.argv[3]
print "# Import three biological replicates"
print "file1 : ",f1
print "file2 : ",f2
print "file3 : ",f3

# Conduct idr test
print "# Conduct IDR test"
print "# IDR for sample 1 vs 2"
os.system("idr --samples "+f1+" "+f2+" --output-file peak12.tmp.narrowPeak")
print "# IDR for sample 2 vs 3"
os.system("idr --samples "+f2+" "+f3+" --output-file peak23.tmp.narrowPeak")
print "# IDR for sample 3 vs 1"
os.system("idr --samples "+f3+" "+f1+" --output-file peak31.tmp.narrowPeak")

# Combine reproducible peaks
print "Combine reproducible peaks IDR < 0.05"
os.system("python combine_idr_peaks_0.05.py peak12.tmp.narrowPeak peak23.tmp.narrowPeak peak31.tmp.narrowPeak")

# Extract significant peaks from the orignal peak file
print "Extract significant peaks from the original peak file (creating .reproducible"
os.system("python extract_significant_peak.py "+f1+" merged.narrowPeak > "+f1+".reproducible")  
os.system("python extract_significant_peak.py "+f2+" merged.narrowPeak > "+f2+".reproducible")  
os.system("python extract_significant_peak.py "+f3+" merged.narrowPeak > "+f3+".reproducible")  
# End of the script
print "End..."
