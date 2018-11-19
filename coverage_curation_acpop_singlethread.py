#!/usr/bin/env python

import gzip, glob, time, sys, os, numpy
from itertools import *

#reading the line and splitting it up properly
def readline_from_gz(filehandle):
    return filehandle.readline().strip().split("\t")

#No need to read any line or process the data if not using it
#but we do need to move the iterator ahead one step in the file
#so we can read next line
def seekforward(filehandle):
    return next(filehandle)

#I got this from https://stackoverflow.com/questions/29737195/python-jump-to-a-line-in-a-txt-file-a-gzipped-one
#Looks like just seeking forward with next is faster though.
def consume(iterator, n=None):
    '''Advance the iterator n-steps ahead. If n is None, consume entirely.'''
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)

#do some stats
def calc_stats(coverages):
    num_sample = float(len(coverages))
    results = {1:0, 5:0, 10:0, 15:0, 20:0, 25:0, 30:0, 50:0, 100:0}
    for s in coverages:
        for k in results.keys():
            if s>k:
                results[k]+=1
    for k in results.keys():
        results[k] = round(results[k]/num_sample, 4)
    return results



t_zero = time.time()
coverage_files = glob.glob("*.gz")

filedata = [gzip.open(filename, 'r') for filename in coverage_files]

if len(sys.argv) != 2:
    exit("usage: ACpop_curate_coverage_single.py <output_filename>")
part = 0

#output_file = open(sys.argv[1]+".part"+str(part), "w")
output_file_subsample = open(sys.argv[1]+".part"+str(part)+".subsampled", "w")

i=0
num_sample = len(filedata)

#output_file.write("#chrom\tpos\tmean\tmedian\t1\t5\t10\t15\t20\t25\t30\t50\t100\n")
output_file_subsample.write("#chrom\tpos\tmean\tmedian\t1\t5\t10\t15\t20\t25\t30\t50\t100\n")
t_one = time.time()
t_two = time.time()
file_loc = {}
for name in coverage_files:
    file_loc[name] = 0

while True:
    #reading the first line to check if we can skip reading the 299 other files
    firstline = readline_from_gz(filedata[0])
    chromosome = firstline[0]
    position = firstline[1]

    if int(position) % 10 != 0:
        #move forward in the file iterator one step - skip a line if not needed
        for handle in filedata[1:]:
            seekforward(handle)
        continue

    #if we are not skipping the line, then lets read it from the rest of the files
    #I wish I could parallelise this step but I cant seem to get that to work...
    linedata = [readline_from_gz(handle) for handle in filedata[1:]]

    coverages = [int(firstline[2])] #Adding the first coverage value

    for sample in linedata:
        if sample:
            try:
                coverages.append(int(sample[2]))
            except ValueError as e: #when coverage is very high it is listed as X.XXXXeYY which cant be parsed by int()
                try:
                    coverages.append(int(float(sample[2]))) #turns out float can parse that
                except ValueError as e:
                    coverages.append(0)
                    print "BROKEN COVERAGE VALUE", sample
        else:
            print "No more file"
            exit()
        if sample[1] != position:
            print "broken pos"
            exit()
    coverage_mean = round(numpy.mean(coverages), 4)
    coverage_median = numpy.median(coverages)
    xcovs = calc_stats(coverages)

    #output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
    #	chromosome, position, coverage_mean, coverage_median,
    #	xcovs[1], xcovs[5], xcovs[10], xcovs[15], xcovs[20],
    #	xcovs[25], xcovs[30], xcovs[50], xcovs[100]))
    if int(position) % 10 == 0:
        output_file_subsample.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            chromosome, position, coverage_mean, coverage_median,
            xcovs[1], xcovs[5], xcovs[10], xcovs[15], xcovs[20],
            xcovs[25], xcovs[30], xcovs[50], xcovs[100]))

    #just cutting up the output files a little
    i+=1
    if i % 1000000 == 0:
        #output_file.close()
        output_file_subsample.close()
        part+=1
        #output_file = open(sys.argv[1]+".part"+str(part), "w")
        output_file_subsample = open(sys.argv[1]+".part"+str(part)+".subsampled", "w")
        t_one = t_two
        t_two = time.time()
        print round(t_two - t_one, 2), "s"

t_three = time.time()

print round(t_three - t_zero, 2), "s"
