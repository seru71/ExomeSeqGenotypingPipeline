#!/usr/bin/env python
#
#
# The script chops BED intervals in such a way that for each overlapping pair of intervals  
# in the output reciprocal overlap is 100%. Example:
# 
# INPUT:
# A	----------
# B	     ------------
# C                         ------------------
# D                              ---
# OUTPUT:
# A1    -----
# A2         -----
# B1         -----
# B2              -------
# C1                        -----
# C2				 ---
# C3	                            ----------
# D                              ---
#
# Input BED must be sorted by <start> coordinate. Directonality (+/-) is not taken into account.
# As in BED format, intervals are inclusive <start> position and exclusive <end> position, i.e. [start,end)
#
# by Pawel Sztromwasser
# Feb, 2016
#

def get_exon(input_line):
	line_split = input_line.strip().split('\t')
	return {'chrom': line_split[0], 
                'start': int(line_split[1]), 
                'end': int(line_split[2]), 
		'rest': line_split[3:]}


def chop_exons(exons, end_position):
	""" Input is sorted by start """
	
#	print 'chopping '+str(exons)	


	sorted_ends = sorted([e['end'] for e in exons])
	left = exons[0]['start']
	while left <= end_position:

# 		print 'wheee:', left, end_position, sorted_ends

		# we will surely chop the first one
		chop_off = [0]
		
		# how many additional start at the same position?
		for i in range(1, len(exons)):
			if exons[i]['start'] > left:
				break
			chop_off.append(i)
		
		# find right-hand side of the intervals to chop off
		min_end = sorted_ends[0]
		min_start = exons[chop_off[-1]+1]['start'] if len(chop_off) < len(exons) else min_end+1
		
		if min_start < min_end:
			right = min_start
		else:
			right = min_end
			while len(sorted_ends)>0 and sorted_ends[0] == min_end:
				sorted_ends.pop(0)			

		# print chopped intervals and update their start coordinate
		for i in chop_off:
			print '\t'.join([exons[i]['chrom'],
					str(left),
					str(right),
					'\t'.join(exons[i]['rest'])])
			exons[i]['start'] = right

		# update exons list		
		i = 0
		while i < len(exons):
			if exons[i]['start'] >= exons[i]['end']:	exons.pop(i)
			else: i+=1

		# move the left coordinate
		if len(exons) > 0:  # in case we popped everyone
			left = right  #=exons[0]['start']
		else:
			break



import sys

buf = []
buf_chrom = None
buf_end_position = None

in_bed = open(sys.argv[1])
for l in in_bed.xreadlines():

	exon = get_exon(l)

	if buf_chrom == None: 
		buf_chrom=exon['chrom']
		buf_end_position = exon['end']

	# accumulate overlapping exons in buf
	if exon['start'] <= buf_end_position and exon['chrom'] == buf_chrom:
		buf.append(exon)
		if exon['end'] > buf_end_position:
			buf_end_position = exon['end']
	# and if next is not overlapping...
	else:
		chop_exons(buf, buf_end_position)
		buf = [exon]
		buf_chrom = exon['chrom']
		buf_end_position = exon['end']

in_bed.close()

# chop last batch
if buf != []: 
	chop_exons(buf, buf_end_position)
		



