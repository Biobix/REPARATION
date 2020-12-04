#####################################
##	REPARARTION: Ribosome Profiling Assisted (Re-) Annotation of Bacterial genomes
##
##	Copyright (C) 2017 Elvis Ndah
##
##	This program is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##	
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##	contact: elvis.ndah@gmail.com
#####################################


import sys
import re


"""

Generate Ribo-seq occupancy from sam alignment file
output file format
chromosome,start,strand,reads

This script using the p-site assignments based on the work from https://doi.org/10.1186/s12864-016-3278-x and https://doi.org/10.1016/j.celrep.2015.03.014
- occupancy = 1 (option -p = 1): the p-site assignment is based on the paper https://doi.org/10.1186/s12864-016-3278-x
- occupancy = 3 or 5 (option -p = 3 or -p = 5): the p-site assignment is based on the 3 or 5 prime end of the ribo-seq read based on the paper https://doi.org/10.1016/j.celrep.2015.03.014

The output file 


"""

def split_cigar(cigar):
	# function to parse CIGAR string (3M1I3M1D5M) from sam file
    # Seperate CIGAR pattern into a list of list [[code,length], .... [code,length]]
	
    cigar_pat = re.compile(r"\d+[MIDNSHP=X]{1}")
    cigar_list = []

    #split cigar entries
    for centry in cigar_pat.findall(cigar):
        count  = int(centry[:-1])
        symbol = centry[-1]
        cigar_list.append([symbol,count])

    return cigar_list


def process_ribo(inputFile,occupancy,min_read_len,max_read_len,outputFileS,outputFileAS,occupancyFile,prefix_S,prefix_AS,inFile_offset):

    #global count_reads, count_unmapped
    count_reads = 0
    count_unmapped = 0

    flags = {}
    read_table = {}     # dictionary to store ocupancy
    offset = {}

    if occupancy == 1:
        inFile = open(inFile_offset, 'r')
        line = inFile.readline()
        while line != '':
            if not (line.startswith('#') or line.startswith('length') or line.startswith('default')):
                fields = line.split()
                offset[str(fields[0])] = int(fields[1])
                
            line = inFile.readline()
        # end offset file processing

    inFile = open(inputFile, 'r')
    line = inFile.readline()
    while line != '':
        if not (line.startswith('@')):
            fields = line.split()
            read_name = str(fields[0])   # Read name;
            bit_flag = str(fields[1])   # FLAG;
            chromosome = str(fields[2])   # chromosome
            genomic_position = int(fields[3])   # left-most position
            col5 = str(fields[4])   # Map Quality
            cigar = str(fields[5])   # CIGAR string [only M and S observed in bacteria]
            col7 = str(fields[6])   # Reference sequence name of the primary alignment of the NEXT read in the template
            col8 = str(fields[7])   # Position of the primary alignment of the NEXT read in the template
            col9 = str(fields[8])   # signed observed Template LENgth
            col10 = str(fields[9])  # footprint sequence

            length = len(col10)     # footprint length

            if (int(bit_flag) == 4):	# skip unmapped read
                count_unmapped += 1

            elif (int(bit_flag) == 0):	# for plus strand

                five_prime_end = genomic_position
                three_prime_end = five_prime_end + length - 1

                # check for soft clipping in CIGAR
                if "S" in cigar:
                    cigar_list = split_cigar(cigar)
                    flag = 0    # flag to check whether to trim at 5' or 3'
                    clip_3 = 0  # flag to check if there has been any clipping at the 3'end
                    for operation in cigar_list:
                        symbol = operation[0]
                        count = operation[1]

                        if symbol == "M":
                            flag = 1

                        if flag == 0 and symbol == "S":   # if first Soft clipping operation then clip at the 5'
                            five_prime_end = genomic_position + count
                            length = length - count

                        elif flag == 1 and symbol == "S":   # if second soft clipping operation then clip at the 3'
                            length = length - count
                            three_prime_end = five_prime_end + length - 1
                            clip_3 = 1

                    if clip_3 == 0:
                        three_prime_end = five_prime_end + length - 1

                count_reads += 1
                if int(min_read_len) <= length <= int(max_read_len):

                    if occupancy == 3:
                        ribo_position = three_prime_end
                    elif occupancy == 5:
                        ribo_position = five_prime_end
                    elif occupancy == 1:
                        ribo_position = five_prime_end + offset[str(length)]

                    #print str(five_prime_end) + " processed " + str(ribo_position)

                    strand = '+'
                    if read_table.has_key(chromosome):
                        if read_table[chromosome].has_key(ribo_position):
                            if read_table[chromosome][ribo_position].has_key(strand):
                                read_table[chromosome][ribo_position][strand] += 1
                            else:
                                read_table[chromosome][ribo_position][strand] = 1
                        else:
                            read_table[chromosome][ribo_position] = {}
                            read_table[chromosome][ribo_position][strand] = 1
                    else:
                        read_table[chromosome] = {}
                        read_table[chromosome][ribo_position] = {}
                        read_table[chromosome][ribo_position][strand] = 1

            elif (int(bit_flag) == 16):		# for reverse strand

                three_prime_end = genomic_position         # 3' end is POS in sam file for the reverse strand
                five_prime_end = three_prime_end + length - 1

                # check for soft clipping in CIGAR
                if "S" in cigar:
                    cigar_list = split_cigar(cigar)

                    flag = 0    # flag to check whether to trim at 5' or 3'
                    clip_5 = 0  # flag to check if there has been any clipping at the 5'end
                    for operation in cigar_list:
                        symbol = operation[0]
                        count = operation[1]

                        if symbol == "M":
                            flag = 1

                        if flag == 0 and symbol == "S":   # if first Soft clipping operation then clip at the 3'
                            three_prime_end = genomic_position + count
                            length = length - count

                        if flag == 1 and symbol == "S":   # if second soft clipping operation then clip at the 5'
                            length = length - count
                            five_prime_end = three_prime_end + length - 1
                            clip_5 = 1

                    if clip_5 == 0: # no soft clipping has been performed at the 5' end
                        five_prime_end = three_prime_end + length - 1

                count_reads += 1
                if int(min_read_len) <= length <= int(max_read_len):
                    if occupancy == 3:
                        ribo_position = three_prime_end
                    elif occupancy == 5:
                        ribo_position = five_prime_end
                    elif occupancy == 1:
                        ribo_position = five_prime_end - offset[str(length)]

                    strand = '-'
                    if read_table.has_key(chromosome):
                        if read_table[chromosome].has_key(ribo_position):
                            if read_table[chromosome][ribo_position].has_key(strand):
                                read_table[chromosome][ribo_position][strand] += 1
                            else:
                                read_table[chromosome][ribo_position][strand] = 1
                        else:
                            read_table[chromosome][ribo_position] = {}
                            read_table[chromosome][ribo_position][strand] = 1
                    else:
                        read_table[chromosome] = {}
                        read_table[chromosome][ribo_position] = {}
                        read_table[chromosome][ribo_position][strand] = 1

            else:		# for reverse strand
                count_unmapped += 1

        line = inFile.readline()

    print "Total number of mappable reads " + str(count_reads)
    print "Total number of unmapped reads " + str(count_unmapped)

	# create bedgraph Sense
    headerS = 'track type=bedGraph name=\"' + prefix_S + '" description="" visibility=full color=0,0,255 priority=20'
    outFileS = open(outputFileS, 'w')
    outFileS.write(str(headerS) + '\n')

	# create bedgraph Antisense
    headerAS = 'track type=bedGraph name="' + prefix_AS + '" description="" visibility=full color=239,61,14 priority=20'
    outFileAS = open(outputFileAS, 'w')
    outFileAS.write(str(headerAS) + '\n')

	# create occupancy file
    table = open(occupancyFile, 'w')
    table.write('mappable_reads:'+ str(count_reads) + '\n')
    for contig in read_table:
        for start in sorted(read_table[contig].iterkeys()):
            for strand in read_table[contig][start]:
                line = [contig ,str(start), str(strand),str(read_table[contig][start][strand])]
                table.write('\t'.join(line) + '\n')

					# Bedgraph output files
                if str(strand) == '+':
                    bedgraph = [contig ,str(start - 1), str(start),str(read_table[contig][start][strand])]    # Bedgraph is 0-based
                    outFileS.write('\t'.join(bedgraph) + '\n')
                else:
                    bedgraph = [contig ,str(start - 1), str(start),str(read_table[contig][start][strand])]    # Bedgraph is 0-based
                    outFileAS.write('\t'.join(bedgraph) + '\n')


# Run python script
if __name__=='__main__':

	# Mono sample
    inputFile = sys.argv[1]				# sam file
    occupancy = int(sys.argv[2])		# plastid (1), 3' (3) or 5' (5) mapping
    min_read_len = int(sys.argv[3])		# minimun read length allow
    max_read_len = int(sys.argv[4])		# maximum red length
    outputFileS = str(sys.argv[5])		# Bedgraph Sense filename
    outputFileAS = str(sys.argv[6])		# Bedgraph AntiSense filename
    occupancyFile = str(sys.argv[7])	# RPF occupancy file
    prefix_S = str(sys.argv[8])			# prefix sense bedgraph file
    prefix_AS = str(sys.argv[9])	    # prefix antisense bedgraph file
    inFile_offset = str(sys.argv[10])	# plastid estimated offset files

    process_ribo(inputFile,occupancy,min_read_len,max_read_len,outputFileS,outputFileAS,occupancyFile,prefix_S,prefix_AS,inFile_offset)



