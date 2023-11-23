import os
import matplotlib
matplotlib.use('Agg') # Configure matplotlib to use 'Agg' backend for non-GUI environments
from matplotlib import pyplot as plt 

shared={} # Dictionary to store shared genomic segments
firstgen=['P1','P2','P3','P4'] # List of first generation individuals

# Function to check if two segments overlap
def overlap(seg1, seg2):
    if (seg1[0]>=seg2[0] and seg1[0]<=seg2[1]) or (seg2[0]>=seg1[0] and seg2[0]<=seg1[1]):
        return True
    else:
        return False

# Read centromere track data from a file
infile=open("/net/harris/vol1/data/hg38/ucsc_hg38_centromere_tracks.txt")
lines=infile.readlines()
infile.close()

centromere_start, centromere_end={},{}
for chrom in range(1,23):
    centromere_start[chrom]=10**9 # Initialize centromere start coordinates (this is larger than any centromere start coordinate will realistically be)                                           
    centromere_end[chrom]=0 # Initialize centromere end coordinates

# Process centromere track data
for line in lines[1:]:
    s = line.strip('\n').split('\t')
    if s[0].startswith('chr'):
        chrom = int(s[0][3:])  # Extract chromosome number without "chr" prefix
        # Update centromere start and end coordinates
        if int(s[1]) < centromere_start[chrom]:
            centromere_start[chrom] = int(s[1])
        if int(s[2]) > centromere_end[chrom]:
            centromere_end[chrom] = int(s[2])

# Similarly, read and process telomere track data
telomere_start, telomere_end = {}, {}
telomere_file = open("/net/harris/vol1/data/hg38/ucsc_hg38_telomere_tracks.txt")
telomere_lines = telomere_file.readlines()
telomere_file.close()

for line in telomere_lines[1:]:
    s = line.strip('\n').split('\t')
    if s[0].startswith('chr'):
        chrom = int(s[0][3:])  # Extract chromosome number without "chr" prefix
        # Update telomere start and end coordinates
        telomere_start[chrom] = int(s[1])
        telomere_end[chrom] = int(s[2])

output={} # Dictionary for output data

# Process each chromosome
for chrom in range(1,23):
    # Function to simplify segments by merging overlapping segments and removing small gaps    
    def simplify(segs):
        cen1, cen2 = centromere_start[chrom], centromere_end[chrom] #coordinates of centromere
        tel1, tel2 = telomere_start[chrom], telomere_end[chrom]  #coordinates of telomere
        segs.sort() # Sort segments
	# Merge overlapping segments and remove small gaps
        ind = 1
        while ind < len(segs):  # get rid of small gaps between segments
            if overlap([segs[ind - 1][0], segs[ind - 1][1] + 10 ** 7], segs[ind]):
                segs[ind - 1] = [min(segs[ind - 1][0], segs[ind][0]), max(segs[ind - 1][1], segs[ind][1])]
                segs.pop(ind)
            else:
                ind += 1                
        
        ind=0 #remove centromere

	# Remove segments overlapping with centromere or telomere        
        while ind<len(segs):
            if overlap(segs[ind], [cen1, cen2]) or overlap(segs[ind], [tel1, tel2]):
	    # Logic for adjusting segment boundaries based on centromere and telomere overlap
                if cen1<=segs[ind][0]:
                    if cen2>=segs[ind][1]:
                        segs.pop(ind)
                    else:
                        segs[ind][0]=cen2
                else:
                    if cen2>=segs[ind][1]:
                        segs[ind][1]=cen1
                    else:
                        segs.insert(ind,[segs[ind][0],cen1])
                        segs[ind+1][0]=cen2
            ind+=1
        return segs
    
    # Initialize shared dictionary for first generation individuals
    for i in range(len(firstgen)):
        for h in 'ab':
            for j in range(i+1,len(firstgen)):
                for h2 in 'ab':
                    shared[(firstgen[i]+h,firstgen[j]+h2)]=[]
                
    # Read IBD (Identity by Descent) file for the current chromosome
    infile = open("/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/surrogate_accessible_regions/IBD_files/chr" + str(chrom) + "_reformat.ibd")
    lines=infile.readlines()
    infile.close()

    # Process IBD data
    for line in lines[1:]:
        # Logic for processing IBD data and updating shared dictionary
        s=line.strip('\n').split('\t')
        if s[-1][:2] in firstgen and s[-2][:2] in firstgen and not s[-1][:2]==s[-2][:2]:
            if s[-2]<s[-1]:
                shared[(s[-2],s[-1])].append([int(float(s[1])),int(float(s[2]))])
            else:
                shared[(s[-1],s[-2])].append([int(float(s[1])),int(float(s[2]))])

    for sibkey in shared.keys():
        shared[sibkey]=simplify(shared[sibkey])
        
    overlap_class={}

    for sibkey1 in shared.keys():
        for ind1 in range(len(shared[sibkey1])):
            overlap_class[(sibkey1,ind1)]=[]
    for sibkey1 in shared.keys():
        for sibkey2 in shared.keys():
            if not sibkey1==sibkey2:
                for ind1 in range(len(shared[sibkey1])):
                    for ind2 in range(len(shared[sibkey2])):
                        if overlap(shared[sibkey1][ind1], shared[sibkey2][ind2]):
                            overlap_class[(sibkey1,ind1)].append((sibkey2,ind2))
                            overlap_class[(sibkey2,ind2)].append((sibkey1,ind1))

    trio_segs={}
    duo_segs={}

    # Generate output for each first generation individual
    for i in range(len(firstgen)):
        for j in range(i+1,len(firstgen)):
            duo_segs[(firstgen[i],firstgen[j])]=[]
            for k in range(j+1,len(firstgen)):
                trio_segs[(firstgen[i],firstgen[j],firstgen[k])]=[]
                trio_segs[(firstgen[j],firstgen[i],firstgen[k])]=[]
                trio_segs[(firstgen[k],firstgen[i],firstgen[j])]=[]

    ab='ab'
    for i in range(len(firstgen)):
        for j in range(i+1,len(firstgen)):
            for h2 in range(2):
                pair1=(firstgen[i]+'a',firstgen[j]+ab[h2])
                for ind1 in range(len(shared[pair1])):
                    pair2=(firstgen[i]+'b',firstgen[j]+ab[1-h2])
                    for ind2 in range(len(shared[pair2])):
                        if overlap(shared[pair1][ind1],shared[pair2][ind2]):
                            duo_segs[(firstgen[i],firstgen[j])].append([max(shared[pair1][ind1][0],shared[pair2][ind2][0]),min(shared[pair1][ind1][1],shared[pair2][ind2][1])])
                for k in range(j+1,len(firstgen)):
                    for h3 in range(2):
                        pair2=(firstgen[i]+'b',firstgen[k]+ab[h3])
                        for ind1 in range(len(shared[pair1])):
                            for ind2 in range(len(shared[pair2])):
                                if overlap(shared[pair1][ind1],shared[pair2][ind2]):
                                    trio_segs[(firstgen[i],firstgen[j],firstgen[k])].append([max(shared[pair1][ind1][0],shared[pair2][ind2][0]),min(shared[pair1][ind1][1],shared[pair2][ind2][1])])
                pair1=(firstgen[i]+ab[h2],firstgen[j]+'a')
                for k in range(j+1,len(firstgen)):
                    for h3 in range(2):
                        pair2=(firstgen[j]+'b',firstgen[k]+ab[h3])
                        for ind1 in range(len(shared[pair1])):
                            for ind2 in range(len(shared[pair2])):
                                if overlap(shared[pair1][ind1],shared[pair2][ind2]):
                                    trio_segs[(firstgen[j],firstgen[i],firstgen[k])].append([max(shared[pair1][ind1][0],shared[pair2][ind2][0]),min(shared[pair1][ind1][1],shared[pair2][ind2][1])])
                        pair1=(firstgen[i]+ab[h2],firstgen[k]+'a')
                        pair2=(firstgen[j]+ab[h3],firstgen[k]+'b')
                        for ind1 in range(len(shared[pair1])):
                            for ind2 in range(len(shared[pair2])):
                                if overlap(shared[pair1][ind1],shared[pair2][ind2]):
                                    trio_segs[(firstgen[k],firstgen[i],firstgen[j])].append([max(shared[pair1][ind1][0],shared[pair2][ind2][0]),min(shared[pair1][ind1][1],shared[pair2][ind2][1])])

    if chrom==1:
        for thruple in trio_segs.keys():
            output[thruple]=""
        for pair in duo_segs.keys():
            output[pair]=""
                                    
    for thruple in trio_segs.keys():
        trio_segs[thruple]=simplify(trio_segs[thruple])
        for i in range(len(trio_segs[thruple])):
            if trio_segs[thruple][i][0]>0:
                output[thruple]+='chr'+str(chrom)+'\t'+str(trio_segs[thruple][i][0]-1)+'\t'+str(trio_segs[thruple][i][1])+'\n'
            else:
                output[thruple]+='chr'+str(chrom)+'\t'+str(trio_segs[thruple][i][0])+'\t'+str(trio_segs[thruple][i][1])+'\n'
        
    for pair in duo_segs.keys():
        duo_segs[pair]=simplify(duo_segs[pair])
        for i in range(len(duo_segs[pair])):
            if duo_segs[pair][i][0]>0:
                output[pair]+="chr"+str(chrom)+"\t"+str(duo_segs[pair][i][0]-1)+'\t'+str(duo_segs[pair][i][1])+'\n'
            else:
                output[pair]+="chr"+str(chrom)+"\t"+str(duo_segs[pair][i][0])+'\t'+str(duo_segs[pair][i][1])+'\n'

output_directory = '/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/surrogate_accessible_regions/beds'

for thruple in trio_segs.keys():
    output_file = os.path.join(output_directory, "proband_gen_positive_mask", thruple[0] + "_v_surrogate_parents" + thruple[1] + "_" + thruple[2] + ".bed")
    with open(output_file, "w") as outfile:
        outfile.write(output[thruple])

for pair in duo_segs.keys():
    output_file = os.path.join(output_directory, "proband_gen_positive_mask", pair[0] + "_" + pair[1] + "_double_surrogates.bed")
    with open(output_file, "w") as outfile:
        outfile.write(output[pair])
