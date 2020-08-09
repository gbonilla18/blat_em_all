import urllib
import urllib2
import time
import re
import xml.dom.minidom
import itertools
import datetime

def get_page(url,values):
	data = urllib.urlencode(values)
	req = urllib2.Request(url, data)
	response = urllib2.urlopen(req)
	the_page = response.read()
	return the_page

	
# just a generic function to write to a file, will use it for testing	
def write_page(web_page,output_fname):
	output_file=file(output_fname,"w")
	output_file.writelines(web_page)
	output_file.close()

def read_in_probes_file(input_fname):
	probe_list_f=open(input_fname)
	probe_list=probe_list_f.readlines()
	probe_list_f.close()
	probe_data=[]
	for l in probe_list:
		broken_down=l.strip().split("\t")
		probe_data.append(broken_down)
	
	return probe_data
	
def read_in_probes_textarea(input_str):
	probe_list=input_str.split("\n")
	probe_data=[]
	for l in probe_list:
		broken_down=l.strip().split()
		if len(broken_down)>1:
			probe_data.append(broken_down)
	#print probe_data
	return probe_data

def probe_list_to_dict(probe_data):
	original_probe_dict={}
	for p in probe_data:
		original_probe_dict[p[0]]=p[1]
	
	return original_probe_dict
	
def probe_list_to_fasta(probe_data,id_prefix):
	fasta_lines=[">"+id_prefix+p[0]+"\n"+p[1]+"\n" for p in probe_data ]
	fasta_in_one="\n".join(fasta_lines)
	return fasta_in_one

def probe_list_to_fasta_list(probe_data,id_prefix):
	fasta_lines=[">"+id_prefix+p[0]+"\n"+p[1]+"\n" for p in probe_data ]
	#fasta_in_one="\n".join(fasta_lines)
	return fasta_lines

#suffix1 and suffix2 correspond to the distance extension that is used, and they are lists of the same length as probe_data	
def probe_list_to_fasta_w_dist(probe_data,id_prefix,suffix1, suffix2):
	### make sequence 20 bp adding all possible kmers 
	fasta_lines=[]
	for p in range(0,len(probe_data)):
		fasta_lines1=[]
		fasta_lines2=[]
		cseq=probe_data[p][1]
		#print cseq
		lendiff=20-len(cseq)
		#print lendiff
		if lendiff>0:
			## append sequences of lendiff
			#lendiff=lendiff
			allkmers=getAllPossibleKmers(lendiff)
			fasta_lines1=[">"+id_prefix+probe_data[p][0]+"-seql_"+str(km)+"-dist1_"+str(suffix1[p]+lendiff)+"-dist2_"+str(suffix2[p])+"\n"+allkmers[km]+probe_data[p][1]+"\n"  for km in range(0,len(allkmers))]
			fasta_lines2=[">"+id_prefix+probe_data[p][0]+"-seqr_"+str(km)+"-dist1_"+str(suffix1[p])+"-dist2_"+str(suffix2[p]+lendiff)+"\n"+probe_data[p][1]+allkmers[km]+"\n"  for km in range(0,len(allkmers))]
		fasta_lines.extend(fasta_lines1)
		fasta_lines.extend(fasta_lines2)
	
	fasta_lines.extend([">"+id_prefix+probe_data[p][0]+"-dist1_"+str(suffix1[p])+"-dist2_"+str(suffix2[p])+"\n"+probe_data[p][1]+"\n" for p in range(0,len(probe_data)) ])
	fasta_in_one="\n".join(fasta_lines)
	
	return fasta_in_one
def probe_list_to_fasta_list_w_dist(probe_data,id_prefix,suffix1,suffix2):
	fasta_lines=[]
	for p in range(0,len(probe_data)):
		cseq=probe_data[p][1]
		fasta_lines1=[]
		fasta_lines2=[]
		#print cseq
		lendiff=20-len(cseq)
		#print lendiff
		if lendiff>0:
			## append sequences of lendiff
			#lendiff=lendiff
			allkmers=getAllPossibleKmers(lendiff)
			fasta_lines1=[">"+id_prefix+probe_data[p][0]+"-seql_"+str(km)+"-dist1_"+str(suffix1[p]+lendiff)+"-dist2_"+str(suffix2[p])+"\n"+allkmers[km]+probe_data[p][1]+"\n"  for km in range(0,len(allkmers))]
			fasta_lines2=[">"+id_prefix+probe_data[p][0]+"-seqr_"+str(km)+"-dist1_"+str(suffix1[p])+"-dist2_"+str(suffix2[p]+lendiff)+"\n"+probe_data[p][1]+allkmers[km]+"\n"  for km in range(0,len(allkmers))]
			fasta_lines.extend(fasta_lines1)
			fasta_lines.extend(fasta_lines2)
	
	fasta_lines.extend([">"+id_prefix+probe_data[p][0]+"-dist1_"+str(suffix1[p])+"-dist2_"+str(suffix2[p])+"\n"+probe_data[p][1]+"\n" for p in range(0,len(probe_data)) ])
	
	return fasta_lines
	
def get_chrom_location(loc_string):
	loc1=loc_string.split(":")
	chr=loc1[0]
	loc2=loc1[1].split("-")
	start=int(loc2[0])
	end=int(loc2[1])
	return [chr,start,end]
	
def get_alignment_info(alignment_res_list):
	alignment_info=[]
	for alignment_str in alignment_res_list:
		alignment_info_split=alignment_str.strip().split()
		# the first element in alignment_info_split is the id, and we want to get the dist1 and dist2 from there
		#print "<Br>alignment info split: "
		#print(alignment_info_split)
		#print "<Br>"
		pn_dist=re.compile("dist1_([0-9]+)-dist2_([0-9]+)")
		dist_from_id=pn_dist.findall(alignment_info_split[0])[0]
		#print dist_from_id
		
		pn_id=re.compile("^(.+?)-")
		id_wo_dist=pn_id.findall(alignment_info_split[0])[0]
		
		
		## returning id without the distance
		alignment_info_temp=[id_wo_dist,dist_from_id[0],dist_from_id[1]]

		#id with dist
		#print "<br> id with dist: " 
		#print "-".join(alignment_info_split[0].split("-")[0:3])
		id_w_dist="-".join(alignment_info_split[0].split("-")[0:3])
		#alignment_info_temp=[id_w_dist,dist_from_id[0],dist_from_id[1]]
		alignment_info_temp.extend(alignment_info_split[1:])
		alignment_info.append(alignment_info_temp)
	return alignment_info		

def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)

def get_dna_from_ucsc(chr_loc,genome):
	dna_url="http://genome.ucsc.edu/cgi-bin/das/"+genome+"/dna?segment="
	new_loc_dna_url=dna_url+chr_loc
	#print new_loc_dna_url
	new_dna_page=get_page(new_loc_dna_url,[])
			
	dom3 = xml.dom.minidom.parseString(new_dna_page)
	dna_str=getText(dom3.getElementsByTagName("DNA")[0].childNodes).strip()
	return dna_str

### this was adapted from something I got from the internet (adapted to work with our version of python)	
def ReverseComplement3(seq):
	nsseq="".join(seq.split())
	for base in nsseq:
		if base not in 'ATCGatcg':
			print "Error: NOT a DNA sequence"
			return None
	seq1 = 'ATCGTAGCatcgtagc'
	
	seq_list=[]
	for i in range(16): 
		if i < 4 or 8<=i<12 :
			seq_list.append((seq1[i],seq1[i+4]))
	seq_dict=dict(seq_list)
	
	#seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
	return "".join([seq_dict[base] for base in reversed(nsseq)])

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)


def getAllPossibleKmers(k):
	bases=['A','T','C','G']
	return [''.join(p) for p in product(bases, repeat=k)]

def do_blat_multiple(sequences,genome):
	
	probe_ids={}
	probes_in_dict={}

	maxi=0
	s=0
	while maxi< len(sequences): 
		maxi=min([25*s+25,len(sequences)])
		cprobe_index=range(25*s+0,maxi)
		#print cprobe_index
		#print "\n"

		csequences=[sequences[i] for i in cprobe_index]

		s=s+1
		sequences_in_one="".join(csequences)
		
		blat_values = {'hgsid' : '441376325_ViPefWQZdY1UkwTASgKSALLXwC4A',
			  #'org' : 'Mouse',
			  'db' : genome,
			  'userSeq': sequences_in_one,
			  'type': '1'
			  }
		first_blat_page=get_page(blat_url,blat_values)
		#write_page(first_blat_page,"first_blat_output"+datetime.datetime.now().strftime("%Y%m%d%H%M%S")+".html")
		
		pn1=re.compile("browser.*HREF.(.*?)>details")
		view_link=pn1.findall(first_blat_page)
		view_link
		
		pn_id=re.compile("details.*>(.*)")
		alignment_matches=pn_id.findall(first_blat_page)

		alignment_info=get_alignment_info(alignment_matches)
		
		#print "alignment info<br>"
		#print(alignment_info)
		## if there is more than one result, save this into the "problematic" list
		number_of_matches=len(view_link)
		#print "number of matches"
		#print number_of_matches
		
		#print "html results"
		
		#print view_link
		## make new id in case there are duplicates... fill in list using id as location 

		#initialize
		for i in range(len(alignment_info)):
			ind_alignment=alignment_info[i]

			id=ind_alignment[0]
			probe_ids[id]=[]
			probes_in_dict[id]=[]

			
		for i in range(len(alignment_info)):
			ind_alignment=alignment_info[i]
			#print ind_alignment
			alignment_url=view_link[i]
			id=ind_alignment[0]
			new_val=len(probe_ids[id])+1
			probe_ids[id].append(new_val)
			probes_in_dict[id].append([ind_alignment,alignment_url])

			
		# ## make new id in case there are duplicates... and store the probe info
		# for k,v in probe_ids.iteritems():
			# number_of_matches_current_probe=len(v)
			# #print "number of matches: "+str(number_of_matches_current_probe)
			# for iv in v:
				# #print str(k)+"_"+str(iv)
				# ##print alignment_info
		
	return probe_ids,probes_in_dict


	
def get_alignments_from_blat(probe_id_dict,probe_blat_dict):
	probe_blat_alignment={}
	
	## make new id in case there are duplicates... and store the probe info
	for k,v in probe_id_dict.iteritems():
		
		probe_blat_alignment[k]=[]
		
		number_of_matches_current_probe=len(v)
		#print "number of matches: "+str(number_of_matches_current_probe)
		
		if number_of_matches_current_probe>1:
			#print("found more than one match")
			multiple_match.append(k)
		
		current_probe_info=probe_blat_dict[k]
		
		#print(len(current_probe_info))
		#print(current_probe_info)
		
		
		# for iv in v:
			# print str(k)+"_"+str(iv)
			# ##print alignment_info
				
		for i in range(len(current_probe_info)):
		
			# print "current probe info============================="
			# print current_probe_info
			blat_link=current_probe_info[i][-1]
			
			### dist to the left
			distl=current_probe_info[i][0][1]
			### dist to the right
			distr=current_probe_info[i][0][2]
			
			## need to get the right distance, in the case that the sequence was shorter than 20 and had additional nucleotides appended
			
			# [
			# 
			#	 [
			#		 ['p2f_match_1_dl_16_dr_15', '16', '15', '50', '1', '50', '50', '100.0%', '13', '+', '110395341', '110395390', '50'],
			#		 '"../cgi-bin/hgc?o=110395340&g=htcUserAli&i=../trash/hgSs/hgSs_genome_3398_c80070.pslx+..%2Ftrash%2FhgSs%2FhgSs_genome_3398_c80070.fa+p2f_match_1_dl_16_dr_15-dist1_16-dist2_15-ext_15&c=chr13&l=110395340&r=110395390&db=mm10&hgsid=507833609_bN7dhs3ohmDva4ECFbjRLIedO1HS"'
			#	 ]
			# 
			#]

			#print blat_link
			next_url=base_url+blat_link.strip('"')[2:]
			results_page=get_page(next_url,[])
			
			pn2=re.compile("SRC.(.*)N.*body")
			mainframe_link=pn2.findall(results_page)
			alignment_url=base_url+mainframe_link[0].strip(' "')[2:]
			alignment_page=get_page(alignment_url,[])
			alignment_page_replaced=alignment_page.replace("\n"," ")
			pn3=re.compile("cDNA.p.*?<PRE>(.*?)</PRE>")
			alignment_diagram=pn3.findall(alignment_page_replaced)[0]
			# print "alignment diagram"
			# print alignment_diagram
			
			pn4=re.compile("Alignment of p.* and (.*)</H2")
			chr_location=pn4.findall(alignment_page_replaced)[0]
			#print chr_location
			
			#print "current probe info"
			#print current_probe_info[i]
			
			#strand=current_probe_info[i][0][7]
			
			#after adding distance to the left and to the right, strand becomes the 9th
			strand=current_probe_info[i][0][9]
			
			
			#print "strand"
			#print strand
			
			#seq_alignment_key=str(new_probe_id)+"."+str(i)
			
			
			### could modify the alignment diagram here to highlight the core sequence only
			### also count mismatches in core sequence
			
			## get dna sequence only
			
			
			#pn_seq_only=re.compile(">([actgACTG]+)<") #missing white space
				
			alignment_diagram_no_space=alignment_diagram.replace(" ","")
			alignment_diagram_seq_only=''.join(i for i in alignment_diagram_no_space if not i.isdigit())
			
			# print "alignment_diagram_seq_only"
			# print alignment_diagram_seq_only
			
			pn_seq_only=re.compile(">\s*([actgACTG]+\s*[actgACTG]*)\s*<")
			
			seq_only="".join(pn_seq_only.findall(alignment_diagram_seq_only))
			#seq_only=seq_only.replace(" ","")
			
			pn_mismatches=re.compile("([actg]+)")
			### trim seq only to remove the extra characters that are not part of the core sequence
			
			#test_str[3:lenstr-5]
			core_sequence=seq_only[int(distl):len(seq_only)-(int(distr))] 
			
			#print "this is the whole alignment diagram"
			#print alignment_diagram
			#print "<br>"
			#print "this is the core sequence"
			#print core_sequence
			#print "<br>"
			#print "seq only"
			#print seq_only
			#print "<br>"
			
			#make a function that counts the letters, in actgACTG and once it reaches the distance where the core sequence starts, add <u> 
			
			mismatches=pn_mismatches.findall(core_sequence)
			
			match_count=len("".join(mismatches))
			
			### change alignment diagram to highlight the core sequence
			#it won't keep the blue color, but the core sequence will be underlined
			
			new_alignment_diagram=seq_only[0:int(distl)]+"<u>"+core_sequence+"</u>"+seq_only[len(seq_only)-int(distr):len(seq_only)]
			
			#print "new alignment"
			#print new_alignment_diagram
			#print "<br>"
			#print "<br>"
			
			probe_blat_alignment[k].append([chr_location,strand,alignment_diagram, new_alignment_diagram,distl,distr,match_count])
			
			#sequence_alignments.append([seq_alignment_key,chr_location,strand,alignment_diagram])
			#sequence_alignments.append((seq_alignment_key,[chr_location,alignment_diagram]))
			#sequence_alignments[seq_alignment_key]=[chr_location,alignment_diagram]
	
	return probe_blat_alignment

#### TODO: finish this later	
#split into the different parts inside the span, and count the characters
# def format_alignment(ucsc_diagram,dist):
	# pn_seq_only=re.compile(">\s*([actgACTG]+)\s*<")
	# seq_matches=pn_seq_only.findall(alignment_diagram)
	
	# ci=0
	
	# new_alignment_diagram=""
	# for cm in seq_matches:
		# new_seq=""
		# for cmcc in cm:
			# if cmcc in "actgACTG":
				# ci=ci+1
				# if ci==dist:
					# ### insert u tag here
				# else:
					# new_seq=new_seq+cmcc
	# return new_alignment_diagram
	
def get_extended_seqs(probe_dict,distances,genome):
	#print probe_dict
	#new_fasta_for_blat_lines=[]
	new_blat_lines=[]
	for k in probe_dict.keys():
	#for k,v in probe_dict.items():
		#loop over multiple matches
		#print k
		v=probe_dict[k]
		mm=v
			
		#print mm
		for im in range(len(mm)):
			#print im
			mmm=mm[im]
			#print "\n\n\ncurrent match: "
			#print mmm
			chr_location=mmm[0]
			[chr,start,end]=get_chrom_location(chr_location)
			
			strand=mmm[1]
			distl=mmm[4]
			distr=mmm[5]
			# if strand=="-":
				# print "neg strand"
			# else:
				# print "pos strand"
			## need to parse location string to get chr, start and end, and store it in data structure
			starts=[start-dist for dist in distances ]
			ends=[end+dist for dist in distances ]
			
			for i in range(len(distances)):
				new_loc_string=chr+":"+str(starts[i])+","+str(ends[i])
				
				#print new_loc_string
				#print "starts "+str(i)+" : "+str(starts[i])
				#print "ends "+str(i)+" : "+str(ends[i])
				
				DNA_seq=get_dna_from_ucsc(new_loc_string,genome)
				
				#p1-dist1_0-dist2_0
				#p25_match1_dist3
				
				
				com_probe_id="".join([k,"_match_"+str(im+1),"_dl_",str(int(distl)+distances[i]),"_dr_",str(int(distr)+distances[i]),"-dist1_",str(int(distl)+distances[i]),"-dist2_",str(int(distr)+distances[i]),"-ext_"+str(distances[i])])
				#com_probe_id="_".join([k,"match"+str(im+1),"dist"+str(distances[i])])
				
				#print "new DNA seq page for "+com_probe_id+"======================="
				#print DNA_seq
				
				### now, need to get the reverse complement if it's neg
				if strand=="-":
					#print "reverse complementing...."
					DNA_seq=ReverseComplement3(DNA_seq)
					#print DNA_seq
				
				# new_fasta_for_blat_lines.append(">"+com_probe_id)
				# new_fasta_for_blat_lines.append(DNA_seq)
				
				new_blat_lines.append([com_probe_id,DNA_seq])
				#print new_blat_lines
				# move this part somewhere else
				# new_fasta_for_blat="\n".join(new_fasta_for_blat_lines)
			
	return 	new_blat_lines	

def get_dict_of_pkeys(dict1,dict2):
	oid_to_pid={}
	for okey in dict1.keys():
		for key in dict2.keys():
			tok=key.split("_")
			cid=tok[0]
			pn2=re.compile(okey+"$")
			m=pn2.findall(cid)
			if(len(m)==1):
				# print tok
				# print "found "+okey
				
				if okey not in oid_to_pid.keys():
					oid_to_pid[okey]=[]
				
				oid_to_pid[okey].append(key)
	return oid_to_pid
	
def get_dict_of_pkeys_round2(dict1,dict2):
	oid_to_pid={}
	for okey in dict1.keys():
		for key in dict2.keys():
			tok=key.split("_")
			cid=tok[0]
			pn2=re.compile(okey+"$")
			m=pn2.findall(cid)
			if(len(m)==1):
				# print tok
				# print "found "+okey
				
				if okey not in oid_to_pid.keys():
					oid_to_pid[okey]=[]
				
				oid_to_pid[okey].append(key)
	return oid_to_pid	

# ####### main
# multiple_match=[]

# base_url="http://genome.ucsc.edu/"
# blat_url="http://genome.ucsc.edu/cgi-bin/hgBlat"

# probe_set=read_in_probes("probe_list_input.tab")
# probes_in_one_fasta=probe_list_to_fasta(probe_set)

# probe_ids,probes_in_dict=do_blat_multiple(probes_in_one_fasta)
# probe_alignment_dict=get_alignments_from_blat(probe_ids,probes_in_dict)
# ## now extend using 3,5,7 and put all the strings in one list

# distances=[3,5,10]

# new_fasta_for_blat=get_extended_seqs(probe_alignment_dict,distances)
# probe_ids_second_round,probes_in_dict_second_round=do_blat_multiple(new_fasta_for_blat)
# probe_alignment_dict_second_round=get_alignments_from_blat(probe_ids_second_round,probes_in_dict_second_round)

# output_html=[]
# for k,v in probe_alignment_dict_second_round.iteritems():
	# output_html.append(v[0][2])
	
	
# "<br><br><br>".join(output_html)
			
#################### works up to here!!!	


#################### back here	
def identify_multiple_matches(probe_dict):
	ids_w_multiple=[]
	for k,v in probe_dict.iteritems():
		if len(v)>1:
			#print "multiple match for probe "+ str(k)
			ids_w_multiple.append(k)
			#print "more info for probe"
			#print probes_in_dict_second_round[k]
			#print probe_dict[k]
	return ids_w_multiple
	
	
# #################################################################### left over
	
	# for k,v in probe_ids.iteritems():
		# number_of_matches=len(v)
		
		# if number_of_matches>1:
			# print "more than one match for probe" + k

	# #### from old ==============================================================
			
		# if number_of_matches>1 :
		# print("found more than one match")
		# multiple_match.append(new_probe_id)
		# #then need to highlight this probe as problematic
	
		# for i in range(number_of_matches):
			# next_url=base_url+view_link[i].strip('"')[2:]
			# results_page=get_page(next_url,[])
			
			# pn2=re.compile("SRC.(.*)N.*body")
			# mainframe_link=pn2.findall(results_page)
			# alignment_url=base_url+mainframe_link[0].strip(' "')[2:]
			# alignment_page=get_page(alignment_url,[])
			# alignment_page_replaced=alignment_page.replace("\n"," ")
			# pn3=re.compile("cDNA.Your.*?<PRE>(.*?)</PRE>")
			# alignment_diagram=pn3.findall(alignment_page_replaced)[0]
			# pn4=re.compile("Alignment of YourSeq and (.*)</H2")
			# chr_location=pn4.findall(alignment_page_replaced)[0]
			
			
			# alignment_details_list=alignment_details[i].split()
			# strand=alignment_details_list[7]
			
			# seq_alignment_key=str(new_probe_id)+"."+str(i)
			# sequence_alignments.append([seq_alignment_key,chr_location,strand,alignment_diagram])
			# #sequence_alignments.append((seq_alignment_key,[chr_location,alignment_diagram]))
			# #sequence_alignments[seq_alignment_key]=[chr_location,alignment_diagram]
	
# ############### end of leftover
	
# def do_blat(new_probe_id,sequence):
	# sequence_alignments=[]
	
	# blat_values = {'hgsid' : '441376325_ViPefWQZdY1UkwTASgKSALLXwC4A',
          # 'org' : 'Mouse',
          # 'db' : 'mm10',
		  # 'userSeq': sequence,
		  # 'type': '1'
		  # }
	# first_blat_page=get_page(blat_url,blat_values)
	# write_page(first_blat_page,"first_blat_output.html")
	
	# pn1=re.compile("browser.*HREF.(.*?)>details")
	# view_link=pn1.findall(first_blat_page)
	# view_link
	
	# ## also need to get the strands
	# pn_details=re.compile("details.*>(.*)")
	# alignment_details=pn_details.findall(first_blat_page)
	
	# #just to test
	# # print alignment_details[0].split()
	# # print "strand"
	# # print alignment_details[0].split()[7]
	
	
	# ######################## 8/31 maybe leave this for second round of blat
	
	
	# ## if there is more than one result, save this into the "problematic" list
	# number_of_matches=len(view_link)
	
	# print new_probe_id+" matches "+str(number_of_matches)+"\n"
	
	# if number_of_matches>1 :
		# print("found more than one match")
		# multiple_match.append(new_probe_id)
		# #then need to highlight this probe as problematic
	
	# for i in range(number_of_matches):
		# next_url=base_url+view_link[i].strip('"')[2:]
		# results_page=get_page(next_url,[])
		
		# pn2=re.compile("SRC.(.*)N.*body")
		# mainframe_link=pn2.findall(results_page)
		# alignment_url=base_url+mainframe_link[0].strip(' "')[2:]
		# alignment_page=get_page(alignment_url,[])
		# alignment_page_replaced=alignment_page.replace("\n"," ")
		# pn3=re.compile("cDNA.[0-9]+.*?<PRE>(.*?)</PRE>")
		# alignment_diagram=pn3.findall(alignment_page_replaced)[0]
		# pn4=re.compile("Alignment of YourSeq and (.*)</H2")
		# chr_location=pn4.findall(alignment_page_replaced)[0]
		
		
		# alignment_details_list=alignment_details[i].split()
		# strand=alignment_details_list[7]
		
		# seq_alignment_key=str(new_probe_id)+"."+str(i)
		# sequence_alignments.append([seq_alignment_key,chr_location,strand,alignment_diagram])
		# #sequence_alignments.append((seq_alignment_key,[chr_location,alignment_diagram]))
		# #sequence_alignments[seq_alignment_key]=[chr_location,alignment_diagram]
	
	# time.sleep( 15 )
# #	return sequence_alignments
	# return (new_probe_id,sequence_alignments)


# # this will be the list where the problematic probes will go	
# multiple_match=[]


# probe_set=read_in_probes("probe_list_input.tab")

# base_url="http://genome.ucsc.edu/"
# blat_url="http://genome.ucsc.edu/cgi-bin/hgBlat"
# dna_url="http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment="

# seq_alignments_first_round=[]
# for individual_probe in probe_set:
	
	# probe_id=individual_probe[0]
	
	# seq_alignments_first_round.append(do_blat(probe_id,individual_probe[1]))
	

# print seq_alignments_first_round

# seq_alignments_fr_dict=dict(seq_alignments_first_round)

# distances=[3,5,10]
# for k,v in seq_alignments_fr_dict.iteritems():
	# #loop over multiple matches
	# for mm in v:
		# print "current match: "
		# print mm
		# chr_location=mm[1]
		# [chr,start,end]=get_chrom_location(chr_location)
		
		# strand=mm[2]
		# if strand=="-":
			# print "neg strand"
		# else:
			# print "pos strand"
		# ## need to parse location string to get chr, start and end, and store it in data structure
		# starts=[start-dist for dist in distances ]
		# ends=[end+dist for dist in distances ]
		
		# for i in range(len(distances)):
			# new_loc_string=chr+":"+str(starts[i])+","+str(ends[i])
			
			# print new_loc_string
			# #print "starts "+str(i)+" : "+str(starts[i])
			# #print "ends "+str(i)+" : "+str(ends[i])
			
			# new_loc_dna_url=dna_url+new_loc_string
			# print new_loc_dna_url
			# new_dna_page=get_page(new_loc_dna_url,[])
			# print new_dna_page
			
			
			
