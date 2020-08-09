#!/usr/bin/python

##import modules for CGI handling 
import cgi, cgitb 
import probe_blat
import urllib

# Create instance of FieldStorage 
form = cgi.FieldStorage() 

# Get data from fields
p_textarea = form.getvalue('probe_text')
p_file = form.getvalue('probe_file')
p_ext=form.getvalue('probe_ext')
p_genome=form.getvalue('genome_text').strip()



print "Content-type: text/html"
print ""
print ""

### main processing
#print p_textarea
#print "more......<br>"
#print "genome <br>"
#print p_genome

probe_set=probe_blat.read_in_probes_textarea(p_textarea)

#print probe_set
original_probe_dict=probe_blat.probe_list_to_dict(probe_set)
#print original_probe_dict
probe_blat.multiple_match=[]

user_dist=p_ext.split(",")


### set distance to extend sequence
try:
	distances=[]
	for n in user_dist:
	   val = int(n)
	   distances.append(val)
	   #print n
	   #print "found int"
except ValueError:
	distances=[10]
	#print("That's not an int!")
   
   
		
hgsid="441376325_ViPefWQZdY1UkwTASgKSALLXwC4A"		
probe_blat.base_url="http://genome.ucsc.edu/"
probe_blat.blat_url="http://genome.ucsc.edu/cgi-bin/hgBlat"


#probe_set=probe_blat.read_in_probes_file("probe_list_input.tab")
le_distances=[0 for i in range(0,len(probe_set))]
re_distances=[0 for i in range(0,len(probe_set))]

probes_in_fasta_list=probe_blat.probe_list_to_fasta_list_w_dist(probe_set,"p",le_distances,re_distances)

probe_ids,probes_in_dict=probe_blat.do_blat_multiple(probes_in_fasta_list,p_genome)

#print "<br>"
#print "probe_ids after first round=========="
#print probe_ids
#print "<br>"
probe_alignment_dict=probe_blat.get_alignments_from_blat(probe_ids,probes_in_dict)
## now extend using 3,5,7 and put all the strings in one list


new_sequence_set=probe_blat.get_extended_seqs(probe_alignment_dict,distances,p_genome)
new_fasta_for_blat=probe_blat.probe_list_to_fasta_list(new_sequence_set,"")
probe_ids_second_round,probes_in_dict_second_round=probe_blat.do_blat_multiple(new_fasta_for_blat,p_genome)
#print "<br>"
#print "probes ids after second round=========="
#print probe_ids_second_round

#print "<br>"
#print "probes in dict after second round=========="
#print probes_in_dict_second_round
##print "<br>"
probe_alignment_dict_second_round=probe_blat.get_alignments_from_blat(probe_ids_second_round,probes_in_dict_second_round)

### identify probes with multiple matches
ids_p_w_multi_sec_round=probe_blat.identify_multiple_matches(probe_alignment_dict_second_round)

ids_p_w_multi_fst_round=probe_blat.identify_multiple_matches(probe_alignment_dict)


### get ids from first and second round based on the first round

first_round_dict=probe_blat.get_dict_of_pkeys(original_probe_dict,probe_ids)

second_round_dict=probe_blat.get_dict_of_pkeys(original_probe_dict,probe_ids_second_round)
	
# print "<br>----------------------"	
# print "second round dict"
# print second_round_dict
# print "<br>----------------------"	


### instead, need a function that would add checkbox and format each item with probe id, chr location and original probe (get from dict), 
### also, want to highlight the probes that have duplicates

# output_html=[]
# for k,v in probe_alignment_dict_second_round.iteritems():
	# output_html.append(v[0][2])

## need a function to generate HTML from dict items

### add ucscchr_loc="%s:%s-%s"%(chr,start,end)
def get_ucsc_link(chr_loc,hgsid,genome):
	base_link="https://genome.ucsc.edu/"
	safe_chr_loc=urllib.quote(chr_loc)
	ucsc_query_str="cgi-bin/hgTracks?db=%s&position=%s"%(genome,safe_chr_loc)
	return base_link+ucsc_query_str



output_html=[]
for k,v in original_probe_dict.iteritems():
	header_str="<div class='entry'><input type='checkbox' name='pid' value='"+k+"'> <h4><span id='pid_"+k+"' class='pid'>Probe "+k+"</span></h4><h4> "+v+"</h4>"
	### get first round matches
	content_str=" <strong> first round match </strong>"
	if k in first_round_dict.keys():
		fr_ids=first_round_dict[k]
		for cid in fr_ids:
			#print cid
			mm=probe_alignment_dict[cid]
			if cid in ids_p_w_multi_fst_round:
				multi_class="multiple"
			else:
				multi_class="single"
				
			for cmm in mm:
				
				content_str=content_str+"<br><span class='"+multi_class+"'>"+cid+"</span><br><a href='"+get_ucsc_link(cmm[0],"",p_genome)+"' target='_blank'>"+" ".join(cmm[0:2])+"<br>"+cmm[2]+"<br>"+cmm[3]+"</a><br>mismatches: "+str(cmm[6])
				
				if cmm[6]>0:
					content_str=content_str+"<br/><span class='multiple'>This alignment contains mismatches. The sequences on the second round will contain the corrected genomic sequence. Note that this is an artifact </span>"
				
		### get second round matches
		content_str=content_str+""
		output_html.append(header_str+content_str)

		content_str="<strong> second round match </strong>"
		if k in second_round_dict.keys():
			sc_ids=second_round_dict[k]
			for cid in sc_ids:
				#print cid
				mm=probe_alignment_dict_second_round[cid]
				if cid in ids_p_w_multi_sec_round:
					multi_class="multiple"
				else:
					multi_class="single"
					
				for cmm in mm:
					
					content_str=content_str+"<br><span class='"+multi_class+"'>"+cid+"</span><br>"+" ".join(cmm[0:2])+"<br>alignment from UCSC: <b>"+cmm[2]+"</b><br>core sequence:<b>"+cmm[3]+"</b><br>distances:<b>"+",".join(cmm[4:6])+"</b><br>mismatches: <b>"+str(cmm[6])+"</b>"
					
			### get second round matches
			content_str=content_str+"</div>"
			output_html.append(content_str)
			
html_res="<br>".join(output_html)	
		
## hidden values for file download
output_download=[]
for k,v in original_probe_dict.iteritems():
	header_str=k+"\t"+v+"\t"
	
	content_str=""
	### get first round matches
	if k in first_round_dict.keys():
		fr_ids=first_round_dict[k]
		for cid in fr_ids:
			#print cid
			mm=probe_alignment_dict[cid]				
			for cmm in mm:
				
				output_download.append(header_str+"\t"+"\t".join(cmm[0:2])+"\t"+cmm[2]+"\t"+cmm[3]+"\t-\t"+str(cmm[6]))
				
		### get second round matches

		if k in second_round_dict.keys():
			sc_ids=second_round_dict[k]
			for cid in sc_ids:
				#print cid
				mm=probe_alignment_dict_second_round[cid]
					
				for cmm in mm:
					output_download.append(header_str+"\t"+"\t".join(cmm[0:2])+"\t"+cmm[2]+"\t"+cmm[3]+"\t"+",".join(cmm[4:6])+"\t"+str(cmm[6]) )
					


## all probes
# print "\n".join(output_download)
file_res="\n".join(output_download)

## all probes with multiple matches first round

## all probes with multiple matches second round
	
html_res="<br><br><br>".join(output_html)

html_w_template='''<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="keywords" content="" />
<meta name="description" content="" />
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>BLAT to the Future!</title>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
<link href='http://fonts.googleapis.com/css?family=Abel' rel='stylesheet' type='text/css'>
<link href='https://fonts.googleapis.com/css?family=Alegreya+Sans+SC:800' rel='stylesheet' type='text/css'>
<link href="/probe_scanner/style.css" rel="stylesheet" type="text/css" media="screen" />
<script src="/probe_scanner/process_probes.js"></script>
</head>
<body>
<div id="wrapper">
	<div id="header-wrapper" class="container">
	<div id="header" class="container">
		<div id="logo">
			<h1><a href="#">Blat 'Em All </a></h1>
		</div>
		<div id="menu">
			
		</div>
	</div>
	</div>
	<!-- end #header -->
	<div id="page">
		<div id="content">
			<div class="post">
				<h2 class="title"><a href="#">Blat 'Em All</a></h2>
				<!--<p class="meta"><span class="date">August 28, 2012</span><span class="posted">Posted by <a href="#">Someone</a></span></p>-->
				<div style="clear: both;">&nbsp;</div>
				<form enctype="multipart/form-data" id="probe_list" action="download_probes.py" method="post" >
					
				<input type="submit" value="Download Selected Sequences" class="download">
				<br><br>
					%s
				<br><br>		
				<input type="submit" value="Download Selected Sequences" class="download">
				<input type="hidden" name="file_output" value="%s">
				</form>
			</div>
			
			<div style="clear: both;">&nbsp;</div>
		</div>
		<!-- end #content -->
		<div id="sidebar">
			<ul>
				<li>
					<div id="search" >
						<form method="get" action="#">
							<div>
								
							</div>
						</form>
					</div>
					<div style="clear: both;">&nbsp;</div>
				</li>
				
			</ul>
		</div>
		<!-- end #sidebar -->
		<div style="clear: both;">&nbsp;</div>
	</div>
	<div class="container"><img src="/probe_scanner/images/img03.png" width="1000" height="40" alt="" /></div>
	<!-- end #page -->
</div>
<div id="footer-content"></div>
<div id="footer">
	<p>Copyright (c) 2012 Sitename.com. All rights reserved. Design by <a href="http://www.freecsstemplates.org">FCT</a>.</p>
</div>
<!-- end #footer -->
</body>
</html>			
'''%(html_res,file_res)


print html_w_template
