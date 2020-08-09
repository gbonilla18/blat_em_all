#!/usr/bin/python

##import modules for CGI handling 
import cgi, cgitb 
import probe_blat

# Create instance of FieldStorage 
form = cgi.FieldStorage() 
selected_probes= form.getlist('pid')

file_output=form.getlist('file_output')


# HTTP Header
#print "Content-Type:application/octet-stream; name=\"probe_list.txt\"\r\n";
print "Content-Disposition: attachment; filename=\"probe_list.txt\"\r\n\n";



#print selected_probes
#for p in selected_probes:
for p in file_output:
	print p	


# Actual File Content will go here.
# fo = open("foo.txt", "rb")

# str = fo.read();
# print str

# # Close opend file
# fo.close()