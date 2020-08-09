#!/usr/bin/python

import cgi, cgitb

print "Content-type: text/html"

print ""
print '''<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="keywords" content="" />
<meta name="description" content="" />
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>BLAT to the Future!</title>
<link href='http://fonts.googleapis.com/css?family=Abel' rel='stylesheet' type='text/css'>
<link href='https://fonts.googleapis.com/css?family=Alegreya+Sans+SC:800' rel='stylesheet' type='text/css'>
<link href="/probe_scanner/style.css" rel="stylesheet" type="text/css" media="screen" />
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
	<div><img src="/probe_scanner/images/img03.png" width="1000" height="40" alt="" /></div>
	</div>
	<!-- end #header -->
	<div id="page">
		<div id="content">
			<div class="post">
				<h2 class="title"><a href="#">Blat 'Em All</a></h2>
				<div style="clear: both;">&nbsp;</div>
				<div class="entry">
					<p>Provide a list of probes, with ID and DNA-sequence separated by tab </p>
					<form enctype="multipart/form-data" action="process_probe_input.py" method="post" >
						<p>
						<textarea name="probe_text" cols="40" rows="4">
						Paste the Stellaris output here
						</textarea>
						</p>
						<!--
						<p>... or upload a tab-separated file</p>
						<p>
						<input type="file"  name="probe_file"  />
						</p>-->

						<p>Enter the bp by which you want to extend the sequences, separated by comma</p>
						<p> for example, enter 3,5 to extend by 3 and 5 bp

						</p>
						<textarea name="probe_ext" cols="40" rows="4">
						Enter distances here
						</textarea>
						</p>

						<p>Enter the genome you want to BLAT against</p>
						<p> For example, mm9, mm10, rn6, etc.

						</p>
						<textarea name="genome_text" cols="40" rows="4">
						mm9
						</textarea>
						</p>

						<input type="submit" value="BLAT" />

					</form>

				</div>
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
'''
