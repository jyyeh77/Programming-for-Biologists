#!/usr/bin/python

'''
AthBlast.cgi is a CGI script that outputs a HTML form supporting an user-friendly BLAST interface for searches against Arabidopsis nucleotide and protein databases.  It is capable of running all five BLAST alignments, BlastN, BlastP, BlastX, TBlastN and TBlastX.

This script contains a variety of functions to validate DNA and protein identity for either uploaded files or sequences entered into
a text box.  There is also a function check_FASTA() that verifies that an uploaded file is in FASTA format.  These functions are
used in the HTML form to minimize user error when selecting the appropriate databases and alignment types for their query sequence.
'''

import cgi
import cgitb; cgitb.enable()
import subprocess
import os
from Bio.Blast import NCBIXML
from getExpPlot import *
##########################################################################################################################

# FUNCTIONS

def blast_n(input_file, upload_file = False):
	'''
	Runs BlastN alignment in the command-line.
	Outputs XML file blastn.xml.
	'''
	
	# calls BlastN on input text or FASTA format file
	blast_output = subprocess.call(['blastn', '-query', input_file, '-db', 'Arabidopsis.nucl',
	'-outfmt', '5', '-out', 'blastn.xml'])
	
def blast_p(input_file, upload_file = False):
	'''
	Runs BlastP alignment in the command-line.
	Outputs XML file blastp.xml.
	'''
	
	# calls BlastP on input text or FASTA format file
	blast_output = subprocess.call(['blastp', '-query', input_file, '-db', 'Arabidopsis.prot', 
		'-outfmt', '5', '-out', 'blastp.xml'])

def blast_x(input_file, upload_file = False):
	'''
	Runs BlastX alignment in the command-line.
	Outputs XML file blastx.xml.
	'''
	
	# calls BlastX on input text or uploaded FASTA file
	blast_output = subprocess.call(['blastx', '-query', input_file, '-db', 'Arabidopsis.prot', 
		'-outfmt', '5', '-out', 'blastx.xml'])

def tblast_n(input_file, upload_file = False):
	'''
	Runs TBlastN alignment in the command-line.
	Outputs XML file tblastn.xml.
	'''
	
	# calls TBlastN on input text or uploaded FASTA file
	blast_output = subprocess.call(['tblastn', '-query', input_file, '-db', 'Arabidopsis.nucl',
        	'-outfmt', '5', '-out', 'tblastn.xml'])

def tblast_x(input_file, upload_file = False):
	'''
	Runs  TBlastX alignment in the command-line.
	Outputs XML file tblastx.xml.
	'''
	
	# calls TBlastX on input text or FASTA format file
	blast_output = subprocess.call(['tblastx', '-query', input_file, '-db', 'Arabidopsis.nucl',
        '-outfmt', '5', '-out', 'tblastx.xml'])


	


def blast_parse(input_xml):
	'''
	Opens input XML file for ONE query and reads it using NCBIXML module.
	Stores in a dictionary the name of the query, the subject name,
	the e-value, the percent identity match of the hit, the start of the query,
	the end of the query, and the start and end of the subject.

	Extracts this information from all the high-scoring pairs from the alignment.
	Returns list of dictionaries containing alignment info from each hit.
	'''
	
	# opens input XML file and parses it 	
	blast_record = NCBIXML.read(open(input_xml, "r"))
	
	# stores the query name	
	queryname = blast_record.query
	
	# stores # of bp in query for % identity calculation
	queryletters = blast_record.query_letters

	# creates list to store alignment information for query	
	blast_data = []
	
	# iterates through all hits
	for alignment in blast_record.alignments:
		
		# iterates through all high-scoring pairs for each hit			
		for hsp in alignment.hsps:
			
			# calculates percent identity match to query
			percent_id = (100.0 * hsp.identities)/queryletters
			
			# stores alignment information for each HSP in dictionary form
			rows = {
                        "QueryName": queryname,
                        "SubjectName": alignment.hit_def,
			"E-value": hsp.expect,
			"PercentIdentity": percent_id,
			"Query_start": hsp.query_start,
			"Query_end": hsp.query_end,
			"Subject_start": hsp.sbjct_start,
			"Subject_end": hsp.sbjct_end,
	                }
			blast_data.append(rows)
	
	return blast_data

def check_DNA(input_sequence):
	'''
	Checks if input sequence is comprised of A, T, C, G nucleotides.
	Is case-insensitive and can account for spaces in sequence.
	'''

	# Possible DNA bases
	DNA = 'ATCGatcg '
	
	# Sets default flag to True
	flag = "True"
	# Checks if every base in input_sequence is a DNA nucleotide
	for base in input_sequence:
		# sets flag to False if base not DNA nucleotide
		if base not in DNA:
			flag = "False"
	return flag

def check_DNA_FASTA(fasta_file):
	''' 
	Checks if sequence in FASTA format with header is comprised of
	A, T, C, G nucleotides, is case-insensitive and accounts for spaces.
	'''
	
	# Possible DNA bases
	DNA = 'ATCGatcg ' 
	
	# Creates list to hold sequence in FASTA file
	base_sequence = []

	# opens FASTA file for reading
	file = open(fasta_file)

	# skips header
	header = file.next()

	# iterates through each line of the sequence
	for line in file:
		# removes newline breaks
		newline = line.strip()
		
		# stores lines of sequence into list
		base_sequence.append(newline)

	# stores sequence in FASTA file as continuous string
	str_sequence = "".join(base_sequence)

	# sets default flag to True
	flag = "True"
	
	# iterates through each base in sequence, sets flag to False if not in DNA
	for base in str_sequence:
		if base not in DNA:
			flag = "False"
	file.close()
	return flag
	

def check_protein(input_sequence):
	'''
	Checks if input sequence is comprised of amino acids.
	Returns False if input sequence contains letters outside of
	20 possible amino acids.
	Overlaps with DNA nucleotides, so cannot be used to distinguish
	DNA and protein sequences.
	'''

	# changes all amino acids in sequence to upper case
	sequence = input_sequence.upper()	

	# Possible amino acids
	Amino_acids = "ACDEFGHIKLMNPQRSTVWY "
	flag = "True"

	# Checks if units in input sequence are amino acids
	for amino_acid in sequence:
		
		if amino_acid not in Amino_acids:
			flag = "False"
	
	return flag

def check_protein_FASTA(fasta_file):
	'''
	Checks if sequence in FASTA file format is comprised of amino acids.
	Returns False if input sequence contains characters not in possible amino acids.
	Overlaps with DNA nucleotides, cannot be used to distinguish DNA and protein sequences.
	'''
	
	# Possible amino acids
	Amino_acids = "ACDEFGHIKLMNPQRSTVWY "

	# sets default flag to True
	flag = "True"

	# creates list to store amino acids in sequence to be checked
	aa_sequence = []

	# opens FASTA file for reading
        file = open(fasta_file)

	# skips header
        header = file.next()

	# iterates through lines in sequence	
        for line in file:
		# removes newline breaks at the end of each line
                newline = line.strip()
		# adds each sequence line to list
                aa_sequence.append(newline)
	
	# joins protein sequence into continuous string
        str_sequence = "".join(aa_sequence)
	
	# checks if each character in sequence is an amino acid
	for aa in str_sequence:
		if aa not in Amino_acids:
			flag = "False"
	file.close()
	return flag

def check_FASTA(upload):
	'''
	Checks if uploaded file is in FASTA format, checking for presence
	of a ">" at the start.  Returns False if otherwise.
	'''
	
	# opens file for reading
	file = open(upload, "r")

	# reads each line in file
	lines = file.readline()

	# creates empty string to store flag value
	flag = ""

	# returns False if beginning of file doesn't start with >
	if not lines.startswith(">"):
		flag = "False"
	return flag

#print check_protein_FASTA("query.fa")
#check_DNA_FASTA("query.fa")
#print check_DNA_FASTA("query.fa")
#with open("query.fa") as file:
 #       contents = file.read()

#if check_DNA(contents) is "True":
#        print "no"
#else:
#        print "noo"
	
#try:	
#	import msvcrt
#	msvcrt.setmode(0, os.0_BINARY)
#	msvcrt.setmode(1, os.0_BINARY)
#except ImportError:
#	pass
###########################################################################################################################3
# MAIN

# uses FieldStorage module to extract HTML values
args = cgi.FieldStorage()

# sets variable for text inputted in text box
input_txt = args.getvalue('input_seq')

# stores uploaded FASTA file
uploadedfile = args.getvalue('file')

# stores option for BlastN alignment
blastn = args.getvalue('BlastN')

# stores option for BlastP alignment
blastp = args.getvalue('BlastP')

# stores option for BlastX alignment
blastx = args.getvalue('BlastX')

# stores option for TblastN alignment
tblastn = args.getvalue('TblastN')

# stores option for TblastX alignment
tblastx = args.getvalue('TblastX')

# stores checkbox option for nucleotide database
nucdata = args.getvalue('nuc_db')

# stores checkbox option for protein database
protdata = args.getvalue('prot_db')
db_flag = ""

# writes text in text box to local file
if input_txt:
	filename = "query.fa"
	output = open(filename, 'wb')
	output.write(input_txt)
	output.close()

# writes uploaded FASTA file to local file
if uploadedfile:
	filename = "query.fa"
	FUPOUT = open(filename, 'wb')
	FUPOUT.write(uploadedfile)
	FUPOUT.close()

# runs BlastN on local sequence file and parses alignment information in list of dictionaries
if blastn:
	blast_n("query.fa")
	table_data = blast_parse("blastn.xml")
# runs BlastP on local file, parses alignment information and stores in list of dictionaries
elif blastp:
	blast_p("query.fa")		
	table_data = blast_parse("blastp.xml")
# runs BlastX on local file, parses alignment information and stores in list of dictionaries
elif blastx:
	blast_x("query.fa")
	table_data = blast_parse("blastx.xml")
# runs TBlastN on local file, parses alignment and stores information 
elif tblastn:
	tblast_n("query.fa")
	table_data = blast_parse("tblastn.xml")
# runs TBlastX on local file, parses alignment and stores information
elif tblastx:
	tblast_x("query.fa")
	table_data = blast_parse("tblastx.xml")
	
###############################################################################################################################
# HTML FORM

print "Content-type: text/html \r\n\r\n"
print "<html>"
print "<head>"
print "<title>BLAST Interface</title>"
print "</head>"
print "<body>"
print '''<h2 style="background-color:palevioletred">Welcome to BLAST Interface </h2>'''

# creates text box for input sequence
print '''
<form action="AthBlast.cgi" method="post">
Input sequence: <input type=\"text\" name=\"input_seq\"><br>
<p>Enter ONLY sequences with NO spaces.</p>
<p>Make sure to press submit to enter your sequence into the database!</p>
<p><input type="submit" value="Submit"></p>
</form>
'''

# allows for file upload
print '''
<form enctype="multipart/form-data" action="AthBlast.cgi" method="post">
<p>File: <input type=\"file\" name="file"></p>
<p><input type="submit" value="Upload"></p>
</form>
'''

# checks if sequence is inputted through text box
if input_txt:
	
	# opens local file containing sequence in text box
	file = open("query.fa")
        contents = file.read()
	file.close()
        
	# checks that input sequence is either DNA or protein
	if check_DNA(contents) is not "False" or check_protein(contents) is not "False":
		print "<p>You have entered the following sequence as your query.</p>"
		# prints sequence submitted through text box
		print input_txt
		print "<p>Please select <b>ONE</b> of the following databases to query your sequence against after pressing submit.</p>"
		# creates checkboxes for selection of either the Arabidopsis nucleotide or protein database
		print '''
		<form action="AthBlast.cgi" method="post">
		Database type:
		<input type="checkbox" name="nuc_db" value="on"/>Arabidopsis nucleotide database
		<input type="checkbox" name="prot_db" value="on"/>Arabidopsis protein database
		<input type="submit" value="Select Database"/>
		</form>
		'''
			
	else:
		# error message if sequence in text box is neither DNA or protein
		print "<p>Please enter a valid DNA sequence or amino acid sequence.</p>"

# checks that input sequence is submitted as an upload
elif uploadedfile:
	# checks that uploaded file is in FASTA format
	if check_FASTA("query.fa") is not "False":
		
		print "<p>Your file has been uploaded, here is your query sequence.</p>"
		# prints out input sequence in uploaded file
		print "<div>{0}</div>".format(open(r"query.fa", "rb").read())
		print "<p>Please select one of the following databases to query your sequence against!.</p>"
		# creates checkboxes for selection of either the Arabidopsis nucleotide or protein database
		print '''
                	<form action="AthBlast.cgi" method="post">
                	Database type:
                	<input type="checkbox" name="nuc_db" value="on"/>Arabidopsis nucleotide database
                	<input type="checkbox" name="prot_db" value="on"/>Arabidopsis protein database
                	<input type="submit" value="Select Database"/>
                	</form>
                	'''

	else:
		# if check_FASTA returns False, prints this error message
		print "<p>Please enter a valid DNA or amino acid sequence in FASTA format!!.</p>"


# Checks if Arabidopsis nucleotide database is selected
# Allows user to select BlastN/TblastX options for nucleotide query vs. nucleotide database
# Selection of TblastN if protein query vs. nucleotide database
if nucdata:
	file = open("query.fa")
	contents = file.read()
	
	# checks that uploaded file contains DNA sequence and is a FASTA file
	if check_FASTA("query.fa") is not "False" and check_DNA_FASTA("query.fa") is "True":
		print "<p>You are querying a nucleotide sequence against a nucleotide database.</p>"
		# outputs only the BLastN option if sequence in upload is nucleotide
		print "<p>Please click the BlastN or TBlastX alignment to proceed.</p>"
		print '''
		<form action="AthBlast.cgi" method="post">
		<p><input type="submit" value="BlastN" name="BlastN"></p>
		<p><input type="submit" value="TblastX" name="TblastX"></p>
		</form>
		'''
	# checks that sequence from text box is DNA
	elif check_DNA(contents) is "True":
		print "<p>You are querying a nucleotide sequence against a nucleotide database.</p>"
	       	# outputs only the BlastN or TblastX options if sequence from text box is nucleotide
		print "<p>Please click the BlastN or TBlastX alignment to proceed.</p>"
        	print '''
                <form action="AthBlast.cgi" method="post">
                <p><input type="submit" value="BlastN" name="BlastN"></p>
                <p><input type="submit" value="TblastX" name="TblastX"></p>
                </form>
                '''

	# checks that input text contains protein sequence
	elif check_protein(contents) is "True" and check_DNA(contents) is "False":
		print "<p>You are querying a protein sequence against a nucleotide database.</p>"
		# outputs only the TBlastN option if sequence from text box is protein
		print "<p>Please click the TblastN alignment to proceed.</p>"
                print '''
                <form action="AthBlast.cgi" method="post">
		<p><input type="submit" value="TblastN" name ="TblastN"></p>
		</form>
		'''
	
	# checks that uploaded file contains protein sequence and is in FASTA format
	elif check_FASTA("query.fa") is not "False" and check_DNA_FASTA("query.fa") is "False" and check_protein_FASTA("query.fa") is "True":

		print "<p>You are querying a protein sequence against a nucleotide database.</p>"
               	# outputs only the TBlastN option if sequence in uploaded file is protein
		print "<p>Please click the TblastN alignment to proceed.</p>"
                print '''
                <form action="AthBlast.cgi" method="post">
                <p><input type="submit" value="TblastN" name ="TblastN"></p>
                </form>
                '''

# checks if Arabidopsis protein database is selected
elif protdata:
	# opens local file for reading
	with open("query.fa") as file:
		contents = file.read()
	
	# confirms sequence in uploaded file is in FASTA format and is protein
	if check_FASTA("query.fa") is not "False" and check_protein_FASTA("query.fa") is "True" and check_DNA_FASTA("query.fa") is "False":
		print "<p>You are querying a protein sequence against a protein database.</p>"
		print '''
		<form action="AthBlast.cgi" method="post">
		<p><input type="submit" value="BlastP" name="BlastP"></p>
		</form>
		'''
	
	# confirms sequence from text box is protein
	elif check_protein(contents) is "True" and check_DNA(contents) is "False":
		print "<p>You are querying a protein sequence against a protein database.</p>"
                print '''
                <form action="AthBlast.cgi" method="post">
                <p><input type="submit" value="BlastP" name="BlastP"></p>
                </form>
                '''
	
	# confirms sequence in uploaded file is in FASTA format and is DNA
	elif check_FASTA("query.fa") is not "False" and check_DNA_FASTA("query.fa") is "True":
		print "<p>You are querying a nucleotide sequence against a protein database.</p>"
                print "<p>Please click the BlastX alignment to proceed.</p>"
		print '''
		 <form action="AthBlast.cgi" method="post">
                <p><input type="submit" value="BlastX" name="BlastX"></p>
                </form>
		'''
	
	# confirms that sequence from text box is DNA
	elif check_DNA(contents) is "True":
		print "<p>You are querying a nucleotide sequence against a protein database.</p>"
                print "<p>Please click the BlastX alignment to proceed.</p>"
                print '''
                 <form action="AthBlast.cgi" method="post">
                <p><input type="submit" value="BlastX" name="BlastX"></p>
                </form>
                '''
# if user selects both the nucleotide and protein databases, outputs this error message
elif nucdata and protdata:
	print "<p>Please select <b>ONLY ONE database</b>!!!!</p>"

# Messages to confirm which type of Blast has been performed on input sequence
if blastn:
	print "BlastN has been performed on your input sequence!"
	print "blastn.xml is now in your local directory."
elif blastp:
	print "BlastP has been performed on your input sequence!"
	print "blastp.xml is now in your local directory."
elif blastx:
       print "BlastX has been performed on your input sequence!"
       print "blastx.xml is now in your local directory."
elif tblastn:
	print "TBlastN has been performed on your input sequence!"
	print "tblastn.xml is now in your local directory."
elif tblastx:
	print "TBlastX has been performed on your input sequence!"
	print "tblastx.xml is now in your local directory."


# prints alignment information in tabular format for selected alignment option 
if blastn or blastp or blastx or tblastn or tblastx:	
	# checks that hits exist for input query sequence
	if len(table_data) > 0:
		print '''
        	<h2>Alignment Information</h2>
		<p>Hyperlinked genes have expression values under different experimental conditions displayed in bar plot format.</p>
        	<table border="1">'''
		# prints column indices for table
	        print '''<tr><th>Query Name</th><th>Subject Name</th><th>Percent Identity</th><th>E-value</th>
                <th>Query Start</th><th>Query End</th><th>Subject Start</th><th>Subject End</th></tr>'''
		# iterates through the dictionaries corresponding to each hit
		for hit in table_data:
			# prints name of search query
			print "<td>%s</td>" % (hit['QueryName'])
			# calls query_genes() from getExpPlot.py to create expression value bar plot for hits
			query_genes("ArabidopsisExpValues.db", hit['SubjectName'])
			
			# creates hyperlink to expression values graph for hit if hit exists in expression values database
			if os.path.isfile(hit['SubjectName']) is True:
				# creates hyperlink for subject or hit
				print '''<td><a href="%s">%s</a></td>''' % (hit['SubjectName'], hit['SubjectName']) 
			else:
				# does not create hyperlink if hit does not exist in expression values database
				print "<td>%s</td>" % (hit['SubjectName'])
			# prints rest of alignment information in table
			print "<td>%s</td>" % (hit['PercentIdentity'])
			print "<td>%s</td>" % (hit['E-value'])
			print "<td>%s</td>" % (hit['Query_start'])
			print "<td>%s</td>" % (hit['Query_end'])
			print "<td>%s</td>" % (hit['Subject_start'])
			print "<td>%s</td></tr>" % (hit['Subject_end'])
		print "</table>"
 	else:
		# error message if there are no hits from an alignment
		print "<p>No hits sorry :(</p>"



print "</body>"
print "</html>"
