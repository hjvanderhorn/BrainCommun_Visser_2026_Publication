#!/usr/bin/env python3

from os.path import expanduser
import requests
import os, sys,ssl
from collections import OrderedDict
import argparse
import fileinput
import subprocess

def url_consulta(key,value,web):
	data = {key : value}
	response = (requests.post(web, data,verify=False)).text
	os.system("clear")
	return response

def process_species(show,species):
	count=0
	for element in show:
		if element: #avoid empty elements
			element=element.split(",")
			kingdom=element[0]
			scientific=element[5]
			db=element[3]
			index_bowtie=element[4]
			if db:
				species[db]=[scientific,"OFF",index_bowtie,kingdom]
				count=count+1
				species = OrderedDict(sorted(species.items(), key=lambda t: t[1]))
	return (species,count)


php="https://bioinfo2.ugr.es/sRNAtoolboxDB/populate.php"
php2="https://bioinfo2.ugr.es/sRNAtoolboxDB/update.php"

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument ('--dbdir','-db',type=str ,help ='Database dir for sRNAtoolbox',default="/opt/sRNAtoolboxDB/")

args = parser.parse_args()


if not os.path.exists(args.dbdir):
	os.system("mkdir "+args.dbdir)

	#Add to path
homeDir = expanduser("~")
filename = homeDir+"/.bashrc"
text_to_search = "toolboxDB=/opt/sRNAtoolboxDB"
text_to_replace = "toolboxDB="+args.dbdir

with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
	for line in file:
		if "toolboxDB=" in line:
			lineNew = "toolboxDB="+args.dbdir+"\n"
			print(line.replace(line,lineNew),end='')
		else:
			print(line,end='')
		
species = 'GRCh38_p13'

response=url_consulta("getspecies",species,php)
response=response.split("\n")

log = []
enlace={}
selection = {}

for element in response:
	element = element.replace("/opt/sRNAtoolboxDB", "/opt/sRNAtoolboxDB")
	if element:
			element=element.split("=")
			sp=element[0]
			element=element[1].split(";")
			file=element[0]
			file=file.split("/")
			file=file[-1]				
			log_element=sp+"="+element[1]+file
			enlace[element[0]]=element[1]
			selection[species]=enlace
			if not log_element in log:
				log.append(log_element)



for sample in selection:
    enlaces = selection[sample]
    count = 0
    total = len(enlaces)
    for elemento in enlaces:
        if elemento:
            output = enlaces[elemento].replace("/opt/sRNAtoolboxDB/", args.dbdir)
            input_url = elemento
            if "mysql" in input_url:
                continue
            if not os.path.isdir(output):
                os.makedirs(output)  # Prefer os.makedirs to create intermediate directories if needed
            count += 1
            print(f"Downloading {elemento}, file {count} of {total}")
            cmd = f"wget --no-check-certificate -c -S -r -np -nH --cut-dirs=3 -R index.html* {input_url} -P {output}"
            try:
                subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"Downloaded {sample}, file {count} of {total} successfully")
            except subprocess.CalledProcessError as e:
                print(f"Error downloading {sample}, file {count} of {total}: {e}")
			


