# Config file for CRUX_release_V1		09-01-2017
# Written by Emily Curd (eecurd@g.ucla.edu) and Gaurav Kandlikar (gkandlikar@ucla.edu)
# Developed at UCLA for the University of California Conservation Consortium's CALeDNA Program



#############################
# Paths to programs / load programs
#############################

#entrez-qiime
ENTREZ_QIIME="${DB}/entrez_qiime.py" 						#path to python script. see documentation for how to obtain this script

#ecoPCR
ecoPCR="${DB}/ecoPCR/src/ecoPCR"							#path to executable. see documentation for how to obtain this program,

#Load / run BLASTn
MODULE_SOURCE="" 											#if none, leave empty <- for HPC
LOAD_BLAST="" 												#if none, leave empty <- for HPC
BLASTn_CMD="${DB}/ncbi-blast-2.6.0+/bin/blastn" 			#either the path to the blastn executable or just blastn if it is loaded or already in your path


#Load / run Qiime
MODULE_SOURCE="source /u/local/Modules/default/init/bash" 	#if none, leave empty <- for HPC
LOAD_QIIME="module load qiime" 								#or what ever code is used to load qiime in a bash shell (e.g. on a mac it might be "macqiime")


###########################
# Paths to reference databases
###########################

#Obitools database folder. Should contain one or more subfolders that contain obitools formatted datbases (e.g. ~/crux_db/Obitools_databases contains EMBL_6162017_std_other, EMBL_6162017_std_pro, EMBL_6162017_std_vert, EMBL_6162017_std_plant.  Each of these subfolders contains the database files with extensions: .sdx, .rdx, and .tdx)
OBI_DB="${DB}/Obitools_databases" 							#see documentation for how to build this library

#NCBI taxonomy dump
TAXO="${DB}/TAXO" 											#see documentation for how to obtain this directory

#NCBI key for accession number to taxonomy assignment
A2T="${DB}/accession2taxonomy/nucl_gb.accession2taxid" 		#see documentation for how to obtain this file




