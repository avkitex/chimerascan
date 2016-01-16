#!/usr/bin/python

#########################################################################################
#	Autor Avktiex, kit.iz.179@gmail.com
#	Copyright BOSTONGENE 2016
#	
#	A programm for converting GTF and UCSC_tables HG19 files to chimerascan-0.5.0-compatable format
#	
#	Download GTF from Ensembl website ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
#	Add header to it #seqname	source	feature	start	end	score	strand	frame	attribute
#	
#	Download custom tracks from UCSC tables browser
#	(all from knownGene except cdsStart, cdsEnd, proteinID, alignID; geneSymbol from hg19.kgXref; clusterId	from hg19.knownIsoforms; value from hg19.knownToEnsembl)
#	The final file will have the header:
#	#hg19.knownGene.name	hg19.knownGene.chrom	hg19.knownGene.strand	hg19.knownGene.txStart	hg19.knownGene.txEnd	hg19.knownGene.exonCount	hg19.knownGene.exonStarts	hg19.knownGene.exonEnds	hg19.kgXref.geneSymbol	hg19.knownIsoforms.clusterId	hg19.knownToEnsembl.value
#
#	RUN: program -gtf GTFFile -ucsc hg19_ucsc_annotation.txt -o output.txt
##########################################################################################

from __future__ import print_function
import argparse, sys


def getGtfLine(gtfFile):
	line = gtfFile.readline()
	hList = []
	hList.append('seqname')
	hList.append('source')
	hList.append('feature')
	hList.append('start')
	hList.append('end')
	hList.append('score')
	hList.append('strand')
	hList.append('frame')
	hList.append('attribute')
	if line[0] != '#':
		print("There is no annotation in gtf. Try 'seqname	source	feature	start	end	score	strand	frame	attribute'")
		sys.exit(1)
	header = line[1:].strip().split('\t')
	for hItem in hList:
		if hItem not in header:
			print("There is no " + hItem + " in header in gtf file")
			print("There must be all " + ' '.join(hList))
			sys.exit(1)
	line = gtfFile.readline()
	while line:
		entry = dict(zip(header, line.strip().split('\t')))
		attributes = entry['attribute'][:-1].strip().split(';')
		for attr in attributes:
			k = attr.strip().split(' ')
			#print(k)
			k[1] = k[1].strip('"')
			entry[k[0]] = k[1]
		del entry['attribute']
		yield entry
		line = gtfFile.readline()
	gtfFile.close()
	
def getUcscLine(ucscFile):
	line = ucscFile.readline()
	hList = []
	hList.append('hg19.knownGene.name')
	hList.append('hg19.knownGene.chrom')
	hList.append('hg19.knownGene.strand')
	hList.append('hg19.knownGene.txStart')
	hList.append('hg19.knownGene.txEnd')
	hList.append('hg19.knownGene.exonCount')
	hList.append('hg19.knownGene.exonStarts')
	hList.append('hg19.knownGene.exonEnds')
	hList.append('hg19.kgXref.geneSymbol')
	hList.append('hg19.knownIsoforms.clusterId')
	hList.append('hg19.knownToEnsembl.value')

	if line[0] != '#':
		print("There is no annotation in gtf. Try '#hg19.knownGene.name	hg19.knownGene.chrom	hg19.knownGene.strand	hg19.knownGene.txStart	hg19.knownGene.txEnd	hg19.knownGene.exonCount	hg19.knownGene.exonStarts	hg19.knownGene.exonEnds	hg19.kgXref.geneSymbol	hg19.knownIsoforms.clusterId	hg19.knownToEnsembl.value'")
		sys.exit(1)
	header = line[1:].strip().split('\t')
	for hItem in hList:
		if hItem not in header:
			print("There is no " + hItem + " in header in ucsc file")
			print("There must be all " + ' '.join(hList))
			sys.exit(1)
	line = ucscFile.readline()
	while line:
		if line[0] == '#':
			continue
		entry = dict(zip(header, line.strip().split('\t')))
		yield entry
		line = ucscFile.readline()
	ucscFile.close()
	
def outLineResult(conTranscript, delim='\t'):
	fields = [conTranscript['hg19.knownGene.chrom'],
			conTranscript['hg19.knownGene.txStart'],
			conTranscript['hg19.knownGene.txEnd'],
			conTranscript['transcriptSerial'],
			conTranscript['hg19.knownIsoforms.clusterId'],
			conTranscript['hg19.knownGene.strand'],
			conTranscript['hg19.knownGene.exonCount'],
			conTranscript['hg19.knownGene.exonStarts'],
			conTranscript['hg19.knownGene.exonEnds'],
			conTranscript['gtf']['gene_biotype'],
			','.join([conTranscript['hg19.knownToEnsembl.value'],conTranscript['gtf']['transcript_name'], conTranscript['hg19.knownGene.name']]) + ',',
			','.join([conTranscript['gtf']['gene_name'], conTranscript['gtf']['gene_id'], conTranscript['hg19.kgXref.geneSymbol']]) + ',',
			','.join([conTranscript['gtf']['source'], 'ucsc', conTranscript['gtf']['gene_source'], conTranscript['gtf']['transcript_source']]) + ','
	]
	return delim.join(fields)  
	
def main():
	parser = argparse.ArgumentParser(prog='hg_gtd_to_chimera.py', usage='%(prog)s [options]', description='description', epilog="\xa9 BostonGene 2015")

	parser.add_argument('-gtf', '--gtf', metavar='GTF', type = argparse.FileType('r'), help = 'GTF', required=True)
	parser.add_argument('-ucsc', '--ucsc', metavar='UCSC', type = argparse.FileType('r'), help = 'HG_tables', required=True)
	parser.add_argument('-o', '--output', metavar='OUT', type = argparse.FileType('w'), help = 'outFile', required=True)
	args = parser.parse_args()

	ensemblTr = {}
	transcriptSerial = 1
	
	print("Processing GTF file")
	for entry in getGtfLine(args.gtf):
		if entry['feature'] == "transcript":
			print('.', end = '')
			#TODO cheks for all fields
			ensemblTr[entry['transcript_id']] = entry
	print("\nProcessing ucsc file")
	for transcript in getUcscLine(args.ucsc):
		trEnsId = transcript['hg19.knownToEnsembl.value']
		if trEnsId == 'n/a':
			continue
		if trEnsId not in ensemblTr:
			print("\nThere is no transcript " + trEnsId + ' in GTF file')
			continue
		if transcript['hg19.knownGene.chrom'] != 'chr' + ensemblTr[trEnsId]['seqname']:
			print("\nDifferent chromosomes in UCSC and GTF,", transcript['hg19.knownGene.chrom'], 'chr' + ensemblTr[trEnsId]['seqname'])
			continue
		if transcript['hg19.kgXref.geneSymbol'] != ensemblTr[trEnsId]['gene_name']:
			print("\nDifferent geneName(symbol) in UCSC and GTF,", transcript['hg19.kgXref.geneSymbol'], ensemblTr[trEnsId]['gene_name'])
			#continue
		if transcript['hg19.knownGene.strand'] != ensemblTr[trEnsId]['strand']:
			print("\nDifferent strand in UCSC and GTF,", transcript['hg19.kgXref.geneSymbol'], ensemblTr[trEnsId]['strand'])
			continue
		print('.', end = '')
		transcript['gtf'] = ensemblTr[trEnsId]
		transcript['transcriptSerial'] = str(transcriptSerial)
		transcriptSerial += 1
		print(outLineResult(transcript), file=args.output)
	print("\nDone")
	

	args.output.close()
	
main()
	