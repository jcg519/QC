#/usr/bin/python
#!-*- coding:utf-8 -*-

'''
23 mofang OR values QC (Hom,LD and so on remove)

'''

import os,sys
import re
import MySQLdb 
import numpy as np
from optparse import OptionParser
 
def parseCommand():
	usage = "usage: ./QC.py -1 input1 -o outputfile"
	version = "%prog 1.0"
	parser = OptionParser(usage = usage, version = version)
	parser.add_option("-1", "--input1", dest = "input1", help = "the rsid input file,rs1000597 CEU rs10033464 CEU rs1004819 CEU rs10050860 CEU rs10086908 CEU")
	parser.add_option("-o", "--output", dest = "output", help = "the output file")
	return parser.parse_args()

def getQC(input1,output):
	LD= "plink" + ".ld"
	CHB="CHB"
	file1=open(input1,'r')
	file2=open(output,'a')
	os.system("plink --bfile " + CHB + " --r2 --ld-snp-list " + input1 + " --ld-window-r2 0.8  --noweb")
	if not os.path.exists(LD):
		print 'PLINK run erro'
	else:
#		ldd="###### pairwise LD values #######"+"\n"
#		file2.write(ldd)
		os.system("awk  '{if($3!=$6)print$3,$6,$7}' " + LD + " > " + output)
	rsids=[]
	for line in file1.readlines():
		rsids.append(re.split('\n',line)[0])
	file1.close()
	for rs in rsids:
		res=re.split(' ',rs)
		rsid=res[0].strip()
		pop=res[1].strip()
		rr=[]
		if pop !='CHB':	
			db = MySQLdb.connect(host='192.168.30.252', port=3306, db='gendb', user='dna', passwd='dna', charset='utf8') 
			# 使用cursor()方法获取操作游标  	
			cursor = db.cursor()
			SQL = "select Genotype1,Genotype2,Genotype3 from my_population_genotype where rsid='%s' and GenotypeOwn='CHB'" % (rsid)
			cursor.execute(SQL)
			result = cursor.fetchall()
			genotype={}
			for value in result:
				for i in range(0,3):
					key = value[i].split(':')
					basekey = key[0].replace(',', '')
					baseValue = key[1].split('|')[1]
					genotype[basekey]=baseValue
			chb_rate=[]
			for i in list(genotype.values()):
				chb_rate.append(float(i))
			ma=max(chb_rate)
			rr.append(chb_rate)		
			if float(ma)>0.95:
				chunhe="###### 基因频率纯合######"+"\n"
				file2.write(chunhe)
				rsidlist=str(rsid)+ " " + str(genotype)+"\n"
				file2.write(rsidlist)

			gene="select Genotype1,Genotype2,Genotype3 from my_population_genotype where rsid='%s' and GenotypeOwn='%s'" % (rsid,pop)
			cursor.execute(gene)
			genevalues = cursor.fetchall()
			genotypes={}
			for value in genevalues:
				for i in range(0,3):
					key = value[i].split(':')
					basekey = key[0].replace(',', '') 
					baseValue = key[1].split('|')[1]
					genotypes[basekey]=baseValue
			other_rate=[]
			for i in list(genotypes.values()):
				other_rate.append(float(i))
			rr.append(other_rate)
			cor=np.corrcoef(rr)			
			cr=cor[0][1]
			if float(cr)<-0.9:
				words="###### 人种基因频率与中国人的基因频率相反 ######"+"\n"
				pearson=str(rsid)+" " +str(genotype)+ " "+str(genotypes)+"\n"
				file2.write(words)
				file2.write(pearson)

		if pop =='CHB':
			chb = MySQLdb.connect(host='192.168.30.252', port=3306, db='gendb', user='dna', passwd='dna', charset='utf8')
			cursor.execute(chb)
			info = cursor.fetchall()
			geno={}
			for value in info:
				for i in range(0,3):
					key = value[i].split(':')
					basekey = key[0].replace(',', '')
					baseValue = key[1].split('|')[1]
					geno[basekey]=baseValue
			han_rate=[]
			for i in list(geno.values()):
				han_rate.append(float(i))
			mam=max(han_rate)
			if float(mam)>0.95:
				ll= "###### 基因频率纯合######"+"\n"
				llist= str(rsid)+ " " +str(genotype)+"\n"
				file2.write(ll)
				file2.write(llist)
#	os.system("plink --bfile " + CHB + " --r2 --ld-snp-list " + input1 + " --ld-window-r2 0.8  --noweb")
#	os.system("awk  '{if($3!=$6)print$3,$6,$7}' " + LD + " >> " + output)	
	
	file2.close()

if __name__=="__main__":
	(options, args) = parseCommand()
	getQC(options.input1,options.output)	 
						
