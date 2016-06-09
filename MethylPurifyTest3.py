#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys,re,os,sys,random,math
from optparse import OptionParser
import logging
from pkg_resources import resource_filename

try:
	import numpy
except ImportError:
	print >> sys.stderr, "No numpy module"

try:
	from pyfasta import Fasta
except ImportError:
	print >> sys.stderr, "No pyfasta module"

random.seed(1234)

global logfhd

def writelog(path):
	logfhd = open(os.path.join(path, "log"), "w")

	logging.basicConfig(level=20,
			   format='%(levelname)-5s @ %(asctime)s: %(message)s',
			   datefmt='%a, %d %b %Y %H:%M:%S',
			   stream=sys.stderr,
			   filemode="w"
			   )
	error = logging.critical
	warn = logging.warning

	return logfhd

def info(a, path):
	logfhd = writelog(path)
	logging.info(a)
	logfhd.write(a+"\n")
	logfhd.flush()

def EM(READ):
	m1 = random.random()
	m2 = random.random() ## initial values of m1 and m2
	alpha1 = random.random()

	steps = 0 ## number of EM steps
	while True:
		nu1 = 0
		de1 = 0
		nu2 = 0
		de2 = 0
		
		nu_alpha = 0
		for key in READ:
			if READ[key]["All"] >= 1:
				p1 = (m1**READ[key]["Methy"])*((1-m1)**(READ[key]["All"] - READ[key]["Methy"]))
				p2 = (m2**READ[key]["Methy"])*((1-m2)**(READ[key]["All"] - READ[key]["Methy"]))
				Q1 = alpha1*p1 / (alpha1*p1 + (1-alpha1)*p2)
				Q2 = (1-alpha1)*p2 / (alpha1*p1 + (1-alpha1)*p2)
				nu1 += Q1 * READ[key]["Methy"]
				de1 += Q1 * READ[key]["All"]
				nu2 += Q2 * READ[key]["Methy"]
				de2 += Q2 * READ[key]["All"]

				nu_alpha += Q1
        #print de1
        #print de2
        #sys.exit(0)
		if (de1*de2 ==0): ## invalid, return a random result
			return (m1,m2,alpha1)
			break

		m1_new = nu1 / de1
		m2_new = nu2 / de2

		alpha1_new = nu_alpha/len(READ)

		if (abs(m1_new - m1) < 0.01 and abs(m2_new - m2) < 0.01 and abs(alpha1_new - alpha1) < 0.01) or steps > 200: ## num of steps <= 200
			if (abs(m1_new - m1) < 0.01 and abs(m2_new - m2) < 0.01 and abs(alpha1_new - alpha1) < 0.01):

                    # alpha1 is the smaller one
				if (alpha1 > 0.5): # switch
					alpha1 = 1 - alpha1
					temp = m1
					m1 = m2
					m2 = temp
				return (m1,m2,alpha1)
				break
			elif steps > 200:
                    #print "m1,m2,alpha1 not converge in 200 steps by EM calculating!\n"

                    # alpha1 is the smaller one
				if (alpha1 > 0.5): # switch
					alpha1 = 1 - alpha1
					temp = m1
					m1 = m2
					m2 = temp
				return (m1,m2,alpha1)
				break
		m1 = m1_new
		m2 = m2_new
		alpha1 = alpha1_new
		steps = steps + 1

def bootstrap(READ, repeated_times):
	num_of_read = len(READ)

	M1 = [0]*repeated_times
	M2 = [0]*repeated_times
	alpha1 = [0]*repeated_times

	for i in range(repeated_times):
		SAMPLE_key = []
		sample_times = min(1000, num_of_read)

		for j in range(sample_times):
			SAMPLE_key.append(random.sample(READ.keys(),1)[0])

		SAMPLE = {}
		for key in SAMPLE_key:
			SAMPLE[key] = READ[key]
		(m1_s, m2_s, al_s) = EM(SAMPLE)

		alpha1[i] = al_s
		M1[i] = m1_s
		M2[i] = m2_s
		del SAMPLE_key
		del SAMPLE

	var_alpha = numpy.var(alpha1)
	var_M1 = numpy.var(M1)
	var_M2 = numpy.var(M2)

	return (var_alpha, var_M1, var_M2)

def candidate_selection(READ):
	N_methy = 0
	N_unmethy = 0
	N_read = len(READ)

	upper = 0.8
	lower = 0.2
	effective_ratio = 0.9
	minor_ratio = 0.05

	for key in READ:
		if READ[key]["All"] > 0:
			methyLevel = READ[key]["Methy"]*1.0/READ[key]["All"]
			if methyLevel >= upper:
				N_methy += 1
			elif methyLevel <= lower:
				N_unmethy += 1

	if (N_methy+N_unmethy)*1.0/N_read >= effective_ratio and N_methy*1.0/N_read > minor_ratio and N_unmethy*1.0/N_read > minor_ratio:
		return "T"
	else:
		return "F"

def get_mixingRatio(input_file, output_file, bin_length, coverage_cutoff, repeat_times, genome):
	genome_file = genome
	f = Fasta(genome_file)

	IN = open(input_file, "r")
	out = open(output_file, "w")
	out.write("# chr\tstart\tend\tm1\tm2\talpha1\tvar_alpha1\tvar_M1\tvar_M2\tread_count\tcytosine_count_in_CG\n")	

	chr_current = "chr1"
	bin_start = 0
	READ = {}

	for line in IN:
		arr = line.strip().split("\t")
		if arr[0].startswith("@") or not arr:
			continue
		if len(arr[5]) != 4:
			continue

		readName = arr[0]
		chrN = arr[2]
		len_read = int(arr[5][:-1])
		start = int(arr[3])
		end = start + len_read

		read = arr[9]
		refseq = str.upper(f[chrN][start-1:start-1+len_read].__str__())
		strand = arr[-1]

		if chrN != chr_current or start - bin_start > bin_length: ## it
			CGcount = f[chrN][start:start+bin_length].count("cg")+f[chrN][start:start+bin_length].count("CG")
			if CGcount > 10: ## start a new bin
				read_count_cutoff = coverage_cutoff * bin_length * 1.0/len_read
				if len(READ) > read_count_cutoff:
					if candidate_selection(READ) == "T":                     
						(m1,m2,al) = EM(READ)
						(var_alpha,var_M1,var_M2) = bootstrap(READ,repeat_times)
						CGcount = f[chrN][bin_start:bin_start+bin_length].count("cg")+f[chrN][bin_start:bin_start+bin_length].count("CG")
						out.write ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr_current,bin_start,bin_start+bin_length,m1,m2,al,var_alpha,var_M1,var_M2,len(READ),CGcount))


                ## initialization
				READ = {}
				chr_current = chrN
				bin_start = start

				(n_me,n_all) = get_methyl(read,strand,refseq)

				READ[readName] = {}
				READ[readName]["Methy"] = n_me
				READ[readName]["All"] = n_all
			else:
				continue

		else: ## an old bin just add reads

			(n_me,n_all) = get_methyl(read,strand,refseq)

			READ[readName] = {}
			READ[readName]["Methy"] = n_me
			READ[readName]["All"] = n_all

def get_methyl(read, strand, refseq):
	n_me = 0
	n_unme = 0

    #Waston strand

    #------ Watson strand ------#

	if strand == "ZS:Z:++" or strand == "ZS:Z:+-":
		for n in range(len(read) - 1):
			if refseq[n:n+2] == "CG":
				if read[n:n+2] == "CG":
					n_me += 1
				elif read[n:n+2] == "TG":
					n_unme += 1
				else:
                    #print read[n:n+2]
					pass


    #------ Crick strand ------#

	elif strand == "ZS:Z:-+" or strand == "ZS:Z:--":
		for n in range(len(read) - 1):
			if refseq[n:n+2] == "CG":
				if read[n:n+2] == "CG":
					n_me += 1
				elif read[n:n+2] == "CA":
					n_unme += 1
				else:
                    #print read[n:n+2]
					pass

    #print strand,"\n",refseq,"\n",read,"\n","methyl:",n_me,"\t","unmethyl:",n_unme
    #raw_input("Next...")
	return(n_me,n_me + n_unme)

def get_methyl_detail(input_file, input_file2, outputFile, cpg_length, genome):
	""" read first informative bin methylation detail """
	genome_file = genome
	f = Fasta(genome_file)

	IN1 = open(input_file, "r")
	IN2 = open(input_file2, "r")
	out = open(outputFile, "w")
	out.write("#chr\tnumber\tcpg_number\tmethylation_count\n")

	chr_current = "chr10"
	methyl_read_detail = []
	bins_cpg_number = 0	
	
	for line2 in IN2:
		arr2 = line2.strip().split("\t")
		if arr2[0].startswith("X"):
			continue
		Info_bin_start = int(arr2[1])
		bins_cpg_number = 0
		#print "info %d" % (Info_bin_start)

		for line1 in IN1:
			arr1 = line1.strip().split("\t")

			if arr1[0].startswith("@") or not arr1:
				continue
			if len(arr1[5]) != 4:
				continue
		
			chrN = arr1[2]
			len_read = int(arr1[5][:-1])
			start = int(arr1[3])
			#print "seq %d" % (start)
			end = start + len_read
			read = arr1[9]
			refseq = str.upper(f[chrN][start-1:start-1+len_read].__str__())
			strand = arr1[-1]
			#print chrN
			#print chr_current
			if chrN != chr_current:
				chr_current = chrN
			if chrN == chr_current and start < Info_bin_start:
				continue
			if chrN == chr_current and start > Info_bin_start and start - Info_bin_start < cpg_length:
				#print "success"
				methyl_read_detail, read_cpg_number = get_methylread_detail(read, strand, refseq)
				bins_cpg_number = bins_cpg_number + read_cpg_number
				out.write("%s\t%s\t%s\t%s\n" % (chr_current,Info_bin_start,bins_cpg_number,methyl_read_detail))
			if chrN == chr_current and start > Info_bin_start and start - Info_bin_start > cpg_length:
				break

def get_methylread_detail(read, strand, refseq):
	#change every read in informative bins to 0,1 by methylation status
	methyl_list = []
	read_cpg_number = 0
	#------ Watson strand ------#

	if strand == "ZS:Z:++" or strand == "ZS:Z:+-":
		for n in range(len(read) - 1):
			if refseq[n:n+2] == "CG":
				if read[n:n+2] == "CG":
					methyl_list.append(1)
					read_cpg_number = read_cpg_number + 1
				elif read[n:n+2] == "TG":
					methyl_list.append(0)
					read_cpg_number = read_cpg_number + 1
				else:
					pass
	

	#------ Crick strand ------#
	elif strand == "ZS:Z:-+" or strand == "ZS:Z:--":
		for n in range(len(read) - 1):
			if refseq[n:n+2] == "CG":
				if read[n:n+2] == "CG":
					methyl_list.append(1)
					read_cpg_number = read_cpg_number + 1
				elif read[n:n+2] == "CA":
					methyl_list.append(0)
					read_cpg_number = read_cpg_number + 1
				else:
					pass

	return methyl_list, read_cpg_number

def one_For_Each_CGI(fn, cgi, inputFile, outputFile, species):
    """ retain up to one bins for each CGI  """
    Rfile = open(os.path.join(fn, "one_For_Each_CGI.r"),"w")

    Rfile.write("raw_file = '"+inputFile+"'\n")
    Rfile.write("new_file = '"+outputFile+"'\n")
    Rfile.write("dat = read.delim(raw_file)\n")
    Rfile.write("dat_rm_CGI = c()\n")
    Rfile.write("CGI = read.delim(\"%s\", header = F)\n" % cgi)
    Rfile.write("print(head(CGI))\n")
    Rfile.write("total = 0\n")
    Rfile.write("for (i in 1:nrow(CGI)) {\n")
    Rfile.write("\tchr = CGI[i, 1]\n")
    Rfile.write("\tstart = CGI[i, 2]\n")
    Rfile.write("\tend = CGI[i, 3]\n")
    Rfile.write("\tbins = dat[as.character(dat[, 1]) == chr & dat[, 2] >= start & dat[, 2] < end, ]\n")
    Rfile.write("\ttotal = total + nrow(bins)\n")
    Rfile.write("\tif (nrow(bins) > 0) {\n")
    Rfile.write("\t\tdat_rm_CGI = rbind(dat_rm_CGI, bins[which.min(bins[, 'var_M1']), ])\n")
    Rfile.write("\t}\n}\n\n")

    Rfile.write("print(total)\n")
    Rfile.write("colnames(dat_rm_CGI) = colnames(dat)\n")
    Rfile.write("write.table(dat_rm_CGI, file = new_file, sep = '\t',quote = F, row.names = F)\n")

    Rfile.close()

    os.system("Rscript %s " % os.path.join(fn, "one_For_Each_CGI.r"))

def get_bestBins(fn, inputFile, outputFile):
	
    Rfile = open(os.path.join(fn, "get_best_bins.r"),"w")
	
    Rfile.write("raw_file = '"+inputFile+"'\n")
    Rfile.write("new_file = '"+outputFile+"'\n")
    Rfile.write("dat = read.delim(raw_file)\n")
    Rfile.write("dat2 = dat[dat[, 'alpha1'] >= 0.28,]\n")
    Rfile.write("dat3 = dat2[dat2[, 'alpha1'] <= 0.32,]\n")
    Rfile.write("colnames(dat3) = colnames(dat)\n")
    Rfile.write("write.table(c(dat3), file = new_file, sep = '\t',quote = F, row.names = F)\n")

    Rfile.close()
	
    os.system("Rscript %s " % os.path.join(fn, "get_best_bins.r"))
	
def read_mixingRatio(fn,inputFile):
    """ get mixing ratio from bed file """

    Rfile = open(os.path.join(fn, "get_composition.r"),"w")
    Rfile.write("dat = read.delim('"+inputFile+"')\n")
    Rfile.write("dat2 = dat[dat[, 'var_M1'] <= 0.01,]\n")
    Rfile.write("order = order(dat2[,'var_M1'],decreasing=F)\n")
    Rfile.write("x = dat2[order[1:500],'alpha1']\n")

    Rfile.write("factorx <- factor(cut(x, breaks=seq(0,0.5,0.01)))\n")
    Rfile.write("xout <- as.data.frame(table(factorx))\n")
    Rfile.write("peak = xout[which.max(xout[,'Freq']),'factorx']\n")
    Rfile.write("peak = as.character(peak)\n")
    Rfile.write("start = strsplit(substr(peak,2,nchar(peak)-1),',')[[1]][1]\n")
    Rfile.write("end = strsplit(substr(peak,2,nchar(peak)-1),',')[[1]][2]\n")
    Rfile.write("summit = mean(c(as.numeric(start),as.numeric(end)))\n")
    Rfile.write("write.table(c(summit,nrow(dat2)),file = '%s', row.names = F,col.names = F)\n\n" % os.path.join(fn, "alpha1.pred"))

    Rfile.close()

    os.system("Rscript %s" % os.path.join(fn,"get_composition.r"))


def main():	
	usage = "usage: python %prog <-f filename> <-g ref_genome> [...]"
	description = "Select reads from two sam files"

	op = OptionParser(version="%prog 0.1", description=description, usage=usage, add_help_option=False)

	op.add_option("-f", "--filename", dest="filename", type="str",
			     help="The file name of mixing tissue, only accept bam file currently")
	op.add_option("-b","--BinLength",dest="BinLength",type="int",default="300",
			     help="Length of each bin, default is 300")
	op.add_option("-c","--coverage_cutoff",dest="coverage_cutoff",type="int",default="20",
			     help="Lowest coverage cutoff in each bin, default is 20")
	op.add_option("-s","--SamplingTimes",dest="SamplingTimes",type="int",default="50",
                 help="sampling times for bootstraping in each bin, default is 50")
	op.add_option("--species",dest="species",type="str",default = "hg19",
			     help="the ref genome used for mapping, default is 'hg19'")
	op.add_option("-g","--genome",dest="genome",type="str",
			     help="the ref genome fasta used for building your genome index")
	op.add_option("-i","--cpgi",dest="cgi",type="str",
			     help="CpG island bed file")
	(options,args) = op.parse_args()

	if not options.filename:
		op.print_help()
		sys.exit(1)

	filename = options.filename
	bin_length = options.BinLength
	coverage_cutoff = options.coverage_cutoff
	repeat_times = options.SamplingTimes
	species = options.species
	genome = options.genome
	cgi = options.cgi

	if not os.path.exists(cgi):
		print >>sys.stderr, "Not found CpG island files"
		sys.exit(1)

	if not os.path.exists(genome):
		print >> sys.stderr, "Not found fasta files"
		sys.exit(1)

	if not os.path.exists(filename):
		print >> sys.stderr, "Not found bam files"
		sys.exit(1)

	fn = os.path.split(filename)[1].strip(".bam")

	#------- step 1: map to CGI and sort --------#
	if not os.path.isdir(fn):
		os.mkdir(fn)
	path = fn

	if not os.path.exists(os.path.join(fn, fn+".CGI.bam")):
		CMD_map_to_CGI = "samtools view -h -b -L "+ cgi + " " + filename + " > "+ os.path.join(fn, fn+".CGI.bam")
		info("Runing: %s" % CMD_map_to_CGI, path)
        os.system(CMD_map_to_CGI)

	if not os.path.exists(os.path.join(fn, fn+".CGI.sorted.bam")):
		CMD_sort = "samtools sort "+ os.path.join(fn, fn + ".CGI.bam") + " " +" " + os.path.join(fn, fn+".CGI.sorted")
		info("Runing: %s" % CMD_sort, path)
		os.system(CMD_sort)

	if not os.path.exists(os.path.join(fn, fn+".CGI.sorted.sam")):
		CMD_bam_to_sam = "samtools view -h "+ os.path.join(fn, fn+".CGI.sorted.bam") + " > " + os.path.join(fn, fn+".CGI.sorted.sam")
		info("Running: %s" % CMD_bam_to_sam,path)
		os.system(CMD_bam_to_sam)	

	#-------- step 2: get mixingRatio --------#
	Input = os.path.join(fn, fn+".CGI.sorted.sam")
	Output_mixRatio = os.path.join(fn, fn+"Informative_bins.bed")

	if not os.path.exists(os.path.join(fn, fn+"Informative_bins.bed")):
		info("Running get_mixingRatio:\n Input:%s\n Output:%s\n Bin_length:%s\n Coverage:%s\n Repeat_times:%s\n Genome:%s\n" % (Input,Output_mixRatio,bin_length,coverage_cutoff,repeat_times,species), path)
		get_mixingRatio(Input,Output_mixRatio,bin_length,coverage_cutoff,repeat_times,genome)
	
	get_bestBins(fn, Output_mixRatio, Output_mixRatio+".Bestbins")		
	Input1 = os.path.join(fn, fn+".CGI.sorted.sam")
	Input2 = os.path.join(fn, fn+"Informative_bins.bed.Bestbins")
	Output_detail = os.path.join(fn, fn+"methyl_detail.bed")
	get_methyl_detail(Input1, Input2, Output_detail, bin_length, genome)
	one_For_Each_CGI(fn, cgi, Output_mixRatio, Output_mixRatio+".OneForCGI", species)

	read_mixingRatio(fn, Output_mixRatio+".OneForCGI")
	MixFile = open(os.path.join(fn, "alpha1.pred"),"r")
	MixRatio = float(MixFile.readline().strip())
	NumInfoBins = float(MixFile.readline().strip())

	print "The predicted mixing ratio is: ",MixRatio
	print "The number of informative bins is: ",NumInfoBins	

	if NumInfoBins < 400:
		print "The number of informative bins is less than 400, too less for the prediction of mixing ratio, continue"
#        sys.exit(0)
	elif NumInfoBins >= 400 and NumInfoBins < 500:
		print "Warning: the number of informative bins is less than 500, not sufficient to get a reliable estimation!"

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt\n")
		sys.exit(0)





