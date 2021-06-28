import os
import numpy as np
import Levenshtein
import xml.etree.ElementTree as ET
import statistics
from scipy import stats

class Analyze_mutation_epitope:
	def remove_breakLines(self, arq):

		f=open("new_"+arq, "w")
		f.close()
		
		f=open(arq, "r")
		c=0
		for line in f:
			l=line.replace("\n","")
			
			if(l.find(">")!=-1 and c>0):
				with open("new_"+arq, "a") as gf:
					gf.write(seq+"\n")
				seq=l.split(".")[0]+"\n"
			else:
				if(c==0):
					seq=l.split(".")[0]+"\n"
				else:
					seq+=l.replace("J","I")
			c+=1	
		f.close()
		
		with open("new_"+arq, "a") as gf:
			gf.write(seq+"\n")

	def translate(self, seq): 
		table = { 
			'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
			'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
			'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
			'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				  
			'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
			'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
			'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
			'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
			'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
			'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
			'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
			'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
			'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
			'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
			'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
			'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
		} 
		protein ="" 
		for i in range(0, len(seq), 3): 
			codon = seq[i:i + 3] 
			if(codon in table.keys()):
				protein+= table[codon] 
		return protein 

	def get_changes(self, epi, refp):
		f=open("temp.fa","w")
		f.write(">ref\n"+refp+"\n")
		f.write(">query\n"+epi+"\n")
		f.close()

		os.system("mafft temp.fa > temp.aln")

		self.remove_breakLines("temp.aln")
		seql=""
		f=open("new_temp.aln","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")!=-1):
				id=l
			else:
				if(id==">query"):
					seql=l
		f.close()

		pos=[]
		i=0
		for s in seql:
			if(s!="-"):
				pos.append(i)
			i+=1
		subref=refp[pos[0]:pos[-1]+1]
		subseq=seql[pos[0]:pos[-1]+1]
		i=0
		ps=pos[0]
		changes=[]
		for s1 in subref:
			if(s1!=subseq[i]):
				changes.append(ps+"_"+s1+"_"+subseq[i])
			i+=1 
			ps+=1

		return ",".join(changes)

	def run_mutation_analysis(self, epitopes):
		self.remove_breakLines("ref.fasta")
		ref=""
		f=open("new_ref.fasta","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")==-1):
				ref=l
		f.close()
		
		c=0
		cord={}
		f=open("coordinate.csv","r")
		for line in f:
			l=line.replace("\n","").split(",")
			if(c>0):
				cord[l[4]]=[int(l[1]), int(l[2])]
			c+=1
		f.close()
		
		"""self.remove_breakLines("brseqs.fasta")
		br={}
		f=open("new_brseqs.fasta","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")!=-1):
				id=l.split("/")[2]
			else:
				br[id]=l
		f.close()"""
		
		f=open("result.tsv","w")
		f.close()
		c=0
		f=open(epitopes,"r")
		for line in f:
			l=line.replace("\n","").split(",")
			if(c>0):
				p=l[0].split("_")[1]
				idseq=l[0].split("_")[0]
				
				if(p.find("nucleoca")!=-1):
					p='nucleocapsid'
				sref=ref[cord[p][0]:cord[p][1]]
				refp=self.translate(sref)
				#sbr=br[idseq][cord[p][0]:cord[p][1]+1]
				res=0
				changes=""
				if(not l[1] in refp):
					res=1
					changes=self.get_changes(l[1], refp)
				with open("result.tsv","a") as gf:
					gf.write("%s\t%s\t%s\t%i\t%s\n" %(l[0], l[1], refp, res, changes) )
			c+=1
		f.close()

	def comparison_epitopes(self, file_epi):
		# 1->2, 1->3, 2->3
		if(not os.path.isdir("compare_epitopes")):
			os.system("mkdir compare_epitopes")
			
		epis={}
		
		""" #input in separate files
		for i in range(1,4):
			f=open("compare_epitopes/class"+str(i)+".tsv","r")
			for line in f:
				l=line.replace("\n","").split("\t")
				epis["c"+str(i)]=l[1]
			f.close()"""
		
		co=0
		f=open(file_epi,"r")
		for line in f:
			if(co>0):
				l=line.replace("\n","")
				if(l!=''):
					l=l.split(",")
					if(l[2].lower()=="class i"):
						key="c1"
					if(l[2].lower()=="class ii"):
						key="c2"
					if(l[2].lower()=="b-cell"):
						key="c3"
					if(not key in epis.keys()):
						epis[key]=[]
					epis[key].append(l[1])
			co+=1
		f.close()
		
		intersecall=np.intersect1d(epis["c1"],epis["c2"],epis["c3"])
		f=open("compare_epitopes/epitopes_in_three.csv","w")
		for e in intersecall:
			f.write("%s\n" %(e) )
		f.close()
		
		compare_list=["1-2","1-3","2-3"]
		for c in compare_list:
			class1=c.split("-")[0]
			class2=c.split("-")[1]
			eps={}
			for e1 in epis["c"+class1]:
				for e2 in epis["c"+class2]:
					eps[e1+"-"+e2]=Levenshtein.ratio(e1,e2)*100
					
			sorted_list = sorted( eps.items(), key=lambda kv: kv[1] )
			f=open("compare_epitopes/epitopes_by_distance_classes"+class1+"-"+class2+".csv","w")
			for e in sorted_list:
				e1=e[0].split("-")[0]
				e2=e[0].split("-")[1]
				f.write("%s\t%s\t%.2f\n" %(e1, e2, e[1]) )
			f.close()
				
	def validate_epitopes(self, file_epi):
		humanprot=[]
		co=0
		seq=""
		f=open("human_38_proteins.faa","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")!=-1):
				if(co>0):
					humanprot.append(seq)

				seq=""
			else:
				seq+=l
			co+=1
		f.close()
		if(seq!=""):
			humanprot.append(seq)

		epis={}
		co=0
		f=open(file_epi,"r")
		for line in f:
			if(co>0):
				l=line.replace("\n","").split(",")
				epis[l[1]]=0
			co+=1
		f.close()

		for gh in humanprot:
			for ep in epis.keys():
				if(gh.find(ep)!=-1):
					epis[ep]+=1

		epis_temp={}
		f=open("epitopes_temp.faa","w")
		i=1
		for ep in epis.keys():
			if(epis[ep]==0):
				epis_temp["epitope"+str(i)]=ep
				f.write(">%s\n%s\n" %("epitope"+str(i), ep) )
			i+=1
		f.close()

		"""os.system("backtranseq -sequence epitopes_temp.faa -outfile validated_epitopes.fasta")
		os.system("Rscript shortmatch_hg38.R validated_epitopes.fasta")
		falses=[]
		f=open("result_shortmatch.csv","r")
		for line in f:
			l=line.replace("\n","").split(";")
			falses.append(l[-1])
		f.close()"""

		f=open("final_validated_epitopes_human_genome.tsv","w")
		for ep in epis_temp.keys():
			#if(not ep in falses):
			f.write("%s\n" %( epis_temp[ep] ) )
		f.close()   

	def analysis_epitope_from_pdb(self, folder):
		aa={}
		f=open("aacodes.csv","r")
		for line in f:
			l=line.replace("\n","").split(",")
			aa[l[1]]=l[0]
		f.close()

		f=open("epitopes_from_pdb.tsv","w")
		f.write("id,sequence,score\n")
		f.close()

		for f in os.listdir(folder):
			if(f.endswith(".txt")):
				name=f.split("_")[0]
				seq=""

				c=0
				arq=open(folder+f, "r")
				for line in arq:
					l=line.replace("\n","").split("\t")
					if(len(l)==6):
						if(c>0 and len(seq)>2):
							with open("epitopes_from_pdb.tsv","a") as gf:
								gf.write("%s,%s,%.3f\n" %(name+"_3D", seq, statistics.mean(scs)) )
						seq=""
						scs=[]
						
					if(len(l)==7):
						seq+=aa[l[2]]
						scs.append(float(l[5]))
						
					c+=1
				arq.close()

				if(c>0 and len(seq)>3):
					with open("epitopes_from_pdb.tsv","a") as gf:
						gf.write("%s,%s,%.3f\n" %(name+"_3D", seq, statistics.mean(scs)) )
						
	def analysis_conservation(self, file_):
		# cat gisaid_* > all_lineages.fasta
		# makeblastdb -in all_lineages.fasta  -dbtype nucl -input_type fasta -out lineages
		
		f=open("result_conservation.tsv","w")
		f.close()
		
		c=0
		f=open(file_, "r")
		for line in f:
			l=line.replace("\n","").split(",")
			if(c>-1):
				fa=open("temp.fasta","w")
				fa.write(">%s\n%s\n" %(l[0], l[1]) )
				fa.close()
				
				high=0
				mod=0
				low=0
				
				seqs=0
				num_aln=0
				"""if(os.path.isfile("blastdb/"+l[2]+".phr")):
					os.system("blastp -db blastdb/"+l[2].lower()+" -query temp.fasta -task blastp-short -out temp.out -outfmt 6")
				
					out=open("temp.out","r")
					for line_ in out:
						l_=line_.replace("\n","").split("\t")
						if(float(l_[2])==100):
							seqs+=1
						if(float(l_[2])>=90):
							high+=1
						if(float(l_[2])>=70 and float(l_[2])<90):
							mod+=1
						if(float(l_[2])<70 ):
							low+=1
							
						num_aln+=1	
					out.close()
				else:"""
				for b in os.listdir("blastdb"):
					if(b.endswith(".phr")):
						b=b.split(".")[0]
						
						os.system("blastp -db blastdb/"+b+" -query temp.fasta -task blastp-short -out temp.out -outfmt 6")
					
						out=open("temp.out","r")
						for line_ in out:
							l_=line_.replace("\n","").split("\t")
							if(float(l_[2])==100):
								seqs+=1
							if(float(l_[2])>=90):
								high+=1
							if(float(l_[2])>=70 and float(l_[2])<90):
								mod+=1
							if(float(l_[2])<70 ):
								low+=1
								
							num_aln+=1	
						out.close()
						
				if(num_aln==0):
					num_aln=13403
					
				#num_aln=15610
				result = seqs/num_aln
				high=(high/num_aln)*100
				mod=(mod/num_aln)*100
				low=(low/num_aln)*100
				
				print(result*100, high, mod, low, num_aln)
				
				with open("result_conservation.tsv","a") as gf:
					gf.write("%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n" %(l[0], l[1], result*100, high, mod, low) )
			c+=1
		f.close()
		
		os.system("rm temp.fasta")
		os.system("rm temp.out")
		
	def check_overlapping(self, file_):
		f=open("result_overlapping.tsv","w")
		f.close()

		pubepis={}
		if(not os.path.isfile("published_epis_sarscov2.tsv")):
			f=open("published_epis_sarscov2.tsv","w")
			f.close()
			
			count=0
			total=len(os.listdir('iedb_export/'))
			for f in os.listdir('iedb_export/'):
				count+=1
				print("Processing "+str(count)+" of "+str(total))

				g=open('iedb_export/'+f)
				temp=g.read()
				g.close()

				fg=temp.replace('xsi:schemaLocation="http://www.iedb.org/schema/CurationSchema http://beta.iedb.org/schema/Curation.xsd" xmlns="http://www.iedb.org/schema/CurationSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"','')
				g=open('iedb_export/'+f,'w')
				g.write(fg)
				g.close()
				
				if(fg.find("<SourceOrganismId>2697049</SourceOrganismId>")!=-1):
					tree = ET.parse('iedb_export/'+f)
					root = tree.getroot()

					refid=""
					article=""
					for ref in root.findall('Reference'):
						for refi in ref.findall('ReferenceId'):
							refid=refi.text

						for art in ref.findall('Article'):
							for pub in art.findall('PubmedId'):
								article=pub.text

					seqs=[]
					ids=[]
					names=[]
					for ref in root.findall('Reference'):
						for epis in ref.findall('Epitopes'):
							for ep in epis.findall('Epitope'):
								for nam in ep.findall('EpitopeName'):
									names.append(nam.text)

								for idd in ep.findall('EpitopeId'):
									ids.append(idd.text)
								
								for estr in ep.findall('EpitopeStructure'):
									for frag in estr.findall('FragmentOfANaturalSequenceMolecule'):
										for seq in frag.findall('LinearSequence'):
											seqs.append(seq.text)

					if(len(seqs)>0):
						with open("published_epis_sarscov2.tsv","a") as gf:
							gf.write("%s\t%s\t%s\t%s\t%s\n" %(refid, article, ",".join(ids), ",".join(names), ",".join(seqs) ) )
		else:
			f=open("published_epis_sarscov2.tsv","r")
			for line in f:
				l=line.replace("\n","").split("\t")
				ids=l[2].split(",")
				c=0
				for s in l[4].split(","):
					if(not s in pubepis):
						pubepis[s]=ids[c]

					c+=1
			f.close()

		report={}
		details={}
		c=0
		f=open(file_, "r")
		for line in f:
			l=line.replace("\n","").split(",")
			if(c>0):
				epi=l[1]

				info=[]

				similar=0
				cutoff=0.5
				for target in pubepis.keys():
					identity=Levenshtein.ratio(epi, target)
					if(identity>cutoff):
						similar+=1
						info.append(pubepis[target]+"-"+str(identity*100) )

				report[epi]=similar
				details[epi]=info

			c+=1
		f.close()

		sorted_list = sorted( report.items(), key=lambda kv: kv[1] )		
		for epi in sorted_list:
			with open("result_overlapping.tsv","a") as gf:
				gf.write("%s\t%s\t%s\n" %(epi[0], epi[1], ",".join(details[epi[0]]) ) )
				
	def calculate_correlation(self, samples, epitopes):
		proteins={}
		f=open(samples, "r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")!=-1):
				id_=l.replace(">","").split(" ")[0]
			else:
				proteins[id_]=len(l)
		f.close()

		epis={}
		c=0
		f=open(epitopes, "r")
		for line in f:
			if(c>0):
				l=line.replace("\n","").split(",")
				if(not l[2] in epis.keys()):
					epis[l[2]]={}
				
				protein=l[0].split("_")[0]+"_"+l[0].split("_")[1]
				if(not protein in epis[l[2]].keys()):
					epis[l[2]][protein]=[]

				epis[l[2]][protein].append(len(l[1]))
			c+=1
		f.close()
		
		f=open("result_correlation.tsv","w")
		f.close()
		for clas in epis.keys():
			size_prot=[]
			size_epi=[]

			for orf in proteins.keys():
				if(orf in epis[clas].keys()):
					size_prot.append(proteins[orf])
					size_epi.append(statistics.mean(epis[clas][orf]))

			slope, intercept, r_value, p_value, std_err = stats.linregress(size_prot, size_epi)
			with open("result_correlation.tsv","a") as gf:
				gf.write("%s\t%.3f\t%.2f\n" %(clas, r_value, p_value) )

	def get_coordinates_epi_in_protein(self, epi, refp):
		f=open("temp.fa","w")
		f.write(">ref\n"+refp+"\n")
		f.write(">query\n"+epi+"\n")
		f.close()

		os.system("mafft temp.fa > temp.aln")

		self.remove_breakLines("temp.aln")
		seql=""
		f=open("new_temp.aln","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")!=-1):
				id=l
			else:
				if(id==">query"):
					seql=l
		f.close()

		pos=[]
		i=0
		
		for s in seql:
			if(s!="-"):
				pos.append(i)
			i+=1
		start=pos[0]
		end=pos[0]+len(epi)-1
		
		muts=[]
		a=0
		i=start
		for s in refp[start:end]:
			if(epi[a]!=s):
				muts.append(s+str(i+1)+epi[a])
			a+=1
			i+=1
		
		return start, end, ','.join(muts)

	def run_get_coordinates_in_protein(self, epitopes):
		self.remove_breakLines("ref.fasta")
		ref=""
		f=open("new_ref.fasta","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")==-1):
				ref=l
		f.close()
		
		c=0
		cord={}
		f=open("coordinate.csv","r")
		for line in f:
			l=line.replace("\n","").split(",")
			if(c>0):
				cord[l[4]]=[int(l[1]), int(l[2])]
			c+=1
		f.close()
		
		"""self.remove_breakLines("brseqs.fasta")
		br={}
		f=open("new_brseqs.fasta","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")!=-1):
				id=l.split("/")[2]
			else:
				br[id]=l
		f.close()"""
		
		f=open("result_epitope_coordinates.tsv","w")
		f.close()
		c=0
		f=open(epitopes,"r")
		for line in f:
			if(line!=""):
				l=line.replace("\n","").split(",")
				if(c>0):
					if(l[0].find("_")!=-1):
						p=l[0].split("_")[1]
						idseq=l[0].split("_")[0]
					else:
						p=l[0].split("_")[0]
						
					sref=ref[cord[p][0]:cord[p][1]]
					refp=self.translate(sref)
					start, end, muts=self.get_coordinates_epi_in_protein(l[1], refp)
					with open("result_epitope_coordinates.tsv","a") as gf:
						gf.write("%s\t%i\t%i\t%s\t%s\t%s\n" %(p, start, end, l[0], l[1], muts) )
			c+=1
		f.close()

	def analysis_coordinates_in_sars2_genome(self, file_):
		id_=file_.split("/")[-1].replace(".tsv","")
		
		# cat gisaid_* > all_lineages.fasta
		# makeblastdb -in all_lineages.fasta  -dbtype nucl -input_type fasta -out lineages
		
		f=open("result_coordinate_genome-"+id_+".tsv","w")
		f.close()
		
		c=0
		f=open(file_, "r")
		for line in f:
			l=line.replace("\n","").split(",")
			if(c>0):
				fa=open("temp.fasta","w")
				fa.write(">%s\n%s\n" %(l[0], l[1]) )
				fa.close()
				
				os.system("tblastn -db lineages -query temp.fasta -num_alignments 2737 -out temp.out -outfmt 6")
				
				counts={}
				
				seqs=0
				out=open("temp.out","r")
				for line_ in out:
					l_=line_.replace("\n","").split("\t")
					if(not l_[8]+"-"+l_[9] in counts.keys()):
						counts[l_[8]+"-"+l_[9]]=0
					counts[l_[8]+"-"+l_[9]]+=1
						
				out.close()
				
				sorted_list = sorted( counts.items(), key=lambda kv: kv[1] )
				start=0
				end=0
				i=0
				for e in sorted_list:
					if(i==0):
						start=int(e[0].split("-")[0])
						end=int(e[0].split("-")[1])
					i+=1
					
				res="not RBD"
				if(start >= 22587 and end <= 23192):
					res="RBD"
				with open("result_coordinate_genome-"+id_+".tsv","a") as gf:
					gf.write("%s\t%s\t%i\t%i\t%s\n" %(l[0], l[1], start, end, res) )
			c+=1
		f.close()
		
		os.system("rm temp.fasta")
		os.system("rm temp.out")
		
class Bepipred:
	def run(self, folder):
		if(not os.path.isdir(folder+"results")):
			os.system("mkdir "+folder+"results")

		for f in os.listdir(folder):
			if(f.endswith(".fasta")):
				name=f.replace(".fasta","")
				os.system("BepiPred-2.0 -t 0.6 "+folder+f+" > "+folder+"results/"+name+".tsv")

	def analyze(self, folder):
		f=open(folder+"epitopes_b.tsv","w")
		f.close()
		for f in os.listdir(folder+"results"):
			id_=f.replace(".tsv","")
			epis=[]
			means=[]
			arq=open(folder+"results/"+f,"r")
			for line in arq:
				l=line.replace("\n", "")
				if(l.find("#")!=-1):
					seq=""
					scs=[]
					ant="."
					if(l.find("###")!=-1):
						idseq=l.split(" ")[2]
				else:
					if(l.find("\t")!=-1):
						info=l.split("\t")
						if(ant=="E" and info[8]=="."):
							epis.append(idseq+"-"+seq)
							means.append(statistics.mean(scs))
							seq=""
							scs=[]
						
						if(info[8]=="E"):
							seq+=info[2]
							scs.append(float(info[7]))
						ant=info[8]
			arq.close()
			
			if(seq!=""):
				epis.append(idseq+"-"+seq)
				means.append(statistics.mean(scs))
			
			c=0
			for epi in epis:
				idseq=epi.split("-")[0]
				s=epi.split("-")[1]
				with open(folder+"epitopes_b.tsv","a") as gf:
					gf.write("%s\t%s\t%s\t%.3f\n" %(id_, idseq, s, means[c]) )
				c+=1

import sys
a=Analyze_mutation_epitope()
if(sys.argv[1]=="1"):
	print("Analyze mutation details such as positions of AA modification")
	a.run_mutation_analysis(sys.argv[2])
if(sys.argv[1]=="2"):
	print("Compare overlapping and similarity between epitopes of distinct classes (I, II and B-cell)")
	a.comparison_epitopes(sys.argv[2])
if(sys.argv[1]=="3"):
	print("Validate epitopes occurrence in human genome and proteins")
	a.validate_epitopes(sys.argv[2])
if(sys.argv[1]=="4"):
	print("Analyze epitopes provided by pdb structure and Discotope")
	a.analysis_epitope_from_pdb(sys.argv[2])
if(sys.argv[1]=="5"):
	print("Analyze conservation of epitopes in lineages genome sequences and portions of genomes that they were in low, moderated or high coverage")
	a.analysis_conservation(sys.argv[2])
if(sys.argv[1]=="6"):
	print("Check epitopes that are already published in IEDB database")
	a.check_overlapping(sys.argv[2])
if(sys.argv[1]=="7"):
	print("Calculate correlation between protein and epitope sizes")
	a.calculate_correlation(sys.argv[2], sys.argv[3])
if(sys.argv[1]=="8"):
	print("Get epitopes coordinates in reference proteins")
	a.run_get_coordinates_in_protein(sys.argv[2])
if(sys.argv[1]=="11"):
	print("Get epitopes coordinates in genome")
	a.analysis_coordinates_in_sars2_genome(sys.argv[2])

b=Bepipred()
if(sys.argv[1]=="9"):
	print("Run bepipred")
	b.run(sys.argv[2]) # Folder with .fasta files
if(sys.argv[1]=="10"):
	print("Analyze output bepipred")
	b.analyze(sys.argv[2]) # Same Folder with .fasta files and the results/ folder should be there 
	

