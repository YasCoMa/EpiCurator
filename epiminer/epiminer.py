import re, string, unicodedata
import inflect
from nltk import word_tokenize, sent_tokenize
from nltk.corpus import stopwords
from nltk.stem import LancasterStemmer, WordNetLemmatizer
import nltk.data
from geniatagger import GENIATagger

from Bio import Entrez
import xml.etree.ElementTree as ET
import os
import itertools
class Literature_hits:

    def get_pmcHits_for_pair(self, epi):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
         
        #handle = Entrez.esearch(db='pmc', sort='relevance', term="epitope ("+epi+")" )
        if(epi.startswith("HLA")):
            handle = Entrez.esearch(db='pmc', sort='relevance', term=""+epi+" and (sars-cov or covid or virus or viral or coronavirus)" )
        else:
            handle = Entrez.esearch(db='pmc', sort='relevance', term=""+epi )
        results = Entrez.read(handle)
        hits = results['IdList']
        return hits

    # Store the title, abstract and sentences for each article (pmcid1.xml, ...) 
    def get_hits_on_pmc_sentences(self, hits, folder):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
        new_hits=0
        
        id=",".join(hits)
        fetch = Entrez.efetch(db='pmc',resetmode='xml',id=id,rettype='full')
        with open(folder+'tempFile.xml', 'w') as f:
            f.write(str(fetch.read()).replace("b\\'","").replace("\\'","").replace("\\n","").replace("b\'","").replace("\'",""))
        tree = ET.parse(folder+'tempFile.xml')
        root = tree.getroot()
        
        if( not os.path.isdir(folder+"pmc_articles") ):
            os.system("mkdir "+folder+"pmc_articles")
        
        # article file creation
        newfile=[]
        pmids_order=[]
        for article in root.findall('article'):
            for front in article.findall('front'):
                for meta in front.findall('article-meta'):
                    for pmids in meta.findall('article-id'):
                        if(pmids.attrib['pub-id-type']=='pmc'):
                            if(not os.path.isfile(folder+"pmc_articles/"+pmids.text+".xml")):
                                newfile.append(True)
                                pmids_order.append(folder+"pmc_articles/"+pmids.text+".xml")
                                with open(folder+"pmc_articles/"+pmids.text+".xml", "a") as g:
                                    g.write("<article>")
                            else:
                                newfile.append(False)
                                pmids_order.append(folder+"pmc_articles/"+pmids.text+".xml")
                        
        #print(pmids_order)
        #print(len(pmids_order), len(newfile))
        # getting title
        c=0
        for article in root.findall('article'):
            for front in article.findall('front'):
                for meta in front.findall('article-meta'):
                    for group in meta.findall('title-group'):
                        for title in group.findall('article-title'):
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("<title>%s</title>" %(title.text) )
            c+=1

        # getting abstract
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("<abstract>")
            
            c+=1
        c=0
        for article in root.findall('article'):
            for front in article.findall('front'):
                for meta in front.findall('article-meta'):
                    for abs in meta.findall('abstract'):
                        for sec in abs.findall('sec'):
                            for blocks2 in sec.findall('p'):
                                complete_text=ET.tostring(blocks2, encoding="unicode") 
                                if(newfile[c]):
                                    with open(pmids_order[c], "a") as g:
                                        g.write("%s" %(complete_text) )

                        for blocks in abs.findall('p'):
                            complete_text=ET.tostring(blocks, encoding="unicode") 
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("%s" %(complete_text) )
            c+=1
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</abstract>")
            c+=1

        # getting sentences
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("<sentences>")
            
            c+=1
        c=0
        for article in root.findall('article'):
            for body in article.findall('body'):
                for section in body.findall('sec'):
                    complete_text=ET.tostring(section, encoding="unicode") 
                    if(newfile[c]):
                        with open(pmids_order[c], "a") as g:
                            g.write("<textgroup >%s</textgroup>" %( complete_text) )

                for blocks in body.findall('p'):
                    complete_text=ET.tostring(blocks, encoding="unicode") 
                    if(newfile[c]):
                        with open(pmids_order[c], "a") as g:
                            g.write("<textgroup >%s</textgroup>" %(complete_text) )

                    """for section2 in section.findall('sec'):
                        section_type2 = ""
                        if('sec-type' in section2.attrib):
                            section_type2 = "gtype='"+section2.attrib['sec-type']+"'"
                        for blocks2 in section2.findall('p'):    
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("<textgroup %s >%s</textgroup>" %(section_type2, blocks2.text) )

                    section_type = ""
                    if('sec-type' in section.attrib):
                        section_type = "gtype='"+section.attrib['sec-type']+"'"
                    for blocks in section.findall('p'):
                        if(newfile[c]):
                            with open(pmids_order[c], "a") as g:
                                g.write("<textgroup %s >%s</textgroup>" %(section_type, blocks.text) )"""
            c+=1
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</sentences>")
            c+=1

        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</article>")
            c+=1

        return new_hits

    def get_pubmedHits_for_pair(self, epi):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
         
        #handle = Entrez.esearch(db='pubmed', sort='relevance', term="epitope ("+epi+")" )
        if(epi.startswith("HLA")):
            handle = Entrez.esearch(db='pubmed', sort='relevance', term=""+epi+" and (sars-cov or covid or virus or viral or coronavirus)" )
        else:
            handle = Entrez.esearch(db='pubmed', sort='relevance', term=""+epi )
        results = Entrez.read(handle)
        hits = results['IdList']
        return hits

    # Store the title, abstract and sentences for each article (pmcid1.xml, ...) 
    def get_hits_on_pubmed_sentences(self, hits, folder):
        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key="4543094c8a41e6aecf9a1431bff42cfac209"
        new_hits=0
        
        id=",".join(hits)
        fetch = Entrez.efetch(db='pubmed',resetmode='xml',id=id,rettype='full')
        with open(folder+'tempFile.xml', 'w') as f:
            f.write(str(fetch.read()).replace("b\\'","").replace("\\'","").replace("\\n","").replace("b\'","").replace("\'",""))
        tree = ET.parse(folder+'tempFile.xml')
        root = tree.getroot()
        
        if( not os.path.isdir(folder+"pubmed_articles") ):
            os.system("mkdir "+folder+"pubmed_articles")
        
        # article file creation
        newfile=[]
        pmids_order=[]
        for article in root.findall('PubmedArticle'):
            for front in article.findall('MedlineCitation'):
                for pmids in front.findall('PMID'):
                    if(not os.path.isfile(folder+"pubmed_articles/"+pmids.text+".xml")):
                        newfile.append(True)
                        pmids_order.append(folder+"pubmed_articles/"+pmids.text+".xml")
                        with open(folder+"pubmed_articles/"+pmids.text+".xml", "a") as g:
                            g.write("<article>")
                    else:
                        newfile.append(False)
                        pmids_order.append(folder+"pubmed_articles/"+pmids.text+".xml")
                        
        #print(pmids_order)
        #print(len(pmids_order), len(newfile))
        # getting title
        c=0
        for article in root.findall('PubmedArticle'):
            for front in article.findall('MedlineCitation'):
                for meta in front.findall('Article'):
                    for title in meta.findall('ArticleTitle'):
                        if(newfile[c]):
                            with open(pmids_order[c], "a") as g:
                                g.write("<title>%s</title>" %(title.text) )
            c+=1

        # getting abstract
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("<abstract>")
            
            c+=1
        c=0
        for article in root.findall('PubmedArticle'):
            for front in article.findall('MedlineCitation'):
                for meta in front.findall('Article'):
                    for abs_ in meta.findall('Abstract'):
                        for sec in abs_.findall('AbstractText'):
                            complete_text=ET.tostring(sec, encoding="unicode") 
                            if(newfile[c]):
                                with open(pmids_order[c], "a") as g:
                                    g.write("%s" %(complete_text) )
            c+=1
        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</abstract>")
            c+=1

        c=0
        for article in pmids_order:
            if(newfile[c]):
                with open(article, "a") as g:
                    g.write("</article>")
            c+=1

        return new_hits

    # evaluate sentences in the text-groups and save those that have the couple of proteins under evaluation

class Nlp_pre_processing_steps:
    ## BEGIN NLP pre processing steps
    def remove_non_ascii(self, words):
        """Remove non-ASCII characters from list of tokenized words"""
        new_words = []
        for word in words:
            new_word = unicodedata.normalize('NFKD', word).encode('ascii', 'ignore').decode('utf-8', 'ignore')
            new_words.append(new_word)
        return new_words

    def to_lowercase(self, words):
        """Convert all characters to lowercase from list of tokenized words"""
        new_words = []
        for word in words:
            new_word = word.lower()
            new_words.append(new_word)
        return new_words

    def remove_punctuation(self, words):
        """Remove punctuation from list of tokenized words"""
        new_words = []
        for word in words:
            new_word = re.sub(r'[^\w\s]', '', word)
            if new_word != '':
                new_words.append(new_word)
        return new_words

    def replace_numbers(self, words):
        """Replace all interger occurrences in list of tokenized words with textual representation"""
        p = inflect.engine()
        new_words = []
        for word in words:
            if word.isdigit():
                new_word = p.number_to_words(word)
                new_words.append(new_word)
            else:
                new_words.append(word)
        return new_words

    def remove_stopwords(self, words):
        """Remove stop words from list of tokenized words"""
        new_words = []
        for word in words:
            if word not in stopwords.words('english'):
                new_words.append(word)
        return new_words

    def stem_words(self, words):
        """Stem words in list of tokenized words"""
        stemmer = LancasterStemmer()
        stems = []
        for word in words:
            stem = stemmer.stem(word)
            stems.append(stem)
        return stems

    def lemmatize_verbs(self, words):
        """Lemmatize verbs in list of tokenized words"""
        lemmatizer = WordNetLemmatizer()
        lemmas = []
        for word in words:
            lemma = lemmatizer.lemmatize(word, pos='v')
            lemmas.append(lemma)
        return lemmas

    def normalize(self, words):
        words = self.remove_non_ascii(words)
        words = self.to_lowercase(words)
        words = self.remove_punctuation(words)
        words = self.replace_numbers(words)
        words = self.remove_stopwords(words)
        return words

    ## END NLP pre processing steps

from lxml import etree
from bs4 import BeautifulSoup
import re
class Text_processing:
    def cleanhtml(self, raw_html):
        cleanr = re.compile('<.*?>')
        cleantext = re.sub(cleanr, '', raw_html)
        return cleantext

    def organize_sentences_with_bait(self, file_evaluation_pairs, folder): 
        p=Nlp_pre_processing_steps()
        try:
            tokenizer = nltk.data.load('tokenizers/punkt/PY3/english.pickle')
        except:
            import nltk
            nltk.download('punkt')
            tokenizer = nltk.data.load('tokenizers/punkt/PY3/english.pickle')
            
        executable_path = os.path.join(".", "geniatagger-3.0.2", "geniatagger")
        tagger = GENIATagger(executable_path)
        
        lt = Literature_hits()

        parser = etree.XMLParser(recover=True)

        pairs=[]

        if(os.path.isdir(folder+"processed_sentences")):
            # Loading processed pairs
            for f in os.listdir(folder+"processed_sentences"):
                if(f.startswith("scs_")):
                    pr=f.replace("scs_","").replace(".tsv","")
                    pairs.append(pr)
        else:
            os.system("mkdir "+folder+"processed_sentences")

        if(os.path.isdir(folder+"pmc_articles") or os.path.isdir(folder+"pubmed_articles")):
            f=open(folder+file_evaluation_pairs,"r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(int(l[1])>0 or int(l[3])>0):
                    p=l[0]
                    
                    if( not p in pairs ):
                        p=p.lower()

                        identifier_result=p

                        g=open(folder+"processed_sentences/scs_"+identifier_result+".tsv","w")
                        g.write("id_pmc\tsentence location\thighlighted phrase\tprotein entities found\tnumber of protein entities\n")
                        g.close()
                        
                        type_="pmc"
                        pmids=[]
                        if(l[4]!=""):
                            pmids=l[4].split(",")
                        for a in pmids:
                            self.process_paper(folder, a, type_, identifier_result, tokenizer, tagger, lt, parser, p)

                        type_="pubmed"
                        pmids=[]
                        if(l[2]!=""):
                            pmids=l[2].split(",")
                        for a in pmids:
                            self.process_paper(folder, a, type_, identifier_result, tokenizer, tagger, lt, parser, p)

                        if(os.path.getsize(folder+"processed_sentences/scs_"+identifier_result+".tsv")==94):
                            os.system("rm "+folder+"processed_sentences/scs_"+identifier_result+".tsv")
                            
            f.close()
        else:
            print("Error: You have to run step 2 first because the pmc/pubmed papers found in step 1 are not processed and stored")

    def process_paper(self, folder, a, type_, identifier_result, tokenizer, tagger, lt, parser, pair):
        if( os.path.isfile(folder+type_+'_articles/'+a+'.xml') ):
            #print(p1, p2, 'pmc_articles/'+a+'.xml')
            #try:
            # fixing xml file containing possible undesirable characters for xml
            root = etree.parse(folder+type_+'_articles/'+a+'.xml', parser=parser)
            fixed_xml=str(etree.tostring(root)).replace("b\'","").replace("'","")
            h=open(folder+type_+'_articles/'+a+'.xml',"w")
            h.write(fixed_xml)
            h.close()

            tree = ET.parse(folder+type_+'_articles/'+a+'.xml')
            root = tree.getroot()

            # Removing unnecessary tags and their contents
            txt=ET.tostring(root, encoding="unicode")
            txt=txt.replace("\\n","")
            soup = BeautifulSoup(txt, "xml")
            for tag in soup.find_all(['supplementary-material','xref','sup','fig']):
                tag.decompose()
            h=open(folder+type_+'_articles/'+a+'.xml',"w")
            h.write(str(soup.prettify()))
            h.close()
    
            sentence_groups={}

            sentences=[]
            for tg in root.findall('abstract'):
                txt=ET.tostring(tg, encoding="unicode") 
                txt=txt.replace("\\n","").replace("\\","")
                aux=[]
                for t in txt.split(" "):
                    if(t.find("HLA")!=-1):
                        t=t.replace(":","-")
                    aux.append(t)
                txt=' '.join(aux)
                temp=tokenizer.tokenize(self.cleanhtml(txt))
                
                #print(temp)
                for s in temp:
                    tokenized = nltk.word_tokenize(s.replace("\n",""))
                    #normalized = p.normalize(tokenized)
                    sentences.append(" ".join(tokenized).replace("/"," / "))
            sentence_groups["abstract"]=sentences

            sentences=[]
            for scs in root.findall('sentences'):
                for tg in scs.findall('textgroup'):
                    txt=ET.tostring(tg, encoding="unicode") 
                    txt=txt.replace("\\n","").replace("\\","")
                    aux=[]
                    for t in txt.split(" "):
                        if(t.find("HLA")!=-1):
                            t=t.replace(":","-")
                        aux.append(t)
                    txt=' '.join(aux)
                    
                    temp=tokenizer.tokenize(self.cleanhtml(txt))
                    for s in temp:
                        tokenized = nltk.word_tokenize(s.replace("\n",""))
                        #normalized = p.normalize(tokenized)
                        sentences.append(" ".join(tokenized).replace("/"," / "))                      
            sentence_groups["body"]=sentences

            for idst in sentence_groups.keys():
                for s in sentence_groups[idst]:
                    #print(s)
                    count_forbidden_expressions=0
                    
                    if(count_forbidden_expressions==0):
                        c1=0
                        c2=0
                        
                        w=[]
                        pos=[]
                        entities=[]
                        for word, base_form, pos_tag, chunk, named_entity in tagger.tag(s):
                            if(named_entity.find("protein")!=-1):
                                entities.append(word)

                            p1=pair
                            p1=p1.replace(":","-").replace("*","").lower()
                                
                            word=word.lower()
                            
                            condition=(word.find(p1)!=-1)
                            if((named_entity.find("protein")!=-1 or named_entity.find("O")!=-1 or named_entity.find("DNA")!=-1) and condition ):
                                c1+=1
                            
                            w.append(base_form)
                            #pos.append(pos_tag)
                        
                        if(lt!=None):
                            main_condition = (c1>0)

                        if( main_condition ):
                            with open(folder+"processed_sentences/scs_"+identifier_result+".tsv","a") as fg:
                                #fg.write("%s\t%s\t%s\t%s\n" %(a, str(w), str(pos), str(entities) ) )
                                fg.write("%s\t%s\t%s\t%s\t%i\n" %(a, idst, "|".join(w), "|".join(entities), len(entities) ) )
            #except:
            #    print("not ok")
                
class Running_config:

    def run_step1_mode1(self, file_pairs, folder):
        lt=Literature_hits()

        if( not os.path.isdir(folder+"pmc_articles") ):
            os.system("mkdir "+folder+"pmc_articles/*")

        f=open(folder+"literature_evaluation_pairs.tsv","w")
        f.close()
        c=0
        f=open(folder+file_pairs,"r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").split(",")
                if(len(l)>1):
                	l[0]=l[1]
                hits_pmc=lt.get_pmcHits_for_pair(l[0])
                hits_pubmed=lt.get_pubmedHits_for_pair(l[0])
                with open(folder+"literature_evaluation_pairs.tsv","a") as g:
                    g.write("%s\t%i\t%s\t%i\t%s\n" %(l[0], len(hits_pubmed), (",".join(hits_pubmed)), len(hits_pmc), (",".join(hits_pmc)) ) )
            c+=1
        f.close()

    def run_step2_mode1(self, file_evaluation_pairs, folder):
        lit=Literature_hits()
        
        print("Step 2 - Getting articles content...")
        f=open(folder+file_evaluation_pairs,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if( int(l[3])>0):
                hits=l[4].split(",")
                hits_sentences=lit.get_hits_on_pmc_sentences(hits, folder)

            if( int(l[1])>0):
                hits=l[2].split(",")
                hits_sentences=lit.get_hits_on_pubmed_sentences(hits, folder)

        f.close()

    def run_step3_mode1(self, file_evaluation_pairs, folder):
        t=Text_processing()
        t.organize_sentences_with_bait(file_evaluation_pairs, folder)


    def run(self, args):
        if(args.folder!="" and os.path.isdir(args.folder)):
            run=0
            if(args.running_type=="" ):
                run=0
            else:
                if(args.running_type in [0,1,2,3]):
                    run=args.running_type
                else:
                    print("Error: invalid choice")

            if(run==0 or run==1):
                if(args.file_epitopes==""):
                    print("Error: you have to give the file with epitopes")
                else:
                    print("Running step 1")
                    self.run_step1_mode1(args.file_epitopes, args.folder)
                    if(run==0):
                        print("Running step 2")
                        self.run_step2_mode1("literature_evaluation_pairs.tsv", args.folder)
                        print("Running step 3")
                        self.run_step3_mode1("literature_evaluation_pairs.tsv", args.folder)
            else:
                if(args.file_evaluation==""):
                    print("Error: you have to give the evaluation file exported in step 1")
                else:
                    if(args.file_evaluation==""):
                        print("Error: you have to give the file with pairs")
                    else:
                        if(run==2):
                            self.run_step2_mode1(args.file_evaluation, args.folder)
                        else:
                            self.run_step3_mode1(args.file_evaluation, args.folder)

        else:
            print("Error: You have to specify a valid folder to store files")


import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=' EpiMiner - Literature pipeline to find epitopes on pubmed articles', formatter_class=RawTextHelpFormatter)
parser.add_argument("-fo", "--folder", action="store", help="(For both modes) Folder to store the files (use the folder where the required files can be found, ex.: /home/user/experiment/ )\n")

parser.add_argument("-rt", "--running_type", action="store", help="0 (default) - Run all steps\n\
1 - Run step 1 (Get mentions of epitopes in PMC articles)\n\
2 - Run step 2 (Get the PMC or Pubmed files, clean and store them)\n\
3 - Run step 3 (Get the exact sentences where the epitopes were found)", type=int)
parser.add_argument("-fp", "--file_epitopes", action="store", help="File with the epitopes (epitope sequence in first column in tsv format without header)")
parser.add_argument("-fe", "--file_evaluation", action="store", help="File exported after step 1 execution in tsv format\n")

args = parser.parse_args()
r=Running_config()
r.run(args)               
