a=[]
c=0
f=open("data/input.csv")
for line in f:
    if(c>0):
        l=line.replace("\n","").split(",")
        a.append(l[1])
    c+=1
 
f.close()

c=0
f=open("data/literature_evaluation_pairs.tsv")
for line in f:
    l=line.replace("\n","").split("\t")
    with open("_literature_evaluation_pairs.tsv","a") as gf:
        gf.write(a[c]+"\t"+("\t".join(l[1:]))+"\n")
        
    c+=1
