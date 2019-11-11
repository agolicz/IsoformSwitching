import sys
import math

sl=int(sys.argv[2])
cm=[]

for x in range(sl):
    cm.append({"A" : 0, "T" : 0, "G" : 0, "C" : 0})

for l in open(sys.argv[1]):
    l_arr=l.rstrip().split("\t")
    l1=l_arr[0].rstrip().upper()
    if(len(l1.strip('AGCT')) > 0):
        continue
    l1l=list(l1)
    for i in range(len(l1l)):
        cl=l1l[i]
        cm[i][cl]=cm[i][cl]+1

for l in open(sys.argv[1]):
    l_arr=l.rstrip().split("\t")
    l1=l_arr[0].rstrip().upper()
    if(len(l1.strip('AGCT')) > 0):
        continue
    l1l=list(l1)
    l1s=[]
    for i in range(len(l1l)):
        cl=l1l[i]
        n=math.log((float(cm[i][cl])/(cm[i]['A']+cm[i]['T']+cm[i]['G']+cm[i]['C']))/0.25)
        l1s.append(n)
    scr=sum(l1s)
    print l.rstrip()+"\t"+str(scr)
