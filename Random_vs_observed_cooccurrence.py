# -*- coding: utf-8 -*-
"""
@author: Feng Ju
@email: richieju520@gmail.com

"""

print 'This script is written for adding attribute to each NODE in GML files (readable by Gehpi)'
print 'This script also calculate the random and observed incidences of co-occurrence between differnt types of nodes!'
print 'The map file is a tab-delimited file with node ID (col 1) and type of node (col 2)
print 'The gml file is the gml-format network file generated from the R scripts'

while True:
    Parameters=raw_input("Enter parameters [map file](.map) and [GML file] sepeated by Space: ")
    try:
        P1=Parameters.strip().split(' ')[0]
        P2=Parameters.strip().split(' ')[1]
        break
    except:
         print 'errors: invalid input format or not enough input !'
         continue

a={}
for line in open(P1,'r'):
    try:
        a[line.rstrip().split('\t')[0]]=line.rstrip().split('\t')[1]
    except:
        print 'Pls provdied a tab-delimited txt file!'

f=open(P2.replace('.gml','') +'.modified.gml','w')
f1=open(P1 +'.map','w')

b=[]

i=0
for line in open(P2,'r'):
    if 'label' not in line:
        f.write(line)
    else:
        i+=1
        name = line.strip().split('"')[1]
        f.write(line)
        try:
            f1.write(name+'\t'+str(a[name])+'\n')
            f.write('    '+'order'+' '+str(a[name])+'\n')
        except KeyError:
            print name
print i,'items were added into the node!'
print 'OK, work finished!'
f.close()
f1.close()

X1 = P1 +'.map'
X2 = P2.replace('.gml','') +'.modified.gml'

f=open(X2.replace('.gml','')+'_Observed_VS_Random.xls','w')
f1=open(X2.replace('.gml','')+'_edge_properties.xls','w')
f2=open(X2.replace('.gml','')+'_node_properties.xls','w')

a1,a2={},{}
lis=[]
for line in open(X1,'r'):
    a1[line.strip().split('\t')[0]]=line.strip().split('\t')[1]
    lis.append(line.strip().split('\t')[1])
print len(a1), 'Node-affiliation pairs!'

dic1={}
lis_U=list(set(lis))
for item in lis_U:
    dic1[item]=lis.count(item)
      
f2.write('\t'.join(['id','nmae','phylum','degree'])+'\n')
b,c,d,e,g = [],[],[],[],[]
for line in open(X2,'r'):
    if '    id ' in line:
        f2.write(line.strip().split(' ')[1]+'\t')
        b.append(line.strip().split(' ')[1])
    elif '    name' in line:
        f2.write(line.strip().split(' ')[1].replace('"','')+'\t')
        c.append(line.strip().split(' ')[1].replace('"',''))
        f2.write(a1[line.strip().split(' ')[1].replace('"','')]+'\t')
    elif '    degree' in line:
        f2.write(line.strip().split(' ')[1]+'\n')
    elif 'source' in line:
        d.append(line.strip().split(' ')[1])
    elif 'target' in line:
        e.append(line.strip().split(' ')[1])
    elif 'weight' in line:
        g.append(line.strip().split(' ')[1])  
    else:
        continue

print len(b), len(c),'nodes!'
print len(d), len(e),'edges!'
    
for i in range(len(b)):
    a2[b[i]]=c[i]

j,k=0,0
n=[]
p=[]

f1.write('Source'+'\t'+'Phylum'+'\t'+'Target'+'\t'+'Phylum'+'\t'+'Weight'+'\n')
for m in range(len(g)):
    if a1[a2[d[m]]]==a1[a2[e[m]]]:
        f1.write(str(a2[d[m]])+'\t'+str(a1[a2[d[m]]])+'\t'+str(a2[e[m]])+'\t'+str(a1[a2[e[m]]])+'\t'+str(g[m])+'\n')
        j+=1
        n.append(a1[a2[d[m]]]+'__'+a1[a2[e[m]]])
    else:
        f1.write(str(a2[d[m]])+'\t'+str(a1[a2[d[m]]])+'\t'+str(a2[e[m]])+'\t'+str(a1[a2[e[m]]])+'\t'+str(g[m])+'\n')
        k+=1
        p.append(a1[a2[d[m]]]+'__'+a1[a2[e[m]]])

print j, 'Internal-type cooccurence!'
print k, 'External-type cooccurence!'

n1=list(set(n))
n1.sort()
dic2={}
for item in n1:
    dic2[item]=str(n.count(item))

p1=list(set(p))
p1.sort()
for item in p1:
    i1=p.count(item)
    i2=p.count(item.split('__')[1]+'__'+item.split('__')[0])
    dic2[item]=i1+i2
    if i2!=0:
        p1.remove(item.split('__')[1]+'__'+item.split('__')[0])
    else:
        continue

f.write(str(len(c))+' '+'nodes and'+' '+ str(len(e))+' '+'edges in the network!'+'\n')
f.write(str(j)+' '+'edges with internal-type cooccurence!'+'\n')
f.write(str(k)+' '+'edges with external-type cooccurence!'+'\n')

f.write('N1__N2'+'\t'+'N1-freq'+'\t'+'N2-freq'+'\t'+'Edges'+'\t'+'Random'+'\t'+'Observed'+'\t'+'O/R-ratio'+'\n')
for key in n1:
    i1=dic1[key.split('__')[0]]
    i2=dic1[key.split('__')[1]]
    if key.split('__')[0]==key.split('__')[1]:
        random=100*float(i1*(i2-1))/(len(c)*(len(c)-1))
    else:
        random=2*100*float(i1*i2)/(len(c)*(len(c)-1))
    observed=100*float(dic2[key])/len(e)
    ratio=observed/random
    f.write('\t'.join([key, str(i1), str(i2), str(dic2[key]), str(random), str(observed), str(ratio)])+'\n')
for key in p1:
    i1=dic1[key.split('__')[0]]
    i2=dic1[key.split('__')[1]]
    if key.split('__')[0]==key.split('__')[1]:
        random=100*float(i1*(i2-1))/(len(c)*(len(c)-1))
    else:
        random=2*100*float(i1*i2)/(len(c)*(len(c)-1))
    observed=100*float(dic2[key])/len(e)
    ratio=observed/random
    f.write('\t'.join([key, str(i1), str(i2), str(dic2[key]), str(random), str(observed), str(ratio)])+'\n')

f.write('\n')
f.write('Node-affiliation'+'\t'+'Nodes_freq'+'\n')

list_s = sorted(dic1, key=dic1.__getitem__, reverse = True)
for key in list_s:
    f.write(key+'\t'+str(dic1[key])+'\n')
print 'OK, finished!'
