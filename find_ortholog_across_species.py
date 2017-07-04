import glob

def covert_ortholog_list(reference):
    dict1={}
    dict3={}
    taxa=[]
    out=open('species_across_ortholog.txt','w')
    out.write(reference.split("_exon_sorted_2.fas")[0]+'\t')
    homo_list=[]

    
    for filename in glob.glob('*_candidate.txt'):
        dict3[filename]=[]
        for record in open(filename):
            dict3[filename].append(record.split('\t')[0])
            if record.split('\t')[0] not in homo_list:
                homo_list.append(record.split('\t')[0])

                
    for filename in glob.glob('*_candidate.txt'):
            out.write(filename.split('_two_way_candidate.txt')[0].split(reference.split("exon_sorted_2.fas")[0])[1]+'\t')
            for homo in homo_list:
                if homo not in dict3[filename]:
                        if dict1.has_key(homo):
                            dict1[homo].append('None')
                        else:
                            dict1.setdefault(homo,['None'])
            for record in open(filename):
                        if int(record.split('\t')[-1].strip())==1:
                            newrecord=str(record.split('\t')[1])+'|Plus'
                            if dict1.has_key(record.split('\t')[0]):
                                dict1[record.split('\t')[0]].append(newrecord)
                            else:
                                dict1.setdefault(record.split('\t')[0],[newrecord])
                        else:
                            newrecord=str(record.split('\t')[1])+'|Minus'
                            if dict1.has_key(record.split('\t')[0]):
                                dict1[record.split('\t')[0]].append(newrecord)
                            else:
                                dict1.setdefault(record.split('\t')[0],[newrecord])                       
     

    dict2={}
    out.write('\n')
    for key in dict1:
        if (len(dict1[key])-dict1[key].count('None'))>=2:
            if dict2.has_key(key.split('|')[2]):
                dict2[key.split('|')[2]].append((key,dict1[key]))
            else:
                dict2[key.split('|')[2]]=[(key,dict1[key])]

    for chromosome in sorted(dict2.iteritems(),key=lambda e:e[0]):
        for record in sorted(chromosome[1],key=lambda e:int(e[0].split('|')[4])):
            out.write(str(record).replace('(','').replace(')','').replace('\'','').replace(' ','').replace(',','\t').replace('[','').replace(']','')+'\n')
    out.close()
        
