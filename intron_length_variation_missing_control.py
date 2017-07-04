def distance(exon_up,exon_end):
    if int(exon_end.split('|')[4]) <int(exon_up.split('|')[4]):
        distance=abs(int(exon_end.split('|')[5])-int(exon_up.split('|')[4]))
    else:
        distance=abs(int(exon_end.split('|')[4])-int(exon_up.split('|')[5]))
    return distance
def gname_ctrl(last_line,line):
    name1= last_line.split()[0].split('|')[0]
    name2= line.split()[0].split('|')[0]
    if (name1==name2):
        return True
    else:
        return False
def intron_length_ctrl(last_line,line,intron_length_ratio):
    intron_length=abs(int(line.split()[0].split('|')[4])-int(last_line.split()[0].split('|')[5]))
    exon_length=abs(int(line.split()[0].split('|')[5])-int(last_line.split()[0].split('|')[4]))-intron_length
    if exon_length*1.0/intron_length<=intron_length_ratio:
        return True
    else:
        return False
        
def miss_ratio(last_line,line,missing_ratio):
    taxa_num1=len(last_line.split())
    taxa_num2=len(line.split())
    miss1=last_line.split().count('None')
    miss2=line.split().count('None')
    if (float(miss1)/taxa_num1<missing_ratio) and (float(miss2)/taxa_num2<missing_ratio):
        return True
    else:
        return False



def candidates_filter(taxon_coverage,intron_length,intron_length_ratio):
    seqlen_variation_ratio=0.5
    missing_ratio=1-taxon_coverage
    c=0    
    last_line=''
    out=open('exon_pairs.txt','w')
    for line in open('species_across_ortholog.txt'):
        if c==0:
            c=1
            out.write(line)
        else:
            if last_line=='':
                last_line=line
            else:
                var=[]
                p=[]
                last=[]
                now=[]
                if gname_ctrl(last_line,line) ==True:
                    if miss_ratio(last_line,line,missing_ratio)==True:
                        if intron_length_ctrl(last_line,line,intron_length_ratio)==True:
                            for i in range(len(line.split())):
                                if (line.split()[i]!='None') and (last_line.split()[i]!='None'):
                                    last.append(last_line.split()[i])
                                    now.append(line.split()[i])
                                    var.append(distance(last_line.split()[i],line.split()[i]))
                                else:
                                    last.append('None')
                                    now.append('None')
                                    var.append('None')
                            t=0
                            d=0
                            a=0
                            for i in range(len(var)):
                                if var[i]!='None':
                                    t+= var[i]
                                    d+=1
                            if d!=0:
                                m=float(t)/d                            
                                for i in range(len(var)):
                                    if var[i]!='None':
                                        if abs(float(var[i])-m)/m > seqlen_variation_ratio:
                                            last[i]='None'
                                            now[i]='None'
                                        else:
                                            a+=1
                                if a>=3 and m<intron_length:
                                    for i in range(len(now)):
                                        out.write(str(last[i])+'&'+str(now[i])+'\t')
                                    out.write('\n')
                                    last_line=line
                                else:
                                    last_line=line
                            else:
                                last_line=line
                        else:
                            last_line=line
                    else:
                        last_line=line
                else:
                    last_line=line
    out.close()
                            
                            
            

