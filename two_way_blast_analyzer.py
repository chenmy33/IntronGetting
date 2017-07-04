import glob,os





def find_orthlog_each_species(filename,score_ratio=1.5,exon_length_vaiation=1.5):
        path_out='Exon_length_uncertainty_ortholog'
        if not os.path.exists(path_out):
                os.makedirs(path_out)
        c=0
        query=[]
        ortholog=''    
        dict1={}
        outp1=filename.replace('_result.txt','_ortholog.txt')
        outp2='Uncertain_'+filename.replace('_result.txt','_ortholog.txt')
        out1=open(outp1,'w')
        out2=open(os.path.join(path_out,outp2),'w')
        for i in open(filename):
            if c==0:
                c=1
            else:
                        if  i.split('\t')[0].strip() not in query:
                            suspects=[]
                            dict1[i.split('\t')[0]]=suspects
                            query.append(i.split('\t')[0].strip())                    
                            suspects.append(i.split('\t'))
                        else:
                            suspects.append(i.split('\t'))
                

        for key in dict1:
            if len(dict1[key])>1:
               list1s=dict1[key]
               list1s=sorted(list1s,key=lambda list1:list1[-2])
               if float(list1s[-1][-2])*1.0/float(list1s[-2][-2])>score_ratio:
                       ortholog = list1s.pop()
                       if (float(ortholog[9])/float(ortholog[10])>exon_length_vaiation) or ((float(ortholog[10])/float(ortholog[9])>exon_length_vaiation)):
                           out2.write('\t'.join(ortholog))
                       else:                                                 
                           out1.write('\t'.join(ortholog))
            if len(dict1[key])==1:
                       ortholog =dict1[key][0]
                       if (float(ortholog[9])/float(ortholog[10])>exon_length_vaiation) or ((float(ortholog[10])/float(ortholog[9])>exon_length_vaiation)):
                           out2.write('\t'.join(ortholog))
                       else:                                                 
                           out1.write('\t'.join(ortholog))
        out1.close()
        out2.close()
        return filename.replace('_result.txt','_ortholog.txt')

def find_two_way_blast_ortholog(filename,ref_name):
        out_path1='two_way_blast_orthologs'
        out_path2='one_way_blast_orthologs'
        if not os.path.exists(out_path1):
                os.makedirs(out_path1)        
        if not os.path.exists(out_path2):
                os.makedirs(out_path2)
        
        outp1=filename.split('_exon_sorted')[0]+'_'+filename.split("=}")[1].split('_exon_sorted')[0]+'_two_way_candidate.txt'
        outp2=filename.split('_exon_sorted')[0]+'_'+filename.split("=}")[1].split('_exon_sorted')[0]+'_one_way_candidate.txt'
        out1=open(os.path.join(out_path1,outp1),'w')
        out2=open(os.path.join(out_path2,outp2),'w')
        dict_f={}
        dict_r={}
        rname= filename.split("=}")[1].strip("_ortholog.txt")+"=}"+ref_name+'_ortholog.txt'
        r_handle=open(rname)
        for line in r_handle:
                dict_r[line.split()[0]]=line.split()[1]
        f_handle=open(filename)
        for line in f_handle:
                if (dict_r.has_key(line.split()[1])) and (dict_r[line.split()[1]]== line.split()[0]):
                    out1.write(line.split()[0]+'\t'+line.split()[1]+'\t'+line.split()[-1]+'\n')
                    del dict_r[line.split()[1]]
                else:
                    out2.write(line)
        for key in dict_r:
                out2.write(str(key)+'\t'+dict_r[key]+'\n')
        out1.close()
        out2.close()       
       
        
            
            
