from Bio import SeqIO
from Bio.Seq import Seq
import os,glob
def reverse(i,seq):
        if i==1:
            return seq
        if i==-1:
            temp=Seq(seq).reverse_complement()
            return str(temp)
def retrieve_seqs(genome_dir,out_dir):
    exon_dict={}#sequeces pool
    chromosome_dict={}#chromosome pool:{taxa:{chromosome:gene_location}}
    a=0
    b=0
    exon_file='exon_pairs.txt'
    for line in open(exon_file):
        if b==0:
            b=1
            taxon=line.split()
            for i in range(len(taxon)):
                chromosome_dict[taxon[i]]={}

        else:
            for i in range(len(line.split())):
                if line.split()[i]!= 'None&None':
                    if chromosome_dict[taxon[i]].has_key(line.split()[i].split('|')[2]):
                            chromosome_dict[taxon[i]][line.split()[i].split('|')[2]].append(line.split()[i])
                    else:
                       chromosome_dict[taxon[i]][line.split()[i].split('|')[2]]= [line.split()[i]]
    for filename in glob.glob(os.path.join(genome_dir,'*.fa')):
        records=SeqIO.parse(filename,'fasta')
        for record in records:
            a=0
            count=0
            chromosome=str(record.description).split()[0]
            taxa_name= os.path.basename(filename).split('.')[0]
            if taxa_name not in chromosome_dict:
                break
            else:
                if chromosome not in chromosome_dict[taxa_name].keys():
                    continue
                if count == len(chromosome_dict[taxa_name].keys()):
                    del chromosome_dict[taxa_name]
                    break
                else:
                    record_seq=str(record.seq)
                    for key in chromosome_dict:
                        if a!=0:
                            break
                        else:
                            if key== taxa_name:
                                for sub_key in chromosome_dict[key]:
                                    if a!=0:
                                        break
                                    else:
                                        if sub_key==chromosome:
                                            count+=1
                                            for genename in chromosome_dict[taxa_name][sub_key]:
                                                if genename!= 'None&None':
                                                    exon1_start=int(genename.split('&')[0].split('|')[4])
                                                    exon2_start=int(genename.split('&')[1].split('|')[4])
                                                    if exon1_start<exon2_start:
                                                        seq_Start=int(genename.split('&')[0].split('|')[4])-1
                                                        seq_End=int(genename.split('&')[1].split('|')[5])
                                                        taxa_seq=record_seq[seq_Start:seq_End]
                                                    else:
                                                        seq_Start=int(genename.split('&')[1].split('|')[4])-1
                                                        seq_End=int(genename.split('&')[0].split('|')[5])
                                                        taxa_seq=record_seq[seq_Start:seq_End]
                                                    if exon_dict.has_key(taxa_name):
                                                        exon_dict[taxa_name].append((genename,taxa_seq))
                                                    else:
                                                        exon_dict.setdefault(taxa_name,[(genename,taxa_seq)])
                                                else:
                                                    continue

                                            a=1
                            
                     

                                    
    taxon=[]
    a=0
    seq_dict={}


    outp1=open('reference_missing_exon.txt','w')

    for line in open(exon_file):
        if a==0:
            a=1
            taxon=line.split()
        else:
            if line.split()[0]!= 'None&None':
                reference_direction= int(line.split()[0].split('&')[1].split('|')[-1])
                exon_info= line.split()[0]
                genename=line.split()[0].split('&')[0].split('|')[0]+'_'+line.split()[0].split('&')[0].split('|')[1]+'_'+line.split()[0].split('&')[1].split('|')[1].split('_')[1]
                for info in exon_dict[taxon[0]]:
    ##                        print info[0]
                            if info[0]== exon_info:
                                seq= info[1]
                                real_seq=reverse(reference_direction,seq)
                                seq_dict.setdefault(genename,[(taxon[0],real_seq)])
                                break

                        
                for i in range(1,len(line.split())): 
                    if line.split()[i]!= 'None&None':
                        sself_direction= [int(line.split()[i].split('&')[1].split('|')[-2]),str(line.split()[i].split('&')[1].split('|')[-1])]
                        if (sself_direction[0]==reference_direction) and (sself_direction[1]=='Plus'):
                                self_direction=reference_direction
                        if (sself_direction[0]==reference_direction) and (sself_direction[1]=='Minus'):
                                self_direction=reference_direction*(-1)
                        if (sself_direction[0]==reference_direction*(-1)) and (sself_direction[1]=='Plus'):
                                self_direction=reference_direction*(-1)
                        if (sself_direction[0]==reference_direction*(-1)) and (sself_direction[1]=='Minus'):
                                self_direction=reference_direction
                        exon_info= line.split()[i]
                        for info in exon_dict[taxon[i]]:

                                if info[0]== exon_info:
                                    seq2= info[1]
                                    real_seq2=reverse(self_direction,seq2)
                                    if seq_dict.has_key(genename):
                                        seq_dict[genename].append((taxon[i],real_seq2))
                                    break

            else:
                outp1.write(line)
                
    for key in seq_dict:
        filename=str(key)+'.fas'
        outp=open(os.path.join(out_dir,filename),'w')
        for seq in seq_dict[key]:
            outp.write('>'+str(seq[0])+'\n'+str(seq[1])+'\n')
        outp.close()
            
    

                    


            
