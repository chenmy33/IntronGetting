from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
import glob,os



def blast_result_filter(filename,iden_value=0.7,cover_value=0.6):
    path_out='Coverage%.2f_Identity%.2f' %(cover_value,iden_value)
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    out=open(path_out+'\\'+filename.replace('.xml','_result.txt'),'w')
    out.write('query_id	sbjct_id	query_start	query_end	sbjc_start	sbjct_end	query_seq	sbjct_seq	identities	query_seq_length	sbject_seq_length	align_length\tcoverage	score	frame\n')
    f_blast_records = NCBIXML.parse(open(filename))
    for f_blast_record in f_blast_records:
        for alignment in f_blast_record.alignments:
                for HSP in alignment.hsps:
                    identities=HSP.identities*1.0/HSP.align_length
                    coverage=HSP.align_length*1.0/f_blast_record.query_length
                    if identities>iden_value and coverage>cover_value:                    
                            out.write(str(f_blast_record.query)+'\t'+str(alignment.hit_def)+'\t'+str(HSP.query_start)+'\t'+str(HSP.query_end)\
                                      +'\t'+str(HSP.sbjct_start)+'\t'+str(HSP.sbjct_end)+'\t'+str(HSP.query)+'\t'+str(HSP.sbjct)+'\t'+str(identities)\
                                      +'\t'+str(f_blast_record.query_length)+'\t'+str(alignment.length)+'\t'+str(HSP.align_length)+'\t'+str(coverage)+'\t'+str(HSP.score)+'\t'+str(str(HSP.frame).split(',')[1][:-1])+'\n')
    out.close()

