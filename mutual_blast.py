import multiprocessing,subprocess,os

def blast(query, db, evalue="1e-5",num_threads=multiprocessing.cpu_count()-1):
    result_dir=os.path.join(os.path.dirname(query), "xml")
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    out=os.path.join(result_dir, os.path.basename(query)+"=}"+os.path.basename(db)+'.xml')
    command="blastn -query %s -db %s -outfmt 5 -max_target_seqs 20 -evalue 1e-5 -num_threads 6 -out %s" %(query,db,out)
    os.system(command)
