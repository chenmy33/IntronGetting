import glob,os,sys,sort_exon,mutual_blast,blast_result_filter,two_way_blast_analyzer, find_ortholog_across_species,intron_length_variation_missing_control,retrieve_seq
def get_EPIC(exon_dir,genome_dir,reference_name,cpu_threads):

    
    # Step1: Sort exons of each species according to their genomic position.
    #head information format:Gene stable ID|Gene name|Chromosome/scaffold name|Exon stable ID|Genomic coding start|Genomic coding end|strand
    out_dir=os.path.join(os.getcwd(),"Results")
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    for filename in glob.glob(os.path.join(exon_dir,"*_exon")):
        sorting_name=filename[:-5]+'_exon_sorted_1.fas'
        sorted_name=filename[:-5]+'_exon_sorted_2.fas'
        sort_exon.sort(filename,sorting_name,5)
        sort_exon.sort(sorting_name,sorted_name,4)

        
    ### Step2: Find orthologous exons using mutual best hit (MBH) in BLAST, selecting one species as the reference species.e.g. we selected Homo_sapiens as reference. Note: Add BLAST to your PATH environment variable
    reference=reference_name+"_exon_sorted_2.fas"
    for filename in glob.glob(os.path.join(exon_dir,'*exon_sorted_2.fas')):
        if not os.path.exists(filename+".nhr"):
            os.system("makeblastdb -in %s -dbtype nucl" %(os.path.abspath(filename))) #creat BLAST databases
        else:
            pass
        
    for filename in glob.glob(os.path.join(exon_dir,'*exon_sorted_2.fas')):
        if filename != os.path.join(exon_dir,reference):
            mutual_blast.blast(os.path.join(exon_dir,reference),os.path.abspath(filename),cpu_threads)
            mutual_blast.blast(os.path.abspath(filename),os.path.join(exon_dir,reference),cpu_threads)


    ### Step3: Covert BLAST results to tabular and filter low quality results. You can set coverage and identity cutoff to remove low quality results.
    os.chdir(os.path.join(exon_dir,"xml"))
    coverage=0.6
    identity=0.7
    for filename in glob.glob("*.xml"):
        blast_result_filter.blast_result_filter(filename,identity,coverage)

        
    # Step4: Analyze result files and find two-way BLAST ortholog exons.
    os.chdir(os.path.join(exon_dir,"xml","Coverage%.2f_Identity%.2f" %(coverage,identity)))
    first_second_score_ratio=1.5
    exon_length_vaiation=1.5
    for filename in glob.glob('*_result.txt'):
        p=two_way_blast_analyzer.find_orthlog_each_species(filename,first_second_score_ratio,exon_length_vaiation)
    for filename in glob.glob('*_ortholog.txt'):
        if filename.startswith(reference):
            two_way_blast_analyzer.find_two_way_blast_ortholog(filename,reference)


    ### Step5: Find ortholog exons across species.            
    os.chdir(os.path.join(exon_dir,"xml","Coverage%.2f_Identity%.2f" %(coverage,identity),"two_way_blast_orthologs"))
    find_ortholog_across_species.covert_ortholog_list(reference)


    ### Step6: Filter candidates by taxon coverage, intron_length, intron/exon length ratio.
    taxon_coverage=0.8
    intron_length=10000
    intron_length_ratio=1
    intron_length_variation_missing_control.candidates_filter(taxon_coverage,intron_length,intron_length_ratio)


    ### Step7: Retrieve sequences.
    ### The step requires large memory as it recruits genome sequences.
    retrieve_seq.retrieve_seqs(genome_dir,out_dir)

if __name__ == "__main__":
    if len(sys.argv) != 5:
	print "usage: python get_EPIC.py exon_dir genome_dir reference_name cpu_threads"
	sys.exit()
    exon_dir = sys.argv[1]
    reference_name = sys.argv[3]
    genome_dir = sys.argv[2]
    cpu_threads=sys.argv[4]
    get_EPIC(exon_dir,genome_dir,reference_name,cpu_threads)
