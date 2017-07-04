from Bio import SeqIO
import glob
def sort(exonsFasta, sortedFasta,sortpara):   
    exonsDict = {}
    for record in SeqIO.parse(exonsFasta,"fasta"):
        try:
            chromosome = record.id.split("|")[2] #chromosome information
            if chromosome not in exonsDict:
                exonsDict.update({chromosome:[record]})
            else:
                exonsDict[chromosome].append(record)
        except IndexError:
            pass
##    print exonsDict
    with open(sortedFasta, "w") as sortedHandle:
        geneDict = {}
        for chromosome, records in sorted(exonsDict.items(),key = lambda ch: ch[0]):
            exonLocation = lambda r: (int(r.id.split('|')[sortpara]),(int(r.id.split('|')[5]) - int(r.id.split('|')[4]) + 1)*(-1)) #exon start information
            lastStart = 0
            lastEnd = 0
            for record in sorted(records, key = exonLocation): #sorted by exon start
                headList = record.id.split("|")
                gene = headList[1]
                exonStart = headList[4]
                exonEnd = headList[5]
                if sortpara== 5:
                    if lastEnd == exonEnd:
                        continue
                    else:
                        lastStart = exonStart
                        lastEnd = exonEnd
                    record.description = ""
                    SeqIO.write(record,sortedHandle,"fasta")
                elif sortpara==4:
                    if lastStart == exonStart:
                        continue
                    else:
                        lastStart = exonStart
                        lastEnd = exonEnd
                    if gene not in geneDict:
                        geneDict.update({gene:1})
                    else:
                        geneDict[gene] += 1
                    if gene:
                        headList[1] = gene+'_'+str(geneDict[gene]) #give exons names
                        record.id = "|".join(headList)
                    record.description = ""
                    SeqIO.write(record,sortedHandle,"fasta")
                else:
                    raise Exception("Unknown sortpara: "+sortpara)
 

if __name__ == "__main__":
    exonsFasta = "human_exon.fas"
    sortedFasta = "human_exon_sorted_1.fas"
    sort(exonsFasta, sortedFasta)
    
