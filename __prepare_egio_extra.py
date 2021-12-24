import argparse
import os
import numpy as np
### extract information of protein-coding transcript from cDNA, CDS, GTF
### extract exon sequence in CDS region
### two files will be generated: species_exon_seqinfo.tab and species_traninfo.tab

def dropduplicated(listdat):
    if len(listdat) == 0:
        unilist = []
    if len(listdat) > 0:
        unilist = [listdat[0]]
        for i in range(0,len(listdat)):
            if not(listdat[i] in unilist):
                unilist.append(listdat[i])
    return unilist

def takeorf(elem):
    return [elem[2],elem[3]]

#gen=ug;nuseq=cdnaDt[tranidtmp];orfpos=orfDt[tranidtmp];cdscom=hstcom.iloc[i]["Exoncom"];chrtmp=hstcom.iloc[i]["Chrom"];strtmp=strtmp            
def extract_CDS(gen,nuseq,orfpos,cdscom,chrtmp,strtmp):
    orfst, orfen = orfpos
    if orfst < 0:
        orfst = 1
    orfseq = nuseq[(orfst-1):orfen]
    cds_pos = cdscom.split("|")
    cds_pos_int = [0 for i in range(0,len(cds_pos))]
    cdslen = [0 for i in range(0,len(cds_pos))]
    cdssum = [0 for i in range(0,len(cds_pos))]
    for i in range(0,len(cds_pos)):
        cds_pos_int[i] = [int(cds_pos[i].split(":")[0]), int(cds_pos[i].split(":")[1])]
        cdslen[i] = cds_pos_int[i][1] - cds_pos_int[i][0] + 1
        cdssum[i] = sum(cdslen)
    ##==============
    arr = np.array(cdssum)
    stidx = min(np.where(arr>=orfst)[0].tolist())
    enidx = min(np.where(arr>=orfen)[0].tolist())
    ##==============
    deltst = cdssum[stidx] - orfst
    delten = cdssum[enidx] - orfen
    ##--------------
    if strtmp == "+":
        cds_pos_int[stidx][0] = cds_pos_int[stidx][1] - deltst
        cds_pos_int[enidx][1] = cds_pos_int[enidx][1] - delten
    if strtmp == "-":
        cds_pos_int[stidx][1] = cds_pos_int[stidx][0] + deltst
        cds_pos_int[enidx][0] = cds_pos_int[enidx][0] + delten
    
    CDS = cds_pos_int[stidx:(enidx+1)]
    ##==============
    seqpos = []
    seqlen = [] 
    st = 0
    for i in range(0,len(CDS)):
        en = st + CDS[i][1] - CDS[i][0]
        seqpos.append([st,en])
        seqlen.append(en-st+1)
        st = en + 1

    if sum(seqlen) != (orfen-orfst+1):
        raise SyntaxError("extract_CDS: nucleotide loss or gain during sequence split")

    CDSm = []
    for i in range(0,len(CDS)):
        seqtmp = orfseq[seqpos[i][0]:(seqpos[i][1]+1)]
        if (CDS[i][1] - CDS[i][0] + 1) != len(seqtmp):
            raise SyntaxError("extract_CDS: cds sequence does not map to cds size")
        CDSm.append([gen,chrtmp,str(CDS[i][0]),str(CDS[i][1]),strtmp,seqtmp])
    
    return CDSm

def formatCDS(cdsnum, species):
    compen = "".join(["0" for i in range(0,10-len(str(cdsnum)))])
    cdsname = species + compen + str(cdsnum)
    return cdsname

def prepare_geio_extra(gtf, cdna, cds, species):
    outpath = os.getcwd() + "/extrainfo"
    ##===============================================
    species = str(species)

    print("extract orf position based on cDNA and CDS sequence")
    
    tran_name = []
    tran_cdna = []
    tf = open(str(cdna))
    line = tf.readline()
    transnum = 0
    seq = ''
    while line:
        if not line or line.startswith('#'):
            line = tf.readline()
            continue
        else:
            if ">" in line:
                if transnum > 0:
                    tran_name.append(tranid)
                    tran_cdna.append(seq)

                transnum = transnum + 1
                
                attr = line.split(" ")
                tranidtmp = [x for x in attr if ">" in x]

                tranid = tranidtmp[0].split(">")[1].split(".")[0].replace("\"","").replace("\n","")
                seq = ''
            else:
                seq = seq + line.replace("\n","")

        line = tf.readline()
    
    if transnum > 0:
        tran_name.append(tranid)
        tran_cdna.append(seq)
    
    tf.close()
    
    cdnaDt = {}
    for k, v in zip(tran_name,tran_cdna):
        cdnaDt[k] = v

    print(str(len(cdnaDt))+' of mRNA sequence is read')
    ##===============================================

    tran_name_cds = []
    tran_cds = []
    tc = open(str(cds))
    line = tc.readline()
    transnum = 0
    seq = ''
    while line:
        if not line or line.startswith('#'):
            line = tc.readline()
            continue
        else:
            if ">" in line:
                if transnum > 0:
                    tran_name_cds.append(tranid)
                    tran_cds.append(seq)

                transnum = transnum + 1
                
                attr = line.split(" ")
                tranidtmp = [x for x in attr if ">" in x]
                
                tranid = tranidtmp[0].split(">")[1].split(".")[0].replace("\"","").replace("\n","")

                record = 1
                seq = ''
            else:
                seq = seq + line.replace("\n","")
        
        line = tc.readline()
    
    if transnum > 0:
        tran_name_cds.append(tranid)
        tran_cds.append(seq)
    
    tc.close()

    cdsDt = {}
    for k, v in zip(tran_name_cds,tran_cds):
        cdsDt[k] = v

    print(str(len(cdsDt))+' CDS sequence is read')
    ##===============================================
    
    orfpos = []
    for i in range(0,len(tran_name_cds)):
        orfforstart = cdnaDt[tran_name_cds[i]].find(cdsDt[tran_name_cds[i]][0:(len(cdsDt[tran_name_cds[i]])//3)])
        orfrevstart = cdsDt[tran_name_cds[i]].find(cdnaDt[tran_name_cds[i]][0:(len(cdsDt[tran_name_cds[i]])//3)])
        orfstart = orfforstart+1 if orfforstart >= 0 else orfrevstart*-1
        orfend = orfstart + len(cdsDt[tran_name_cds[i]])-1 if orfforstart >= 0 else orfrevstart*-1+len(cdsDt[tran_name_cds[i]])
        orfpos.append([orfstart,orfend])

    orfDt = {}
    for k, v in zip(tran_name_cds,orfpos):
        orfDt[k] = v    
    
    ##===============================================
    ## ==== extract transcript information ====
    print("extract transcript information")
    tg = open(str(gtf))
    line = tg.readline()

    tranheader = ["EnsemblG", "GeneName", "GeneType", "EnsemblT", "TranType", "Orf", "Chrom", "Strand", "Exoncom"]
    hst = []
    hst.append(tranheader)
    
    starttmp = 0
    chrtmp = "0"
    strandname = "0"

    genetmp = ""
    isonamels = []
    isotypels = []
    isoexcols = []

    unigen = []
    geneisols = []
    genestrls = []
    genechrls = []
    while line:
        if not line or line.startswith('#'):
            line = tg.readline()
            continue
        else:
            try:
                chrom, source, feature, left, right, score, strand, frame, values = line.split('\t')
            except ValueError:
                continue
            if feature == "gene":
                genetest = 1
                
                if starttmp != 0:
                    geneisols.append(isolstmp)

                isolstmp = []

                if starttmp != 0 and (strandname == "+" or strandname == "-"):
                    if strandname == "-":
                        transpos.sort(reverse = True)
                    if strandname == "+":
                        transpos.sort(reverse = False)
                    transposstr = []
                    for i in range(0,len(transpos)):
                        transposstr.append(str(transpos[i][0]) + ":" + str(transpos[i][1]))
                    trantmp = [gid, gnm, gty, tid, typ, orf, chromosome, strandname, "|".join(transposstr)]
                    hst.append(trantmp)
                    
                    isonamels.append(tid)
                    isotypels.append(typ)
                    isoexcols.append("|".join(transposstr))

                attr = values.split("; ")
                gidtmp = [x for x in attr if "gene_id" in x]
                gnmtmp = [x for x in attr if "gene_name" in x]
                gtytmp = [x for x in attr if "gene_biotype" in x]
                
                gid = gidtmp[0].split(" ")[1].replace("\"","").replace("\n","")
                ############################
                unigen.append(gid)
                genestrls.append(strand)
                genechrls.append(chrom)
                ############################
                gty = gtytmp[0].split(" ")[1].replace("\"","").replace("\n","")
                if not(gnmtmp):
                    gnm = "."
                else:
                    gnm = gnmtmp[0].split(" ")[1].replace("\"","").replace("\n","")
                
            if feature == "transcript":
                
                if starttmp != 0 and (strandname == "+" or strandname == "-") and genetest != 1:
                    if strandname == "-":
                        transpos.sort(reverse = True)
                    if strandname == "+":
                        transpos.sort(reverse = False)
                    transposstr = []
                    for i in range(0,len(transpos)):
                        transposstr.append(str(transpos[i][0]) + ":" + str(transpos[i][1]))
                    trantmp = [gid, gnm, gty, tid, typ, orf, chromosome, strandname, "|".join(transposstr)]
                    hst.append(trantmp)

                    isonamels.append(tid)
                    isotypels.append(typ)
                    isoexcols.append("|".join(transposstr))

                genetest = 0
                attr = values.split("; ")
                orftmp = [x for x in attr if "orf" in x]
                tidtmp = [x for x in attr if "transcript_id" in x]
                typtmp = [x for x in attr if "transcript_biotype" in x]
                
                
                tid = tidtmp[0].split(" ")[1].replace("\"","").replace("\n","")
                isolstmp.append(tid)

                if tid in orfDt:
                    orf = "|".join([str(orfDt[tid][0]),str(orfDt[tid][1])])
                else:
                    orf = ".|."

                typ = typtmp[0].split(" ")[1].replace("\"","").replace("\n","")
                
                strandname = strand
                chromosome = chrom
                transpos = []
                starttmp = starttmp + 1

            if feature == "exon":                
                transpos.append([int(left), int(right)])
                
            line = tg.readline()

            if chrtmp != chrom:
                chrtmp = chrom
                
    if strandname == "-":
        transpos.sort(reverse = True)
    if strandname == "+":
        transpos.sort(reverse = False)
    transposstr = []
    for i in range(0,len(transpos)):
        transposstr.append(str(transpos[i][0]) + ":" + str(transpos[i][1]))    
    trantmp = [gid, gnm, gty, tid, typ, orf, chromosome, strandname, "|".join(transposstr)]
    hst.append(trantmp)

    isonamels.append(tid)
    isotypels.append(typ)
    isoexcols.append("|".join(transposstr))

    geneisols.append(isolstmp)

    tg.close()

    ##===============================================    !!!!!!!!!!!!!!!!!!
    typeDt = {}
    for k, v in zip(isonamels,isotypels):
        typeDt[k] = v
    
    excoDt = {}
    for k, v in zip(isonamels,isoexcols):
        excoDt[k] = v

    geneisoDt = {}
    for k, v in zip(unigen,geneisols):
        geneisoDt[k] = v

    genestrDt = {}
    for k, v in zip(unigen,genestrls):
        genestrDt[k] = v

    genechrDt = {}
    for k, v in zip(unigen,genechrls):
        genechrDt[k] = v
    ##===============================================    !!!!!!!!!!!!!!!!!!
    #print("write transcript information:")

    outtsqlpath = outpath + "/" + species + ".tran"
    ff = open(outtsqlpath,'w')
    for i in hst:
        ff.write("\t".join(i)+"\n")
    ff.close()    
    
    ##===============================================  
    ## ==== get final transcriptome fasta based on gtf ====

    testgennum = [ int((i)*0.1*len(unigen)) for i in range(1,11)]

    print("extract coding exon sequence from "+ str(len(unigen)) + " genes")

    cdsheader = ["EnsemblG","ID","Chrom","Start","End","Strand","Seq"]
    CDSinfostore = []
    CDSinfostore.append(cdsheader)
    countg = 0
    counte = 0
    counttest = 0
    for ug in unigen:
        countg = countg + 1
        ugiso = geneisoDt[ug]
        
        uftmp = []
        for iso in ugiso:
            if typeDt[iso] == "protein_coding" or typeDt[iso] == "nonsense_mediated_decay":
                if iso in orfDt:
                    uftmp = uftmp + extract_CDS(ug,cdnaDt[iso],orfDt[iso],excoDt[iso],genechrDt[ug],genestrDt[ug])
        
        CDSinfotmp = dropduplicated(uftmp)

        if genestrDt[ug] == "+":
            CDSinfotmp.sort(key = takeorf, reverse = False)
        else:
            CDSinfotmp.sort(key = takeorf, reverse = True)

        for ii in range(0,len(CDSinfotmp)):
            counte = counte + 1
            CDSinfotmp[ii].insert(1,formatCDS(counte,species))

        CDSinfostore = CDSinfostore + CDSinfotmp
        
        if countg in testgennum:
            counttest = counttest + 1
            print(".........." + str(counttest*10) + "%")
    ##===============================================

    outpathCDS = outpath + "/" + species + ".exon"
    ff = open(outpathCDS,'w')
    for i in CDSinfostore:
        ff.write("\t".join(i)+"\n")
    ff.close()    
    
    ##===============================================
    print("\ntransform CDS-Exon into fasta")
    nchrrow = 60
    cdsfasta = []
    for i in range(1,len(CDSinfostore)):
        orfseq = CDSinfostore[i][6]
        transname = ">" + CDSinfostore[i][1]
        cdsfasta.append(transname + "\n")

        nrowtmp = len(orfseq) // nchrrow
        comp = len(orfseq) % nchrrow
        if comp == 0:
            nrow = nrowtmp
        else:
            nrow = nrowtmp + 1
        for i in range(0,nrow):
            if i == nrow:
                cdsfasta.append(orfseq[(i*nchrrow):len(orfseq)] + "\n")
            if i < nrow:
                cdsfasta.append(orfseq[(i*nchrrow):((i+1)*nchrrow)] + "\n")

    outtfapath = outpath + "/" + species + ".exonfasta"
    ff = open(outtfapath,'w')
    for i in cdsfasta:
        ff.write(i)
    ff.close()

## main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract information of protein coding transcipt from cDNA, CDS fasta and GTF files')
    parser.add_argument('--gtf',
        nargs='?',
        help='transcriptome gtf file')
    parser.add_argument('--cdna',
        nargs='?',
        help='cdna fasta file')
    parser.add_argument('--cds',
        nargs='?',
        help='cds fasta file')
    parser.add_argument('--species',
        nargs='?',
        help='species name')

    args = parser.parse_args()
    prepare_geio_extra(args.gtf, args.cdna, args.cds, args.species)

#cdna = "/Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/cdna.fa"
#cds =  "/Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/cds.fa"
#gtf= "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/#0mmu/Mus_musculus.GRCm38.102.gtf"
#outpath =  "/Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/"

#cd /Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/
#python prepare_egio_extra.py --gtf Mus_musculus.GRCm38.102.gtf --cdna cdna.fa --cds cds.fa --species mmu

#cd /Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/
#python /Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/prepare_egio_extra.py --gtf /Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/#0mml/Macaca_mulatta.Mmul_10.102.gtf --cdna /Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/#0mml/Macaca_mulatta.Mmul_10.cdna.all.fa --cds /Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/#0mml/Macaca_mulatta.Mmul_10.cds.all.fa --species mml

#python prepare_egio_extra.py --gtf Mus_musculus.GRCm38.102.gtf --cdna cdna.fa --cds cds.fa --species mmu
