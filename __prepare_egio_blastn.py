import argparse
import os
import pandas as pd

def takenamescore(elem):
    return [elem[0],elem[4]]

def summaryblastn(species1, blast1, exonanno1, species2, blast2, exonanno2, orthogen, coverage):
    ##===============================================
    print("\n\n\norganize orthologous gene pairs")
    homogen = pd.read_table(str(orthogen),header=0,sep='\t')
    homogpair = ["." for i in range(0,len(homogen))]
    for i in range(0,len(homogen)):
        homogpair[i] = ":".join([homogen.iloc[i][str(species1)],homogen.iloc[i][str(species2)]])
    
    ##==============================================================
    print("organize cds region information")

    hsacds = []
    tf = open(str(exonanno1))
    line = tf.readline()   ## header
    line = tf.readline()
    while line:
        hsacds.append(line.split("\t"))
        line = tf.readline()
    tf.close()

    hsacdsb = ["." for i in range(0,len(hsacds))]
    hsacdsg = ["." for i in range(0,len(hsacds))]
    hsacdsl = [0 for i in range(0,len(hsacds))]
    for i in range(0,len(hsacds)):
        hsacdsb[i] = hsacds[i][1]
        hsacdsg[i] = hsacds[i][0]
        hsacdsl[i] = int(hsacds[i][4])-int(hsacds[i][3])+1
    
    hsacdsgDt = {}
    for k, v in zip(hsacdsb,hsacdsg):
        hsacdsgDt[k] = v
    hsacdslDt = {}
    for k, v in zip(hsacdsb,hsacdsl):
        hsacdslDt[k] = v
    ##----------------------------------------------------------------
    nhsacds = []
    tf = open(str(exonanno2))
    line = tf.readline()   ## header
    line = tf.readline()
    while line:
        nhsacds.append(line.split("\t"))
        line = tf.readline()
    tf.close()

    nhsacdsb = ["." for i in range(0,len(nhsacds))]
    nhsacdsg = ["." for i in range(0,len(nhsacds))]
    nhsacdsl = [0 for i in range(0,len(nhsacds))]
    for i in range(0,len(nhsacds)):
        nhsacdsb[i] = nhsacds[i][1]
        nhsacdsg[i] = nhsacds[i][0]
        nhsacdsl[i] = int(nhsacds[i][4])-int(nhsacds[i][3])+1
    
    nhsacdsgDt = {}
    for k, v in zip(nhsacdsb,nhsacdsg):
        nhsacdsgDt[k] = v
    nhsacdslDt = {}
    for k, v in zip(nhsacdsb,nhsacdsl):
        nhsacdslDt[k] = v

    ##==============================================================
    print("filter BLASTN results within orthologous gene pairs")
    print("analysis species1")

    hsaquery = []
    tf = open(str(blast1))
    line = tf.readline()
    while line:
        hsaquery.append(line.split("\t"))
        line = tf.readline()
    tf.close()
    
    #hsaquery.columns = ["query","subject","iden","alilen","mismlen","gap","queryst","queryen","subst","suben","eval","score"]

    hsaquery2list = ["." for _ in range(0,len(hsaquery))]
    countlsr = 0
    for i in range(0,len(hsaquery)):
        if ":".join([hsacdsgDt[hsaquery[i][0]],nhsacdsgDt[hsaquery[i][1]]]) in homogpair:
            hsaquery2list[countlsr] = "|".join([hsaquery[i][0],hsaquery[i][1],hsaquery[i][2],hsaquery[i][3],hsaquery[i][11]])
            countlsr = countlsr + 1
            
        if i % 100000 == 0 and i > 0:
            print(str(i) + " records are analyzed")
    print(str(i) + " records are analyzed")

    hsaquery2list = hsaquery2list[0:countlsr]
    ##----------------------------------------------------------------
    print("analysis species2")

    nhsaquery = []
    tf = open(str(blast2))
    line = tf.readline()
    while line:
        nhsaquery.append(line.split("\t"))
        line = tf.readline()
    tf.close()

    #nhsaquery.columns = ["query","subject","iden","alilen","mismlen","gap","queryst","queryen","subst","suben","eval","score"]
    
    nhsaquery2list = ["." for _ in range(0,len(nhsaquery))]
    countlsr = 0
    for i in range(0,len(nhsaquery)):
        if ":".join([hsacdsgDt[nhsaquery[i][1]],nhsacdsgDt[nhsaquery[i][0]]]) in homogpair:
            nhsaquery2list[countlsr] = "|".join([nhsaquery[i][1],nhsaquery[i][0],nhsaquery[i][2],nhsaquery[i][3],nhsaquery[i][11]])
            countlsr = countlsr + 1

        if i % 100000 == 0 and i > 0:
            print(str(i) + " records are analyzed")
    print(str(i) + " records are analyzed")

    nhsaquery2list = nhsaquery2list[0:countlsr]

    ##----------------------------------------------------------------
    common = list(set(hsaquery2list).intersection(set(nhsaquery2list)))
    ##==============================================================
    
    print("filter BLASTN results with reciprocal coverage over threshold")
    
    maplist = [[] for _ in range(0,len(common))]
    for i in range(0,len(common)):
        tmp = common[i].split("|")
        maplist[i] = [tmp[0],tmp[1],float(tmp[2]),float(tmp[3]),float(tmp[4])]

    maplist2 = [[] for _ in range(0,len(maplist))]
    countlsr = 0
    for i in range(0,len(maplist)):
        hsacover = maplist[i][2]*maplist[i][3]/hsacdslDt[maplist[i][0]]/100
        nhsacover = maplist[i][2]*maplist[i][3]/nhsacdslDt[maplist[i][1]]/100
        if hsacover >= coverage and nhsacover >= coverage:
            maplist2[countlsr] = maplist[i]
            countlsr = countlsr + 1

    maplist = maplist2[0:countlsr]
    maplist.sort(key = takenamescore, reverse = True)
    ##==============================================================

    print("filter BLASTN results with reciprocal best hit\n\n\n")
    
    store = [[] for _ in range(0,len(maplist))]
    countlsr = 0
    store[0] = maplist[countlsr]
    countlsr = countlsr + 1
    for i in range(1,len(maplist)):
        if maplist[i][0] != maplist[i-1][0]:
            store[countlsr] = maplist[i]
            countlsr = countlsr + 1
        else:
            if maplist[i][4] == maplist[i-1][4]:
                store[countlsr] = maplist[i]
                countlsr = countlsr + 1

    tarmapping = store[0:countlsr]
    tarmapping.sort(key = takenamescore)

    ##==============================================================

    outtfapath = os.getcwd() + "/extrainfo/blastn_mapping.tab"
    ff = open(outtfapath,'w')
    ff.write("\t".join([str(species1),str(species2)+"\n"]))
    for i in tarmapping:
        ff.write("\t".join([i[0],i[1]+"\n"]))
    ff.close()
    ##****************************************************************
    ##****************************************************************
    ##****************************************************************

## main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract exon mappings from reciprocal BLASTN')
    parser.add_argument('--species1',
        nargs='?',
        help='name of species1 in orthologous gene annotation table')
    parser.add_argument('--blast1',
        nargs='?',
        help='balstn file using species1 as query')
    parser.add_argument('--exonanno1',
        nargs='?',
        help='exon annotation file of species1')
    parser.add_argument('--species2',
        nargs='?',
        help='name of species2 in orthologous gene annotation table')
    parser.add_argument('--blast2',
        nargs='?',
        help='balstn file using species2 as query')
    parser.add_argument('--exonanno2',
        nargs='?',
        help='exon annotation file of species2')
    parser.add_argument('--orthogen',
        nargs='?',
        help='orthologous gene annotation of species1 and species2: please follow the species order in the file: species1, then species2')       
    parser.add_argument('--coverage',
        nargs='?',
        type=float,
        default=0.8,
        help='reciprocal coverage')

    args = parser.parse_args()
    summaryblastn(args.species1, args.blast1, args.exonanno1, args.species2, args.blast2, args.exonanno2, args.orthogen, args.coverage)


#blast1 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/2blastn/ptr_hsa.blast"
#exonanno1 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0hsa/mysql_exon_seqinfo.tab"

#blast2 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/2blastn/hsa_ptr.blast"
#exonanno2 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0ptr/mysql_exon_seqinfo.tab"

#orthogen = "/Users/jeffma/MJFPhD/HsSpDB/step6_homologs/2homologs_gene_blastp_inparanoid8.dat"

#coverage = 0.8
#outpath =  "/Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/"
#summaryblastn(blast1, exonanno1, blast2, exonanno2, orthogen, coverage, outpath)

#cd /Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/
#blast1 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/2blastn/ptr_hsa.blast"
#exonanno1 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0hsa/mysql_exon_seqinfo.tab"

#blast2 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/2blastn/hsa_ptr.blast"
#exonanno2 = "/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0ptr/mysql_exon_seqinfo.tab"

#orthogen = "/Users/jeffma/MJFPhD/HsSpDB/step6_homologs/2homologs_gene_blastp_inparanoid8.dat"


#outpath =  "/Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/"
#cd /Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/
#python __prepare_egio_blastn.py --species1 hsa --blast1 hsa-ptr --exonanno1 hsa.exon --species2 ptr --blast2 ptr-hsa --exonanno2 ptr.exon --orthogen homolog_gene.dat --coverage 0.8