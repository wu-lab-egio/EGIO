import argparse
import string
import ctypes
import os
import copy
import math
import pandas as pd
from multiprocessing import Pool
##======================================================##
##                      basic function                  ##
##======================================================##
##======================================================##

align= ctypes.cdll.LoadLibrary(os.getcwd() + '/___pairwisealign.so')
align.pairwisealign.argtypes = [ctypes.c_char_p]
align.pairwisealign.restype = ctypes.c_char_p

def which(judge):
    backdata = []
    for i in range(0,len(judge)):
        if judge[i]:
            backdata.append(i)
    return backdata

def listwhich(listdata,deal,value):
    judge = []
    for i in list(range(0,len(listdata))):
        if deal == "==":
            if listdata[i] == value:
                judge.append(True)
            else:
                judge.append(False)
        if deal == "!=":
            if listdata[i] == value:
                judge.append(False)
            else:
                judge.append(True)
        if deal == ">":
            if listdata[i] <= value:
                judge.append(False)
            else:
                judge.append(True)
        if deal == ">=":
            if listdata[i] < value:
                judge.append(False)
            else:
                judge.append(True)
        if deal == "<":
            if listdata[i] >= value:
                judge.append(False)
            else:
                judge.append(True)
        if deal == "<=":
            if listdata[i] > value:
                judge.append(False)
            else:
                judge.append(True)           
    jindex = which(judge)
    return(jindex)

def intersect(ls1, ls2):
    store = []
    for i in range(0,len(ls1)):
        if ls1[i] in ls2:
            store.append(ls1[i])
    return store

def union(ls1, ls2):
    store = []
    if not(isinstance(ls1,list)) and not(isinstance(ls2,list)):
        store.append(ls1)
        if not(ls2 in store):
            store.append(ls2)
    elif isinstance(ls1,list) and not(isinstance(ls2,list)):
        store.append(ls2)
        if len(ls1) > 0:
            for j in range(0,len(ls1)):
                if ls1[j] in store:
                    continue
                else:
                    store.append(ls1[j])
    elif isinstance(ls2,list) and not(isinstance(ls1,list)):
        store.append(ls1)
        if len(ls2) > 0:
            for j in range(0,len(ls2)):
                if ls2[j] in store:
                    continue
                else:
                    store.append(ls2[j])
    elif isinstance(ls2,list) and isinstance(ls1,list):
        if len(ls1) == 0 and len(ls2) > 0:
            for j in range(0,len(ls2)):
                if ls2[j] in store:
                    continue
                else:
                    store.append(ls2[j])
        elif len(ls1) > 0 and len(ls2) == 0:
            for j in range(0,len(ls1)):
                if ls1[j] in store:
                    continue
                else:
                    store.append(ls1[j])
        elif len(ls1) > 0 and len(ls2) > 0:
            for j in range(0,len(ls1)):
                if ls1[j] in store:
                    continue
                else:
                    store.append(ls1[j])
            for j in range(0,len(ls2)):
                if ls2[j] in store:
                    continue
                else:
                    store.append(ls2[j])
    store.sort()
    return store

def typeinto(listd,colname,astype=int):
    listdat = listd["dat"]
    listname = listd["coln"]
    coli = listwhich(listname,"==",colname)[0]
    for i in range(0,len(listdat)):
        listdat[i][coli] = astype(listdat[i][coli])
    info = {"dat":listdat,"coln":listname}
    return info

def listtypeinto(listdat,astype=int):
    newst = []
    for i in range(0,len(listdat)):
        newst.append(astype(listdat[i]))
    return newst

def unique(listdat,order=True):
    store = []
    if len(listdat) > 0:
        for i in range(0,len(listdat)):
            if not(listdat[i] in store):
                store.append(listdat[i])
        if order:
            store.sort(reverse = False)
    return store

def selectcol(listdatd,target):
    store =[]
    listdat = listdatd["dat"]
    colname = listdatd["coln"]
    seli = listwhich(colname,"==",target)[0]
    for i in range(0,len(listdat)):
        store.append(listdat[i][seli])
    return store

def listasint(listdat,deal="min"):
    ## only be used in homo-colinearity test
    ## transform list in list to int
    if not(isinstance(listdat,list)):
        store = [listdat]
    else:
        store = []
        for i in range(0,len(listdat)):
            if isinstance(listdat[i],list):
                if deal == "min":
                    store.append(min(listdat[i]))
                if deal == "max":
                    store.append(max(listdat[i]))
                if deal == "all":
                    for j in range(0,len(listdat[i])):
                        store.append(listdat[i][j])
            else:
                store.append(listdat[i])
    return store

def selectele(listdat,targeti):
    store = []
    if len(targeti) > 0:
        for i in range(0,len(targeti)):
            store.append(listdat[targeti[i]])
    return store

def selectrow(listdatd,targeti):
    store = []
    if len(targeti) > 0:
        for i in range(0,len(targeti)):
            store.append(listdatd["dat"][targeti[i]])
    return {"dat":store,"coln":listdatd["coln"]}

def selectlsls(listd,rowi,coln="all"):
    listdat = listd["dat"]
    if coln == "all":
        coln = listd["coln"]
    coli = []
    for i in range(0,len(coln)):
        coli.append(listwhich(listd["coln"],"==",coln[i])[0])
    new_l = []
    for i in range(0,len(rowi)):
        new_r = []
        for j in range(0,len(coli)):
            new_r.append(listdat[rowi[i]][coli[j]])
        new_l.append(new_r)
    listsum = {"dat":new_l,"coln":coln}
    return listsum

def listcal(listdat,calcu,value):
    if len(listdat) == 0:
        listtmp = []
    else:
        listtmp = []
        for i in range(0,len(listdat)):
            if calcu == "+":
                listtmp.append(listdat[i] + value)
            if calcu == "-":
                listtmp.append(listdat[i] - value)
            if calcu == "*":
                listtmp.append(listdat[i] * value)
            if calcu == "/":
                listtmp.append(listdat[i] / value)
            if calcu == "%":
                listtmp.append(listdat[i] % value)
            if calcu == "^":
                listtmp.append(listdat[i] ** value)
    return listtmp

def listlistcal(list1,calcu,list2):
    if len(list1) != len(list2):
        raise SyntaxError("listlistcal: the two lists do not have the same length")
    listtmp = [0 for _ in range(0,len(list1))]
    for i in range(0,len(list1)):
        if calcu == "+":
            listtmp[i] = list1[i] + list2[i]
        if calcu == "-":
            listtmp[i] = list1[i] - list2[i]
        if calcu == "*":
            listtmp[i] = list1[i] * list2[i]
        if calcu == "/":
            listtmp[i] = list1[i] / list2[i]
        if calcu == "%":
            listtmp[i] = list1[i] % list2[i]
        if calcu == "^":
            listtmp[i] = list1[i] ** list2[i]
    return listtmp

def listdrop(listdat,dropi):
    if len(dropi) == 0:
        store = listdat
    else:
        store = []
        for i in range(0,len(listdat)):
            if not(i in dropi):
                store.append(listdat[i])
    return store

def listddrop(listdatd,dropi):
    store = []
    for i in range(0,len(listdatd["dat"])):
        if not(i in dropi):
            store.append(listdatd["dat"][i])
    return {"dat":store,"coln":listdatd["coln"]}

def listdinsert(listdatd,insi,element):
    if len(listdatd["dat"]) == 0:
        listdatd["dat"].insert(insi,element)
    else:
        if len(element) != len(listdatd["dat"][0]):
            raise SyntaxError("listdinsert: insert element do not map to the data dimension (column)")
        listdatd["dat"].insert(insi,element)
    return {"dat":listdatd["dat"],"coln":listdatd["coln"]}

def transpose(listdat):
    new_list = []
    for i in range(0,len(listdat[0])):
        listdat1 = []
        for j in range(0,len(listdat)):
            listdat1.append(listdat[j][i])
        new_list.append(listdat1)
    return new_list

def listsort(listd,by,increase):
    if len(by) > 2:
        raise SyntaxError("listsort: the max element to sort is 2")
    listdat = listd["dat"]
    if len(by) == 1:
        by1 = listwhich(listd["coln"],"==",by[0])[0]
        def takeone(elem):
            return elem[by1]
        if increase:
            listdat.sort(key =takeone, reverse = False)
        else:
            listdat.sort(key =takeone, reverse = True)
    
    else:
        by1 = listwhich(listd["coln"],"==",by[0])[0]
        by2 = listwhich(listd["coln"],"==",by[1])[0]
        def taketwo(elem):
            return elem[by1], elem[by2]
        if increase:
            listdat.sort(key =taketwo, reverse = False)
        else:
            listdat.sort(key =taketwo, reverse = True)    
    listsum = {"dat":listdat,"coln":listd["coln"]}
    return listsum

def dropduplicated(listdat):
    if len(listdat) == 0:
        unilist = []
    if len(listdat) > 0:
        unilist = [listdat[0]]
        for i in range(0,len(listdat)):
            if not(listdat[i] in unilist):
                unilist.append(listdat[i])
    return unilist

def setdiff(list1,list2):
    difflist=[]
    for i in range(0,len(list1)):
        if not(list1[i] in list2):
            difflist.append(list1[i])
    return difflist

def transform(data):
    store = []
    for i in range(0, len(data)):
        rowstore = []
        for j in range(0, len(data.iloc[0])):
            rowstore .append(str(data.iloc[i][j]))
        store.append("#".join(rowstore))
    ffstr = "@".join(store)
    return ffstr

def aslsls(datastr,datcolumns):
    ## transform string as dataframe
    datastr = datastr.replace("^",";") # replace "^" to ";", which ";" is not allowed during php calling py
    rowdata = datastr.split("@")       # seperate row information, which is seperated by "@"
    rowlen = len(rowdata)              # seperate col information, which is seperated by "#"
    datainfo = []                                 # set a list to store seperated information
    for i in list(range(0,rowlen)):               # traversal to seperate each row to seperated columns
        datainfo.append(rowdata[i].split("#"))    # seperate row to seperated columns
    if len(datainfo[0]) != len(datcolumns):
        raise SyntaxError("aslsls: the length of colname names is not match to the data dimension")
    datsum = {"dat":datainfo,"coln":datcolumns}
    return datsum

##======================================================##
##                  translate from cDNA                 ##
##======================================================##

def translate(cdna):
    ## translate from codon to AA
    nt_to_aa = {'AAA':'K','AAT':'N','AAG':'K','AAC':'N','ATA':'I','ATT':'I','ATG':'M','ATC':'I','AGA':'R','AGT':'S','AGG':'R','AGC':'S','ACA':'T','ACT':'T','ACG':'T','ACC':'T','TAA':'','TAT':'Y','TAG':'','TAC':'Y','TTA':'L','TTT':'F','TTG':'L','TTC':'F','TGA':'','TGT':'C','TGG':'W','TGC':'C','TCA':'S','TCT':'S','TCG':'S','TCC':'S','GAA':'E','GAT':'D','GAG':'E','GAC':'D','GTA':'V','GTT':'V','GTG':'V','GTC':'V','GGA':'G','GGT':'G','GGG':'G','GGC':'G','GCA':'A','GCT':'A','GCG':'A','GCC':'A','CAA':'Q','CAT':'H','CAG':'Q','CAC':'H','CTA':'L','CTT':'L','CTG':'L','CTC':'L','CGA':'R','CGT':'R','CGG':'R','CGC':'R','CCA':'P','CCT':'P','CCG':'P','CCC':'P'}
    cdnalen = len(cdna)
    proseq = ""
    for i in range(0,cdnalen // 3):
        try:
            aminoa = nt_to_aa[cdna[i*3:(i*3+3)]]
        except:
            aminoa = "X"
        proseq = proseq + aminoa
    return proseq

def extract_CDSMatrix_from_Exon(orfpos,exoncom,strtmp):
    ## extract CDS component from exon component
    orfst, orfen = orfpos.split("|")
    orfst, orfen = int(orfst), int(orfen)
    if orfst < 0:
        orfst = 1

    cds_pos = exoncom.split("|")
    cds_pos_int = [0 for i in range(0,len(cds_pos))]
    cdslen = [0 for i in range(0,len(cds_pos))]
    cdssum = [0 for i in range(0,len(cds_pos))]
    for i in range(0,len(cds_pos)):
        cds_pos_int[i] = [int(cds_pos[i].split(":")[0]), int(cds_pos[i].split(":")[1])]
        cdslen[i] = cds_pos_int[i][1] - cds_pos_int[i][0] +1
        cdssum[i] = sum(cdslen)
    ##==============
    stidx = min(listwhich(cdssum,">=",orfst))
    deltst = cdssum[stidx] - orfst
    ##--------------
    enidx = min(listwhich(cdssum,">=",orfen))
    delten = cdssum[enidx] - orfen
    ##--------------
    if strtmp == "+":
        cds_pos_int[stidx][0] = cds_pos_int[stidx][1] - deltst
        cds_pos_int[enidx][1] = cds_pos_int[enidx][1] - delten
    if strtmp == "-":
        cds_pos_int[stidx][1] = cds_pos_int[stidx][0] + deltst
        cds_pos_int[enidx][0] = cds_pos_int[enidx][0] + delten
    
    CDS = cds_pos_int[stidx:(enidx+1)]
    return CDS

def tranform_matrix_str(cdscom):
    ## transform matrix to string
    cdsstrtmp = ["." for _ in range(0,len(cdscom))]
    for i in range(0,len(cdscom)):
        cdsstrtmp[i] = str(cdscom[i][0])+":"+str(cdscom[i][1])
    return "|".join(cdsstrtmp)

def extract_pro_seq(cdscom,cdshscds,orfpos):
    ## translate CDS matrix into AA sequence
    orfst, orfen = orfpos.split("|")
    orfst, orfen = int(orfst), int(orfen)
    seq = ""
    if orfst < 0:
        seq = seq + "".join(["N" for _ in range(0,abs(orfst))])
    ##----------------------------    
    for i in range(0,len(cdscom)):
        idx = intersect(cdshscds[(cdshscds.Start==cdscom[i][0])].index.tolist(),cdshscds[(cdshscds.End==cdscom[i][1])].index.tolist())
        seq = seq + cdshscds.loc[idx[0]]["Seq"]
    ##----------------------------
    aaseq = translate(seq)
    return aaseq

def translatecdna(genhscdsdat,cdshscds,strandhsstr):
    ## translate exon com matrix to cds com matrix
    exonnew = []
    for i in range(0,len(genhscdsdat)):
        exoncom = genhscdsdat.iloc[i]["Exoncom"]
        orfpos = genhscdsdat.iloc[i]["Orf"]
        cdscom = extract_CDSMatrix_from_Exon(orfpos,exoncom,strandhsstr)
        proseq = extract_pro_seq(cdscom,cdshscds,orfpos)
        exonnew.append([genhscdsdat.iloc[i]["EnsemblT"],orfpos,tranform_matrix_str(cdscom),proseq])
    exonnew = pd.DataFrame(exonnew)
    exonnew.columns = ["EnsemblT","Orf","CDS","AAseq"]
    return exonnew

##======================================================##
##                   group CDS-exons ***                ##
##======================================================##

def mergeOverlap(datad,grpstr):
    ## organize exon group based on exon region
    if grpstr == 1:
        datad = listsort(datad,by=["start","end"],increase = True)
        sti = listwhich(datad["coln"],"==","start")[0]
        eni = listwhich(datad["coln"],"==","end")[0]
        sqi = listwhich(datad["coln"],"==","seq")[0]
        stff = datad["dat"][0][sti]
        enff = datad["dat"][0][eni]
        seqff = datad["dat"][0][sqi]
        for i in range(1,len(datad["dat"])):
            sttmp = datad["dat"][i][sti]
            entmp = datad["dat"][i][eni]
            seqtmp = datad["dat"][i][sqi]
            if enff >= sttmp and enff < entmp:
                seq_ahead_tmp = list(seqff)[0:(sttmp - stff)]
                seqff = "".join(seq_ahead_tmp) + seqtmp
                stff = min(stff,sttmp)
                enff = max(enff, entmp)
    if grpstr == -1:
        datad = listsort(datad,by=["end","start"],increase = False)
        ## strand = -1, the large position is small
        sti = listwhich(datad["coln"],"==","start")[0]
        eni = listwhich(datad["coln"],"==","end")[0]
        sqi = listwhich(datad["coln"],"==","seq")[0]
        stff = datad["dat"][0][sti]
        enff = datad["dat"][0][eni]
        seqff = datad["dat"][0][sqi]
        for i in range(1,len(datad["dat"])):
            sttmp = datad["dat"][i][sti]
            entmp = datad["dat"][i][eni]
            seqtmp = datad["dat"][i][sqi]
            if stff > sttmp and stff <= entmp:
                seq_ahead_tmp = list(seqff)[0:(enff-entmp)]
                seqff = "".join(seq_ahead_tmp) + seqtmp
                stff = min(stff,sttmp)
                enff = max(enff, entmp)

    return seqff

def regRuler(regidxd,genstr):
    ## estimate exon ruler
    regidxd = typeinto(regidxd,"start",int)
    regidxd = typeinto(regidxd,"end",int)
    regname = regidxd["coln"]
    if genstr == 1:
        regcomd = listsort(regidxd,by=["start","end"],increase = True)
    if genstr == -1:
        regcomd = listsort(regidxd,by=["end","start"],increase = False)
    ## get group of each region
    group = []
    gtmp = 1
    sti = listwhich(regname,"==","start")[0]
    eni = listwhich(regname,"==","end")[0]
    stref = regcomd["dat"][0][sti]
    enref = regcomd["dat"][0][eni]
    for i in range(0,len(regcomd["dat"])):
        sttst = regcomd["dat"][i][sti]
        entst = regcomd["dat"][i][eni]
        if entst < stref or sttst > enref:
            gtmp = gtmp + 1
            group.append(gtmp)
            stref = sttst
            enref = entst
        else:
            group.append(gtmp)
            stref = min(stref, sttst)
            enref = max(enref,entst)
    ## summarize each group
    gid = unique(group)
    gsum = []
    starttmp = selectcol(regcomd,"start")
    endtmp = selectcol(regcomd,"end")
    #regseq = selectcol(regcomd,"seq")
    for i in range(0,len(gid)):
        rgri = listwhich(group,"==",gid[i])
        gst = min(selectele(starttmp,rgri))
        gen = max(selectele(endtmp,rgri))
        if len(rgri) == 1:
            gsum.append([gid[i],gst,gen,"nonunited"])
        else:
            gsum.append([gid[i],gst,gen,"united"])
    gsumname = ["group","start","end","etype"]
    gsumd = {"dat":gsum,"coln":gsumname}
    ## calculate the aligned position
    glentmp = listlistcal(selectcol(gsumd,"end"),"-",selectcol(gsumd,"start"))
    glen = listcal(glentmp,"+",1)
    ginst = [[] for _ in range(0,len(glen))]
    ginst[0] = 1
    ginen = [[] for _ in range(0,len(glen))]
    ginen[0]= glen[0]
    for i in range(1,len(glen)):
        ginst[i] = ginen[i-1] + 1
        ginen[i] = ginen[i-1] + glen[i]
    posinfo=[]
    for i in range(0,len(glen)):
        posinfo.append([glen[i],ginst[i],ginen[i]])
    posinfod = {"dat":posinfo,"coln":["len","instart","inend"]}    
    ## get represent region id and sequence
    gsti = listwhich(gsumname,"==","start")[0]
    geni = listwhich(gsumname,"==","end")[0]
    regi = [[] for _ in list(range(0,len(gsum)))]
    seqtmp = [[] for _ in list(range(0,len(gsum)))]
    regseqi = listwhich(regname,"==","seq")[0]
    #-----------------------
    for i in range(0,len(gsum)):
        j1 = listwhich(starttmp, ">=", gsum[i][gsti])
        j2 = listwhich(endtmp, "<=", gsum[i][geni])
        j = intersect(j1,j2)
        regi[i] = str(gsum[i][gsti])+"-"+str(gsum[i][geni])
        if len(j) == 1:
            seqtmp[i] = regcomd["dat"][j[0]][regseqi]
        else:
            # look for new start, end of the first region
            grpseqtmp = []
            for ji in range(0,len(j)):
                grpseqtmp.append(regcomd["dat"][j[ji]])
            grpseqtmpd = {"dat":grpseqtmp, "coln":["id","start","end","seq"]}
            seqmerge = mergeOverlap(datad = grpseqtmpd, grpstr = genstr)
            seqtmp[i] = seqmerge
    #-----------------------
    addinfo = []
    for i in range(0,len(regi)):
        addinfo.append([regi[i],seqtmp[i]])
    addinfod = {"dat":addinfo,"coln":["id","seq"]}
    ## get final region ruler
    gsum_t = transpose(gsumd["dat"])
    posinfo_t = transpose(posinfod["dat"])
    addinfo_t = transpose(addinfod["dat"])
    #-----------------------
    ruler_t = gsum_t + posinfo_t + addinfo_t
    ruler = transpose(ruler_t)
    rulername = gsumd["coln"] + posinfod["coln"] + addinfod["coln"]
    rulerd = {"dat":ruler,"coln":rulername}
    return rulerd

##======================================================##
##   pairwise alignment and homo colinearity test ***   ##
##======================================================##

def scorealignment(seqa,seqb):
    ## this is a score function, in which continuous gap
    ## will be less punished as a single gap, and the decrease 
    ## is caculated by a decrease factor: α = 0.8
    ## match: 10
    ## mismatch: -10
    ## gap: -8
    ## continuous gap[i]: gap[i-1]*0.8
    if isinstance(seqa,str):
        seqa = list(seqa)
    if isinstance(seqb,str):
        seqb = list(seqb)
    score = 0
    gapscore = -8
    if seqa[0] == seqb[0]:
        score = 10
    else:
        if seqa[0] == "-" or seqb[0] == "-":
            score = -8
        else:
            score = -10
    for i in range(1,len(seqa)):
        if seqa[i] == seqb[i]:
            score = score + 10
        else:
            if seqa[i] != "-" and seqb[i] != "-":
                score = score - 10
            elif seqa[i] == "-" and seqb[i] != "-":
                if seqa[i-1] == "-":
                    gapscore = gapscore*0.8
                else:
                    gapscore = -8
                score = score + gapscore
            elif seqa[i] != "-" and seqb[i] == "-":
                if seqb[i-1] == "-":
                    gapscore = gapscore*0.8
                else:
                    gapscore = -8
                score = score + gapscore
    return round(score,3)

def pairwisealign(seqref,seqtst,method,gapopen=-8):
    ## pairalignment: all is global alignment, but with different algorithm
    ## 1. local: based on Smith-Waterman, and fill the the whole sequence with gap
    ## 2. global: Needleman-Wunsch
    ## 3. gap optimise was added to correct gap location caused by penalty score matrix
    if not (method in ["local","global"]):
        raise SyntaxError("pairwisealign: no available method called " + method)
    if len(seqref) == 0 or len(seqtst) == 0:
        raise SyntaxError("pairwisealign: illegal sequence (no sequence at all)")
    ##=====================================##
    ##====== local paiewise alignment =====##
    ##=====================================##
    if isinstance(seqref,list):
        seqref = "".join(seqref)
    if isinstance(seqtst,list):
        seqtst = "".join(seqtst)
    ##------------------------------------
    ##------------------------------------
    hsseqlo = ctypes.c_char_p(bytes(seqref, 'utf-8'))
    ttseqlo = ctypes.c_char_p(bytes(seqtst, 'utf-8'))
    sumalignlotmp = align.pairwisealign(hsseqlo,ttseqlo,method,gapopen)
    sumalignlo = sumalignlotmp.decode().split("|")
    lohs = sumalignlo[0]
    lott = sumalignlo[1]
    ##=====================================##
    ##======= find excluded sequence ======##
    ##=====================================##
    if method == "local":
        alignref = lohs.replace("-","")
        aligntst = lott.replace("-","")
        ##------------------------------------
        refplen = len(seqref)
        tstplen = len(seqtst)
        refalen = len(alignref)
        tstalen = len(aligntst)
        ##-----find ref excluded sequence-----
        if refplen == refalen:
            refah = ""
            refbe = ""
        else:
            refmapst = seqref.find(alignref)
            refmapen = refmapst + refalen
            refah = seqref[0:refmapst]
            refbe = seqref[refmapen:refplen]
        ##----------------------------------
        if tstplen == tstalen:
            tstah = ""
            tstbe = ""
        else:
            tstmapst = seqtst.find(aligntst)
            tstmapen = tstmapst + tstalen
            tstah = seqtst[0:tstmapst]
            tstbe = seqtst[tstmapen:tstplen]
        ##=====================================##
        ##======= get compensate sequence =====##
        ##=====================================##
        if len(refah) > 0 or len(tstah) > 0:
            if len(refah) > len(tstah):
                comah = "".join(list("-" for _ in range(0,len(refah)-len(tstah))))
                tstah = tstah + comah
            if len(refah) < len(tstah):
                comah = "".join(list("-" for _ in range(0,len(tstah)-len(refah))))
                refah = refah + comah
        ##-----------------------------------     
        if len(refbe) > 0 or len(tstbe) > 0:
            if len(refbe) > len(tstbe):
                combe = "".join(list("-" for _ in range(0,len(refbe)-len(tstbe))))
                tstbe = tstbe + combe
            if len(refbe) < len(tstbe):
                combe = "".join(list("-" for _ in range(0,len(tstbe)-len(refbe))))
                refbe = refbe + combe
        ##-----------------------------------
        correctref = refah+lohs+refbe
        correcttst = tstah+lott+tstbe
    else:
        correctref = lohs
        correcttst = lott
    ##----------------------------------------
    if correctref.replace("-","") != seqref.replace("-",""):
        raise SyntaxError("pairwisealign: ref wrong, nuceotide is lost during optimization")
    if correcttst.replace("-","") != seqtst.replace("-",""):
        raise SyntaxError("pairwisealign: tst wrong, nuceotide is lost during optimization")
    ##---------------------------------------- 
    ##=====================================##
    ##========= correct alignment =========##
    ##=====================================##
    listsum = {"corref":correctref,"cortst":correcttst}
    return listsum

def calIdentity(seqref, seqtst, method, gapopen=-8):
    ## calculate sequence identity
    if not (method in ["local","global"]):
        raise SyntaxError("calIdentity: no available method called " + method)
    if len(seqref) == 0 or len(seqtst) == 0:
        raise SyntaxError("calIdentity: illegal sequence (no sequence at all)")    
    lenseq1 = len(seqref)
    lenseq2 = len(seqtst)
    ##================================================
    hsseqlo = ctypes.c_char_p(bytes(seqref, 'utf-8'))
    ttseqlo = ctypes.c_char_p(bytes(seqtst, 'utf-8'))
    sumalignlotmp = align.pairwisealign(hsseqlo,ttseqlo,method,gapopen)
    ##------------------------------------
    sumalignlo = sumalignlotmp.decode().split("|")
    seq1 = sumalignlo[0]
    seq2 = sumalignlo[1]
    ##============== count exact match ================
    matchcount = 0
    for i in range(0,len(seq1)):
        if seq1[i] == seq2[i]:
            matchcount = matchcount + 1
    ##=================== cauculate ===================
    identmp = round(matchcount/len(seq1),2) if len(seq1) > 0 else 0
    iden = [identmp,round(matchcount/lenseq1,2),round(matchcount/lenseq2,2)]
    
    return iden

##======================================================##
##            pairwise alignment of sequence            ##
##======================================================##

def findgrouppos(ex_len, seqlist, gapingroup = False):
    ## find group pos based on group length and aligned sequence
    if not(gapingroup): ## gap should not be included into neither of the groups
        ecrange = []
        count = 0
        posidx = -1
        pos_st = 0
        ex_idx = 0
        startgrp = True
        for i in range(0,len(seqlist)):
            posidx = posidx + 1
            if seqlist[i] != "-":
                if startgrp:
                    pos_st = posidx
                    startgrp = False
                count = count + 1
                if count == ex_len[ex_idx]:
                    pos_en = posidx
                    ecrange.append([pos_st,pos_en])
                    count = 0
                    startgrp = True
                    ex_idx = ex_idx + 1
        if len(ecrange) != len(ex_len) or ecrange[len(ex_len)-1][1]>=len(seqlist):
            raise SyntaxError("findgrouppos: nucleotide loss during group range screening (no gap)")
    else: ## inter-group gap should be included into the second group, which is required in [optimiseLocal] function
        ecrange = []
        count = 0
        posidx = -1
        pos_st = 0
        ex_idx = 0
        for i in range(0,len(seqlist)):
            posidx = posidx + 1
            if seqlist[i] != "-":
                count = count + 1
                if count == ex_len[ex_idx]:
                    pos_en = posidx
                    ecrange.append([pos_st,pos_en])
                    count = 0
                    pos_st = pos_en + 1
                    ex_idx = ex_idx + 1
                    #break
        if ecrange[len(ex_len)-1][1] > (len(seqlist)-1):
            raise SyntaxError("findgrouppos: nucleotide loss during group range screening (with gap)")
        else: ## confirm bottom gap is included into the group
            ecrange[len(ex_len)-1][1] = len(seqlist)-1
    ##=====================================##
    ##======= get ruler-exon region =======##
    ##=====================================##

    return ecrange

def matchhomoexon(homopair,regr,regt,identhres):
    seqi = listwhich(regr["coln"],"==","seq")[0]
    ## match homo exons
    if len(homopair) == 0:
        newmatch = []
    if len(homopair) == 1:
        newmatch = [[[homopair[0][0]],[homopair[0][1]]]]
    else:
        ## hsa as anchor
        newmatchhsa = []
        sttmp = [homopair[0][0]]
        entmp = [homopair[0][1]]
        for i in range(1,len(homopair)):
            if homopair[i][0] != homopair[i-1][0]:
                newmatchhsa.append([sttmp,entmp])
                sttmp = [homopair[i][0]]
                entmp = [homopair[i][1]]
            else:
                sttmp = union(sttmp,homopair[i][0])
                entmp = union(entmp,homopair[i][1])
        newmatchhsa.append([sttmp,entmp])
        ###=============================================================##
        ###===  add a filter to deal with multiple duplicated exons  ===##
        ###=============================================================##
        if len(newmatchhsa) == 1:
            newmatchtmp = newmatchhsa
        else:
            newmatchtmp = []
            maxttgrp = 0
            for i in range(0,len(newmatchhsa)):
                if len(newmatchhsa[i][1]) == 1:
                    newmatchtmp.append(newmatchhsa[i])
                    maxttgrp = newmatchhsa[i][1][0]
                else:
                    idxtmp = listwhich(newmatchhsa[i][1],">",maxttgrp)
                    grptttmp = selectele(newmatchhsa[i][1],idxtmp)
                    if len(grptttmp) == 1:
                        newmatchtmp.append([newmatchhsa[i][0],grptttmp])
                        maxttgrp = grptttmp[0]
                    if len(grptttmp) > 1:
                        grpdiftt = [1]
                        for j in range(1,len(grptttmp)):
                            grpdiftt.append(grptttmp[j]-grptttmp[j-1])
                        choosettidx = listwhich(grpdiftt,">",1)
                        if len(choosettidx) > 0:
                            stopttidx = min(choosettidx)
                        else:
                            stopttidx = len(grpdiftt)
                        newmatchtmp.append([newmatchhsa[i][0],grptttmp[0:stopttidx]])
                        maxttgrp = min(grptttmp[0:stopttidx])
        ###=============================================================##
        ###=============================================================##
        ## ptr as anchor
        if len(newmatchtmp) == 1:
            newmatch1 = newmatchtmp
        else:
            newmatchtmpd = {"dat":newmatchtmp,"coln":["hsag","ptrg"]}
            newmatch1 = []
            sttmp = newmatchtmp[0][0]
            entmp = newmatchtmp[0][1]
            record = True
            for i in range(1,len(newmatchtmp)):
                if newmatchtmp[i][1] != newmatchtmp[i-1][1]:
                    if record:
                        newmatch1.append([sttmp,entmp])
                        sttmp = newmatchtmp[i][0]
                        entmp = newmatchtmp[i][1]
                    if min(newmatchtmp[i][1]) > max(listasint(selectcol(newmatchtmpd,"ptrg")[0:i],"max")):
                        record = True
                        sttmp = newmatchtmp[i][0]
                        entmp = newmatchtmp[i][1]
                    else:
                        record = False
                else:
                    sttmp = union(sttmp,newmatchtmp[i][0])
                    entmp = union(entmp,newmatchtmp[i][1])
            if record:
                #sttmp = newmatchhsa[i][0]
                #entmp = newmatchhsa[i][1]                
                newmatch1.append([sttmp,entmp])
        ###=============================================================##
        ###  add a filter to deal with non-continuous multiple match    ##
        ###=============================================================##
        newmatch = []
        for i in range(0,len(newmatch1)):
            if len(newmatch1[i][0]) == 1 and len(newmatch1[i][1]) == 1:
                newmatch.append(newmatch1[i])
            elif len(newmatch1[i][0]) == 1 and len(newmatch1[i][1]) > 1:
                grpdiftt = [1]
                for j in range(1,len(newmatch1[i][1])):
                    grpdiftt.append(newmatch1[i][1][j]-newmatch1[i][1][j-1])
                choosettidx = listwhich(grpdiftt,">",1)
                if len(choosettidx) > 0:
                    stopttidx = min(choosettidx)
                else:
                    stopttidx = len(grpdiftt)                
                newmatch.append([newmatch1[i][0],newmatch1[i][1][0:stopttidx]])
            elif len(newmatch1[i][0]) > 1 and len(newmatch1[i][1]) == 1:
                grpdifhs = [1]
                for m in range(1,len(newmatch1[i][0])):
                    grpdifhs.append(newmatch1[i][0][m]-newmatch1[i][0][m-1])
                choosehsidx = listwhich(grpdifhs,">",1)
                if len(choosehsidx) > 0:
                    stophsidx = min(choosehsidx)
                else:
                    stophsidx = len(grpdifhs)
                newmatch.append([newmatch1[i][0][0:stophsidx],newmatch1[i][1]])
            else:
                grpdiftt = [1]
                for j in range(1,len(newmatch1[i][1])):
                    grpdiftt.append(newmatch1[i][1][j]-newmatch1[i][1][j-1])
                choosettidx = listwhich(grpdiftt,">",1)
                if len(choosettidx) > 0:
                    stopttidx = min(choosettidx)
                else:
                    stopttidx = len(grpdiftt)
                ##--------------------------------##
                grpdifhs = [1]
                for m in range(1,len(newmatch1[i][0])):
                    grpdifhs.append(newmatch1[i][0][m]-newmatch1[i][0][m-1])
                choosehsidx = listwhich(grpdifhs,">",1)
                if len(choosehsidx) > 0:
                    stophsidx = min(choosehsidx)
                else:
                    stophsidx = len(grpdifhs)
                ##--------------------------------##
                newmatch.append([newmatch1[i][0][0:stophsidx],newmatch1[i][1][0:stopttidx]])
    ##------ optimise reverse transcription mediated exon fusion ------------##
    optimisereverse = []
    for i in range(0,len(newmatch)):
        if (len(newmatch[i][0]) == 1 and len(newmatch[i][1]) > 1) or (len(newmatch[i][0]) > 1 and len(newmatch[i][1]) == 1):
            optimisereverse.append([newmatch[i][0],newmatch[i][1]])
    ##----------- final optimise ------------## colinearity test
    if len(newmatch) > 0:
        for i in range(0,len(newmatch)):
            if len(newmatch[i][0]) > 1:
                new0 = [newmatch[i][0][0]]
                for j in range(1,len(newmatch[i][0])):
                    if newmatch[i][0][j] - newmatch[i][0][j-1] == 1:
                        new0.append(newmatch[i][0][j])
                    else:
                        newmatch[i][0][j]=newmatch[i][0][j-1]
                newmatch[i][0] = new0
            if len(newmatch[i][1]) > 1:
                new1 = [newmatch[i][1][0]]
                for j in range(1,len(newmatch[i][1])):
                    if newmatch[i][1][j] - newmatch[i][1][j-1] == 1:
                        new1.append(newmatch[i][1][j])
                    else:
                        newmatch[i][1][j]=newmatch[i][1][j-1]
                newmatch[i][1] = new1
    newmatchtmp = copy.deepcopy(newmatch)
    ##---------------------------------##
    dropi = []
    for i in range(0,len(newmatch)):
        if newmatch[i] in optimisereverse:
            continue
        if newmatch[i][0] != [0] and newmatch[i][1] != [0]:
            hsseqtmp = ""
            for j in range(0,len(newmatch[i][0])):
                hsseqtmp = hsseqtmp + regr["dat"][listwhich(selectcol(regr,"group"),"==",newmatch[i][0][j])[0]][seqi]
        ##-------------------------------------##
            ttseqtmp = ""
            for j in range(0,len(newmatch[i][1])):
                ttseqtmp = ttseqtmp + regt["dat"][listwhich(selectcol(regt,"group"),"==",newmatch[i][1][j])[0]][seqi]
        if calIdentity(hsseqtmp, ttseqtmp, method="local")[0] < identhres:
            dropi.append(i)
    if len(dropi) > 0:
        newmatch = listdrop(newmatch,dropi)

    ## drop repeat exon induced match  ##
    return newmatch

##======================================================##
##            pairwise alignment of exon                ##
##======================================================##
#score=2; penalty=-2; gap=-1
def alignmentexon(regr, regt, identhres=0.8, coverthres=0.8, microexon=10, score=2, penalty=-2, gap=-1):  ##Neuron-Specific Alternative Splicing of Transcriptional Machineries: Implications for Neurodevelopmental Disorders
    ##-------------------------------------##
    seqi = listwhich(regr["coln"],"==","seq")[0]
    leni = listwhich(regr["coln"],"==","len")[0]
    ##==========================================================##
    ##  pairwise alignment among exons, this step is necessary  ##
    ##  to esitmate 1-N or N-1 homologous exons                 ##
    ##==========================================================##
    #print( "local pairwise alignment of " + str(len(regr["dat"])) + " human and " + str(len(regt["dat"])) + " non-human exon regions")
    ##-------------------------------------##
    seq2 = selectcol(regr,"group")
    seq2.sort(reverse = True)
    ##-------------------------------------##
    seq1 = selectcol(regt,"group")
    seq1.sort(reverse = True)
    ##-------------------------------------##
    ##       build identity matrix         ##
    ##-------------------------------------##
    exnim = [[0.0]*(len(seq2)+1) for _ in range(0,len(seq1)+1)]
    for i in range(0,len(regt["dat"])):
        ri = len(regt["dat"]) - i -1
        hsgrpseq = regt["dat"][ri][seqi]
        for j in range(0,len(regr["dat"])):
            rj = len(regr["dat"]) - j -1
            ttgrpseq = regr["dat"][rj][seqi]
            identmp = calIdentity(hsgrpseq,ttgrpseq,method="local")
            if identmp[1] >= coverthres or identmp[2] >= coverthres:
                exnim[i+1][j+1] = identmp[0]
    ##------------------------------------------##
    ##     build score and arrow matrix         ##
    ## 1. match(identity>=0.8,coverage>=0.8): 6 ##
    ## 2. mismatch: -6                          ##
    ## 3. exon jump, exon split or fusion: -3   ##
    ##------------------------------------------##
    scorem = [[0]*(len(seq2)+1) for _ in range(0,len(seq1)+1)]
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if exnim[i][j] < identhres:
                scorei = penalty  ## mismatch penalty
            else:
                scorei = score   ## match score
            smat = [0,0,0]
            smat[0] = scorem[i-1][j] + scorei + gap ## position penalty
            smat[1] = scorem[i][j-1] + scorei + gap ## position penalty
            smat[2] = scorem[i-1][j-1] + scorei
            scorem[i][j] = max(smat)
    ##==========================================================##
    ##  test whether global or local alignment                  ##
    ##==========================================================##
    if exnim[len(exnim)-1][len(exnim[0])-1] >= identhres:
        ##-------------------------------------##
        ## global alignment  
        h=len(seq1)
        v=len(seq2)
    else:
        ## local alignment
        max_score = scorem[0][0]
        for i in range(0,len(seq1)+1):
            for j in range(1,len(seq2)+1):
                if max_score <= scorem[i][j]:
                    max_score = copy.deepcopy(scorem[i][j])
                    h = i
                    v = j
    ##==========================================================##
    ##==========================================================##
    ##==========================================================##
    ##  dynamic programming to find best exon matches           ##
    ##==========================================================##
    alnseq1 = []
    alnseq2 = []
    ##-------------------------------------##
    seq1_num = h-1
    seq2_num = v-1
    if exnim[h][v] >= identhres:
        alnseq1.append(seq1[seq1_num])
        alnseq2.append(seq2[seq2_num])

    while True:
        if v == 0 or h == 0:
            break

        score0 = scorem[h-1][v]
        score1 = scorem[h][v-1]
        score2 = scorem[h-1][v-1]

        maxscore = max([score0,score1,score2])

        if score2 == maxscore:
            seq1_num -= 1
            seq2_num -= 1
            v -= 1
            h -= 1
            if exnim[h][v] >= identhres:
                alnseq1.append(seq1[seq1_num])
                alnseq2.append(seq2[seq2_num])

        elif score0 == maxscore:
            seq1_num -= 1
            h -= 1
            if exnim[h][v] >= identhres:
                alnseq1.append(seq1[seq1_num])
                alnseq2.append(seq2[seq2_num])

        elif score1 == maxscore:
            seq2_num -= 1
            v -= 1
            if exnim[h][v] >= identhres:
                alnseq1.append(seq1[seq1_num])
                alnseq2.append(seq2[seq2_num])       

    ##==========================================================##
    ##               drop microexon to avoid mismatch           ##
    ##==========================================================##
    ##==========================================================##
    ##  find microexon                                          ##
    ##==========================================================##
    grpi = listwhich(regr["coln"],"==","group")[0]
    leni = listwhich(regr["coln"],"==","len")[0]
    hsreglen = []
    hsreggrp = []
    for i in range(0,len(regr["dat"])):
        hsreggrp.append(regr["dat"][i][grpi])
        hsreglen.append(regr["dat"][i][leni])
    dropgrphs = selectele(hsreggrp,listwhich(hsreglen,"<=",microexon))
    ##------------------------------------------##
    ttreglen = []
    ttreggrp = []
    for i in range(0,len(regt["dat"])):
        ttreggrp.append(regt["dat"][i][grpi])
        ttreglen.append(regt["dat"][i][leni])
    dropgrptt = selectele(ttreggrp,listwhich(ttreglen,"<=",microexon))
    ##==========================================================##
    ##  construct possible orthologous exon relationship        ##
    ##==========================================================##
    homoexonpair = []
    for i in range(0,len(alnseq1)):
        if alnseq2[i] in dropgrphs or alnseq1[i] in dropgrptt:
            continue
        else:
            homoexonpair.append([alnseq2[i],alnseq1[i]])
    ##==========================================================##
    ## optimise multiple to one or one to multiple relationships #
    ##==========================================================##
    # hs as anchor
    newmatchhsa = []
    if len(homoexonpair) == 1:
        newmatchhsa = [[[homoexonpair[0][0]],[homoexonpair[0][1]]]]
    elif len(homoexonpair) > 1:
        sttmp = [homoexonpair[0][0]]
        entmp = [homoexonpair[0][1]]
        for i in range(1,len(homoexonpair)):
            if homoexonpair[i][0] != homoexonpair[i-1][0]:
                newmatchhsa.append([sttmp,entmp])
                sttmp = [homoexonpair[i][0]]
                entmp = [homoexonpair[i][1]]
            else:
                sttmp = union(sttmp,homoexonpair[i][0])
                entmp = union(entmp,homoexonpair[i][1])
        newmatchhsa.append([sttmp,entmp])
    #pt as anchor
    possiblehomo = []
    if len(newmatchhsa) == 1:
        possiblehomo = copy.deepcopy(newmatchhsa)
    elif len(newmatchhsa) > 1:  
        sttmp = newmatchhsa[0][0]
        entmp = newmatchhsa[0][1]
        for i in range(1,len(newmatchhsa)):
            if newmatchhsa[i][1] != newmatchhsa[i-1][1]:
                possiblehomo.append([sttmp,entmp])
                sttmp = newmatchhsa[i][0]
                entmp = newmatchhsa[i][1]
            else:
                sttmp = union(sttmp,newmatchhsa[i][0])
                entmp = union(entmp,newmatchhsa[i][1])
        possiblehomo.append([sttmp,entmp])
    ##==========================================================##
    ##  got final orthologous exons                             ##
    ##==========================================================##
    ffhomopair = []
    for i in range(0,len(possiblehomo)):
        if len(possiblehomo[i][0]) == 1 and len(possiblehomo[i][1]) == 1:
            ffhomopair.append(possiblehomo[i])
        else:
            hsgrp = possiblehomo[i][0]
            ttgrp = possiblehomo[i][1]
            ##---------------------------##
            hsseqtmp = ""
            for hsg in hsgrp:
                hsgi = listwhich(selectcol(regr,"group"),"==",hsg)[0]
                hsseqtmp = hsseqtmp + selectcol(regr,"seq")[hsgi]
            ttseqtmp = ""
            for ttg in ttgrp:
                ttgi = listwhich(selectcol(regt,"group"),"==",ttg)[0]
                ttseqtmp = ttseqtmp + selectcol(regt,"seq")[ttgi]
            ##------------------
            wholealign = pairwisealign(hsseqtmp,ttseqtmp,method="local",gapopen=-3)
            scorewtest = scorealignment(wholealign["corref"],wholealign["cortst"])                    
            ##---------------------------##
            scorelmax = []
            grptmp = []
            for hsg in hsgrp:
                hsgi = listwhich(selectcol(regr,"group"),"==",hsg)[0]
                hsseqtmp = selectcol(regr,"seq")[hsgi]
                for ttg in ttgrp:
                    ttgi = listwhich(selectcol(regt,"group"),"==",ttg)[0]
                    ttseqtmp = selectcol(regt,"seq")[ttgi]
                    grptmp.append([[hsg],[ttg]])
                    localalign = pairwisealign(hsseqtmp,ttseqtmp,method="local",gapopen=-3)
                    scoreltest = scorealignment(localalign["corref"],localalign["cortst"])                                     
                    scorelmax.append(scoreltest)
                
            selecttari = listwhich(scorelmax,">=",scorewtest) ## judge whether the exon map is caused by repeat exon
            if len(selecttari) > 0:
                ffhomopair = ffhomopair + selectele(grptmp,[listwhich(scorelmax,"==",max(scorelmax))[0]])
            else:
                ffhomopair.append(possiblehomo[i])
    homoexongrp = {"dat":copy.deepcopy(ffhomopair),"coln":["hsag","ptrg"]}
    ##==========================================================##
    ##  retrive lost exons in multiple-one or one-multiple      ##
    ##  relationships because of microexon or lower identity    ##
    ##==========================================================##
    if len(homoexongrp["dat"]) > 0:
        for i in range(0,len(homoexongrp["dat"])):
            if len(homoexongrp["dat"][i][0]) > 1:
                homoexongrp["dat"][i][0] = list(range(min(homoexongrp["dat"][i][0]),max(homoexongrp["dat"][i][0])+1))
            if len(homoexongrp["dat"][i][1]) > 1:
                homoexongrp["dat"][i][1] = list(range(min(homoexongrp["dat"][i][1]),max(homoexongrp["dat"][i][1])+1))
    ##==========================================================##
    ##  return final results                                    ##
    ##==========================================================##   
    return homoexongrp

##======================================================##
##     orthology and colinearity test of exon           ##
##======================================================##

def sortgroup(grouppair, anchor):
    ## resort group pairs, this is important to export ordered homo-exon pairs
    ## re-organize group         
    newpaird = {"dat":grouppair,"coln":["hsa","ptr"]}
    if len(anchor) == 0:
        newpair = grouppair
    else:
        newpair = []
        for ii in range(0,len(anchor)):
            ## test whether list-list (1:N, N-1, N-N) exist
            if ii == 0:
                if isinstance(anchor[ii][0],list):
                    hsamax = min(anchor[ii][0])
                    ptrmax = min(anchor[ii][1])
                    hsamin = min(anchor[ii][0])
                    ptrmin = min(anchor[ii][1])
                else:
                    hsamax = anchor[ii][0]
                    ptrmax = anchor[ii][1]
                    hsamin = anchor[ii][0]
                    ptrmin = anchor[ii][1]                        
            #if ii > 0:
            else:
                if isinstance(anchor[ii][0],list):
                    hsamax = min(anchor[ii][0])
                    ptrmax = min(anchor[ii][1])
                else:
                    hsamax = anchor[ii][0]
                    ptrmax = anchor[ii][1]
                ##----------------------------------------
                if isinstance(anchor[ii-1][0],list):
                    hsamin = max(anchor[ii-1][0])
                    ptrmin = max(anchor[ii-1][1])
                else:
                    hsamin = anchor[ii-1][0]
                    ptrmin = anchor[ii-1][1]                
            ##----------------------------------------
            if ii == 0:
                hsii = union(listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"<",hsamax),listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"==",0))
                ttii = union(listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"<",ptrmax),listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"==",0))
                selecti = intersect(hsii, ttii)
                sttmp = selectrow(newpaird,selecti)
                newpair = newpair + listsort(sttmp,by=["hsa","ptr"],increase=True)["dat"]
                newpair = newpair + [anchor[ii]]
            else:
                hsii1 = intersect(listwhich(listasint(selectcol(newpaird,"hsa"),deal="min"),">",hsamin),listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"<",hsamax))
                ttii1 = intersect(listwhich(listasint(selectcol(newpaird,"ptr"),deal="min"),">",ptrmin),listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"<",ptrmax))
                hsii = union(hsii1,listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"==",0))
                ttii = union(ttii1,listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"==",0))
                selecti = intersect(hsii, ttii)
                sttmp = selectrow(newpaird,selecti)
                newpair = newpair + listsort(sttmp,by=["hsa","ptr"],increase=True)["dat"]
                newpair = newpair + [anchor[ii]]
            ##----------------------------------------
            if ii == len(anchor) - 1:
                hsii = union(listwhich(listasint(selectcol(newpaird,"hsa"),deal="min"),">",hsamax),listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"==",0))
                ttii = union(listwhich(listasint(selectcol(newpaird,"ptr"),deal="min"),">",ptrmax),listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"==",0))
                selecti = intersect(hsii, ttii)
                sttmp = selectrow(newpaird,selecti)
                newpair = newpair + listsort(sttmp,by=["hsa","ptr"],increase=True)["dat"]
    ##----------------------------------------                     
    return newpair

def findnonhomoregion(homogrp):
    ## find local groups that are required tooptimise
    if len(homogrp["dat"]) > 1:
        total = []
        localhs = []
        localpt = []
        if homogrp["dat"][0][0] == 0:
            localpt.append(homogrp["dat"][0][1])
        if homogrp["dat"][0][1] == 0:
            localhs.append(homogrp["dat"][0][0])
        for i in range(1,len(homogrp["dat"])):
            if homogrp["dat"][i][0] == 0:
                if homogrp["dat"][i-1][0] == 0 or homogrp["dat"][i-1][1] == 0:
                    localpt.append(homogrp["dat"][i][1])

            elif homogrp["dat"][i][1] == 0:
                if homogrp["dat"][i-1][1] == 0 or homogrp["dat"][i-1][0] == 0:
                    localhs.append(homogrp["dat"][i][0])
            
            elif homogrp["dat"][i][0] != 0 and homogrp["dat"][i][1] != 0:
                if len(localhs) > 0 and len(localpt) > 0:
                    total.append([localhs, localpt])
                localhs = []
                localpt = []

        if len(localhs) > 0 and len(localpt) > 0:
            total.append([localhs, localpt])               

    else:
        total = []

    return unique(total)

def findlocalmatch(regr,regt,cdshs,cdstt,identhres):
    hsgrp = selectcol(regr,"group")
    ttgrp = selectcol(regt,"group")
    ##------------------------------------------------------##
    ##------------------------------------------------------##
    idenml = []  ##
    #etyi = listwhich(regr["coln"],"==","etype")[0]
    epsti = listwhich(regr["coln"],"==","start")[0]
    epeni = listwhich(regr["coln"],"==","end")[0]
    cdsseqi = listwhich(cdshs["coln"],"==","seq")[0]
    for hsg in hsgrp:
        hsgi = listwhich(selectcol(regr,"group"),"==",hsg)[0]
        for ttg in ttgrp:
            ttgi = listwhich(selectcol(regt,"group"),"==",ttg)[0]
            hseidx = intersect(listwhich(selectcol(cdshs,"start"),">=",regr["dat"][hsgi][epsti]),listwhich(selectcol(cdshs,"end"),"<=",regr["dat"][hsgi][epeni]))
            tteidx = intersect(listwhich(selectcol(cdstt,"start"),">=",regt["dat"][ttgi][epsti]),listwhich(selectcol(cdstt,"end"),"<=",regt["dat"][ttgi][epeni]))
            idenme = []
            idenmmax = []            
            for jjj in hseidx:
                hsseqtmp = cdshs["dat"][jjj][cdsseqi]
                for iii in tteidx:
                    ttseqtmp = cdstt["dat"][iii][cdsseqi]
                    #seqref=hsseqtmp; seqtst=ttseqtmp
                    seqidenltmp = calIdentity(hsseqtmp,ttseqtmp,method="local")
                    idenme.append(seqidenltmp)
                    idenmmax.append(round(seqidenltmp[0]*seqidenltmp[1]*seqidenltmp[2],2))
            
            selecttmpi = listwhich(idenmmax,"==",max(idenmmax))[0]
            idenml.append([hsg,ttg] + idenme[selecttmpi])
    ##------------------------------------------------------##
    groupviol = []
    for i in range(0,len(idenml)):
        if idenml[i][2] >= identhres and (idenml[i][3] >= identhres or idenml[i][4] >= identhres):
            groupviol.append([idenml[i][0],idenml[i][1]])
    ##==============================================##
    ## filter
    if len(groupviol) > 1:
        dropi = []
        groupviold = {"dat":groupviol,"coln":["hsag","ptrg"]}
        for i in range(1,len(groupviol)):
            if not(groupviol[i][1] >= max(listasint(selectcol(groupviold,"ptrg")[0:i],"max"))):
                dropi = union(dropi,i)
        if len(dropi) > 0:
            groupviol = listdrop(groupviol,dropi)
    groupviol = unique(groupviol)
    ##==============================================##
    ## find idenm
    idenmld = {"dat":idenml,"coln":["hsag","ptrg","iden","hsacov","ptrcov"]}
    idenmlnew = []
    for i in range(0,len(groupviol)):
        selecti = intersect(listwhich(selectcol(idenmld,"hsag"),"==",groupviol[i][0]),listwhich(selectcol(idenmld,"ptrg"),"==",groupviol[i][1]))[0]
        idenmlnew.append(idenml[selecti])
    idenml = {"dat":copy.deepcopy(idenmlnew),"coln":["hsag","ptrg","iden","hsacov","ptrcov"]}
    ##======= detect possible 1-1 groups, which is the reference list if 1-N, N-1 even N-N failed
    ##------------------------------------------------------##
    ##------------------------------------------------------##
    if len(groupviol) > 0:
        seqi = listwhich(regr["coln"],"==","seq")[0]
        #homopair=groupviol
        groupvio = matchhomoexon(groupviol,regr,regt,identhres)
        ##unique(groupviol) ###*****************************
        unigroupviol = []
        if len(groupvio) > 0:
            for i in range(0,len(groupvio)):
                if not(isinstance(groupvio[i][0],list)):
                    unigroupviol.append([[groupviol[i][0]],[groupviol[i][1]]])
                else:
                    #idenmtest = []
                    idenmtestmax = []
                    testgrp = []
                    for hsgtmp in groupvio[i][0]:
                        for ttgtmp in groupvio[i][1]:
                            testgrp.append([hsgtmp,ttgtmp])
                            selectideni = intersect(listwhich(selectcol(idenml,"hsag"),"==",hsgtmp),listwhich(selectcol(idenml,"ptrg"),"==",ttgtmp))[0]
                            #idenmtest.append(idenml["dat"][selectideni])
                            idenmtestmax.append(round(idenml["dat"][selectideni][2]*idenml["dat"][selectideni][3]*idenml["dat"][selectideni][4],2))    
                    targrp = testgrp[listwhich(idenmtestmax,"==",max(idenmtestmax))[0]]
                    unigroupviol.append([[targrp[0]],[targrp[1]]])
        unigroupviold = {"dat":unigroupviol,"coln":["hsag","ptrg"]}
        ##------------- reorganize ------------
        eletest = []
        for i in range(0,len(groupvio)):
            eletest.append(groupvio[i] in unigroupviol)
        if len(listwhich(eletest,"==",True)) == len(groupvio):
            groupvio = copy.deepcopy(unigroupviol)
        else:
            existst = []
            existen = []
            for j in range(0,len(groupvio)):
                existst = existst + groupvio[j][0]
                existen = existen + groupvio[j][1]
            selecti = []
            for j in range(0,len(unigroupviol)):
                if not(unigroupviol[j][0][0] in existst) and not(unigroupviol[j][1][0] in existen):
                    selecti.append(j)
            if len(selecti) > 0:
                groupvio = groupvio + selectele(unigroupviol,selecti)
                groupvio.sort()
        ##------------- retest ------------
        newhomopair = []
        for i in range(0,len(groupvio)):
            if len(groupvio[i][0]) == 1 and len(groupvio[i][1]) == 1:
                newhomopair.append(groupvio[i])
            else:
                ##===== calculate group to group identity
                identmpl = []
                for j in range(0,len(groupvio[i][0])):
                    for k in range(0,len(groupvio[i][1])):
                        identmp = calIdentity(regr["dat"][listwhich(selectcol(regr,"group"),"==",groupvio[i][0][j])[0]][seqi],regt["dat"][listwhich(selectcol(regt,"group"),"==",groupvio[i][1][k])[0]][seqi],method="local")
                        identmpl.append(identmp[0]*identmp[1]*identmp[2])
                ##----- calculate merged group identity
                hsseq = ""
                for j in range(0,len(groupvio[i][0])):
                    hsseq = hsseq + regr["dat"][listwhich(selectcol(regr,"group"),"==",groupvio[i][0][j])[0]][seqi]
                ttseq = ""
                for j in range(0,len(groupvio[i][1])):
                    ttseq = ttseq + regt["dat"][listwhich(selectcol(regt,"group"),"==",groupvio[i][1][j])[0]][seqi]
                identmp = calIdentity(hsseq,ttseq,method="local")
                ##------ judge
                dropidx = listwhich(identmpl,">",identmp[0]*identmp[1]*identmp[2]) ## judge whether the exon map is caused by repeat exon
                if len(dropidx) == 0:
                    newhomopair.append(groupvio[i])
                else: ## only keep the first one
                    ## find which unigroupviol contains target group, choose the group pair with maximum identity
                    selectunii = []
                    for hsgtmp in groupvio[i][0]:
                        for ttgtmp in groupvio[i][1]:
                            selectideni = intersect(listwhich(selectcol(unigroupviold,"hsag"),"==",[hsgtmp]),listwhich(selectcol(unigroupviold,"ptrg"),"==",[ttgtmp]))
                            selectunii = union(selectunii,selectideni)
                    if len(selectunii) > 0:
                        newhomopair.append(unigroupviold["dat"][selectunii[0]])
    else:
        newhomopair = []
    return newhomopair

def colinearity_test(homopair):
    if len(homopair) == 0:
        passhomo2 =  {"dat":[],"coln":["hsag","ptrg"]}
    else:
        homogdict = {"dat":unique(homopair),"coln":["hsag","ptrg"]}
        homogdict = listsort(homogdict,by=["hsag","ptrg"],increase=True)
        ##===========================================
        ##===========================================
        ##--------------- the start group: homopair[0] --------------
        ##===========================================
        minttg = homogdict["dat"][0][1]
        minhsg = homogdict["dat"][0][0]
        passhomo = {"dat":[[minhsg,minttg]],"coln":["hsag","ptrg"]}
        ttgtmp = homogdict["dat"][0][1]
        ##-- find the hsa group that need to be optimised ---
        unighs = unique(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"hsag"),">=",minhsg)))
        for i in range(1,len(unighs)):
            ptrgtmp = selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i]))
            #minttg = max(selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i])))
            maxttgminus = max(selectele(selectcol(passhomo,"ptrg"),listwhich(selectcol(passhomo,"hsag"),"<",unighs[i])))
            ptrtari = listwhich(ptrgtmp,">",maxttgminus)
            if len(ptrtari) > 0:
                ttgtmptmp = ptrgtmp[min(ptrtari)]
                if len(ttgtmptmp) != 1 or len(ttgtmp) != 1:
                    ttgtmp = setdiff(ttgtmptmp,ttgtmp)
                else:
                    ttgtmp = copy.deepcopy(ttgtmptmp)
                
                if len(ttgtmp) > 0:
                    passhomo["dat"].append([unighs[i],ttgtmp])                

        ##===========================================
        ##===========================================
        ##--------------- the start group: [minhsg,minttg] --------------
        ##===========================================
        minttg = min(selectcol(homogdict,"ptrg"))
        minhsg = min(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"ptrg"),"==",minttg)))
        passhomo1 = {"dat":[[minhsg,minttg]],"coln":["hsag","ptrg"]}
        ttgtmp = copy.deepcopy(minttg)
        ##-- find the hsa group that need to be optimised ---
        unighs = unique(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"hsag"),">=",minhsg)))
        for i in range(1,len(unighs)):
            ptrgtmp = selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i]))
            #minttg = max(selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i])))
            maxttgminus = max(selectele(selectcol(passhomo1,"ptrg"),listwhich(selectcol(passhomo1,"hsag"),"<",unighs[i])))
            ptrtari = listwhich(ptrgtmp,">",maxttgminus)
            if len(ptrtari) > 0:
                ttgtmptmp = ptrgtmp[min(ptrtari)]
                if len(ttgtmptmp) != 1 or len(ttgtmp) != 1:
                    ttgtmp = setdiff(ttgtmptmp,ttgtmp)
                else:
                    ttgtmp = copy.deepcopy(ttgtmptmp)
                
                if len(ttgtmp) > 0:
                    passhomo1["dat"].append([unighs[i],ttgtmp])
        
        ##===========================================
        ## choose the best matches
        ##===========================================
        if len(passhomo1["dat"]) > len(passhomo["dat"]):
            passhomo2 = passhomo1
        else:
            passhomo2 = passhomo
    ##===========================================
      
    return passhomo2

def colinearity_testold(homopair):
    if len(homopair) == 0:
        passhomo2 =  {"dat":[],"coln":["hsag","ptrg"]}
    else:
        homogdict = {"dat":unique(homopair),"coln":["hsag","ptrg"]}
        homogdict = listsort(homogdict,by=["hsag","ptrg"],increase=True)
        ##--------------- find the start group --------------
        ##===========================================
        minttg = homogdict["dat"][0][1]
        minhsg = homogdict["dat"][0][0]
        passhomo = {"dat":[[minhsg,minttg]],"coln":["hsag","ptrg"]}
        ##-- find the hsa group that need to be optimised ---
        unighs = unique(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"hsag"),">",minhsg)))
        for i in range(0,len(unighs)):
            ptrgtmp = selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i]))
            #minttg = max(selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i])))
            maxttgminus = max(selectele(selectcol(passhomo,"ptrg"),listwhich(selectcol(passhomo,"hsag"),"<",unighs[i])))
            ptrtari = listwhich(ptrgtmp,">",maxttgminus)
            if len(ptrtari) > 0:
                passhomo["dat"].append([unighs[i],ptrgtmp[min(ptrtari)]])

        ##===========================================
        minttg = min(selectcol(homogdict,"ptrg"))
        minhsg = min(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"ptrg"),"==",minttg)))
        passhomo1 = {"dat":[[minhsg,minttg]],"coln":["hsag","ptrg"]}
        ##-- find the hsa group that need to be optimised ---
        unighs = unique(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"hsag"),">",minhsg)))
        for i in range(0,len(unighs)):
            ptrgtmp = selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i]))
            #minttg = max(selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i])))
            maxttgminus = max(selectele(selectcol(passhomo1,"ptrg"),listwhich(selectcol(passhomo1,"hsag"),"<",unighs[i])))
            ptrtari = listwhich(ptrgtmp,">",maxttgminus)
            if len(ptrtari) > 0:
                passhomo1["dat"].append([unighs[i],ptrgtmp[min(ptrtari)]])   

        if len(passhomo1["dat"]) > len(passhomo["dat"]):
            passhomo2 = passhomo1
        else:
            passhomo2 = passhomo
    ##===========================================
      
    return passhomo2

def mergehomogroup(homog1,homog2):
    if len(homog1["dat"]) == 0 and len(homog2["dat"]) == 0:
        homogrptmp = {"dat":[],"coln":["hsag","ptrg"]}
    elif len(homog1["dat"]) == 0:
        homogrptmp = copy.deepcopy(homog2)
    elif len(homog2["dat"]) == 0:
        homogrptmp = copy.deepcopy(homog1)
    else:
        selecti = [] ## for homog2
        dropi = []   ## for homog1
        for i in range(0,len(homog2["dat"])):
            rig = homog2["dat"][i][0]
            leg = homog2["dat"][i][1]
            idxtmp = []
            storenon = []
            storeexi = []
            for j in range(0,len(rig)):
                for m in range(0,len(leg)):
                    idxtmptmp = intersect(listwhich(selectcol(homog1,"hsag"),"==",[rig[j]]),listwhich(selectcol(homog1,"ptrg"),"==",[leg[m]]))
                    if len(idxtmptmp) > 0 :
                        idxtmp.append(idxtmptmp[0])
                    ## test if possible to replace
                    if not(([rig[j]] in selectcol(homog1,"hsag")) or ([leg[m]] in selectcol(homog1,"ptrg")) ):
                        storenon.append(True)
                    else:
                        storenon.append(False)
                    ## test if possible to add non-record
                    if not(([rig[j]] in selectcol(homog1,"hsag")) and ([leg[m]] in selectcol(homog1,"ptrg")) ):
                        storeexi.append(True)
                    else:
                        storeexi.append(False)
                    
            if len(idxtmp) > 0 and len(listwhich(storeexi,"==",True)) == len(storeexi)-1:
                dropi.append(idxtmp[0])
                selecti.append(i)
            elif len(idxtmp) == 0 and len(listwhich(storenon,"==",True)) == len(storenon):
                selecti.append(i)
        homog = []
        homog1["dat"] = listdrop(homog1["dat"],dropi)
        homog2["dat"] = selectele(homog2["dat"],selecti)
        homog = homog + homog1["dat"] + homog2["dat"]
        homogrptmp = colinearity_test(homog)

    return homogrptmp

def homo_colinearity_test(species1, species2,regrhs,regrtt,cdshs,cdstt,homocdsdat,identhres=0.8, coverthres=0.8, minexon=10, mapscore=2, misscore=-2, gapscore=-1):
    ## colinearity test with a homologous exon file
    ##=====================================##
    ##============= find index ============##
    ##=====================================##
    hsai = listwhich(homocdsdat["coln"],"==",species1)[0]
    ptri = listwhich(homocdsdat["coln"],"==",species2)[0]
    idsti = listwhich(cdshs["coln"],"==","start")[0]
    ideni = listwhich(cdshs["coln"],"==","end")[0]
    grpi = listwhich(regrhs["coln"],"==","group")[0]
    ##==================================================================##
    ##===== construct orthologous exon matrix based on exon BLASTN =====##
    ##==================================================================##
    if len(homocdsdat["dat"]) > 0:
        #print("homologous exons pass blastN")
        homogroup = []
        for hci in range(0,len(homocdsdat["dat"])):
            hsae = homocdsdat["dat"][hci][hsai]
            ptre = homocdsdat["dat"][hci][ptri]
            ##----------------------------------------
            hsast = cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hsae)[0]][idsti]
            hsaen = cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hsae)[0]][ideni]
            idxhs = intersect(listwhich(selectcol(regrhs,"start"),"<=",hsast),listwhich(selectcol(regrhs,"end"),">=",hsaen))[0]
            hsgrp = regrhs["dat"][idxhs][grpi]
            ##----------------------------------------
            ptrst = cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",ptre)[0]][idsti]
            ptren = cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",ptre)[0]][ideni]
            idxtt = intersect(listwhich(selectcol(regrtt,"start"),"<=",ptrst),listwhich(selectcol(regrtt,"end"),">=",ptren))[0]
            ttgrp = regrtt["dat"][idxtt][grpi]
            ##----------------------------------------
            homogroup.append([[hsgrp],[ttgrp]])
        homogrp1 = {"dat":unique(homogroup),"coln":["hsag","ptrg"]}
        homogrp1 = listsort(homogrp1,by=["hsag","ptrg"],increase=True)
        ##== colinearity_test of conserved exons that are called by BLASTN =##
        homogrp1 = colinearity_test(homogrp1["dat"])
        ##=====================================================================##
    else:
        homogrp1 = {"dat":[],"coln":["hsag","ptrg"]}
    ##==================================================================##
    ## construct orthologous exon matrix based on exon based alignment  ##
    ##==================================================================##
    #print("homologous exons estimated using exon alignment")
    #regr=regrhs; regt=regrtt; identhres=0.6; minexon=2; coverthres=0.7
    homogrp2 = alignmentexon(regr=regrhs, regt=regrtt, microexon=minexon, identhres=identhres, coverthres=coverthres, score=mapscore, penalty=misscore, gap=gapscore)
    homogrp2tmp = copy.deepcopy(homogrp2)
    homogrp2tmptmp = colinearity_test(homogrp2["dat"])
    while homogrp2tmptmp != homogrp2tmp:
        homogrp2tmp = copy.deepcopy(homogrp2tmptmp)
        homogrp2tmptmp = colinearity_test(homogrp2tmptmp["dat"])
    homogrp2 = copy.deepcopy(homogrp2tmptmp)
    ##==================================================================##
    ## correct the dynamic programming method using BLASTN based results##
    ##==================================================================##
    homogrp = mergehomogroup(homogrp1,homogrp2)
    ##==================================================================##      
    ##----------- organize orthologous exon group -------------##
    for i in range(0,len(homogrp["dat"])):
        if len(homogrp["dat"][i][0]) == 1 and len(homogrp["dat"][i][1]) == 1:
            homogrp["dat"][i][0] = homogrp["dat"][i][0][0]
            homogrp["dat"][i][1] = homogrp["dat"][i][1][0]
    anchor = copy.deepcopy(homogrp["dat"])
    ##==========================================================##
    ##============ add non record exonsinto homo list ==========##
    ##==========================================================##
    #print("optimise ruler alignment based on homologous exon group")
    hsgroup = selectcol(regrhs,"group")
    ttgroup = selectcol(regrtt,"group")
    for hsg in hsgroup:
        if not(hsg in listasint(selectcol(homogrp,"hsag"),deal="all")):
            if len(homogrp["dat"])>0 and hsg < max(listasint(selectcol(homogrp,"hsag"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"hsag"),deal="min"),">",hsg)),[hsg,0])
            else:
                homogrp["dat"].append([hsg,0])
    for ttg in ttgroup:
        if not(ttg in listasint(selectcol(homogrp,"ptrg"),deal="all")):
            if len(homogrp["dat"])>0 and ttg < max(listasint(selectcol(homogrp,"ptrg"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"ptrg"),deal="min"),">",ttg)),[0,ttg])
            else:
                homogrp["dat"].append([0,ttg])
    ##==========================================================##
    ## sort group in standard order so that findnonhomoregion   ##
    ## can estimate correct position range                      ##
    ##==========================================================##
    homogrp["dat"] = sortgroup(grouppair=homogrp["dat"], anchor=anchor)
    ##=============================================================##
    ## find non-orthologous exon that are needed violent alignment ##
    ##=============================================================##
    total = findnonhomoregion(homogrp)
    ##-------------------------------------##
    localmatch = []
    if len(total) > 0:
        #print("detecting orthologous group that are partial matched")
        for ii in range(0,len(total)):
            hsgi = []
            for j in total[ii][0]:
                hsgi.append(listwhich(selectcol(regrhs,"group"),"==",j)[0])
            ptgi = []
            for j in total[ii][1]:
                ptgi.append(listwhich(selectcol(regrtt,"group"),"==",j)[0])
            
            #regr=selectrow(regrhs,hsgi); regt=selectrow(regrtt,ptgi);cdshs=cdshs;cdstt=cdstt
            violentalign = findlocalmatch(regr=selectrow(regrhs,hsgi), regt=selectrow(regrtt,ptgi),cdshs=cdshs,cdstt=cdstt,identhres=identhres)
            localmatch = localmatch + violentalign
    ##---------- drop and loop to insert ---------## 
    if len(localmatch) > 0:
        for i in range(0,len(localmatch)):
            if len(localmatch[i][0]) == 1 and len(localmatch[i][1]) == 1:
                localmatch[i][0] = localmatch[i][0][0]
                localmatch[i][1] = localmatch[i][1][0]
        ## drop all the non-orthologous exon that are needed violent alignment
        dropi = []
        for i in range(0,len(total)):
            if isinstance(total[i][0],list):
                for j in total[i][0]:
                    dropi.append(listwhich(selectcol(homogrp,"hsag"),"==",j)[0])
                for j in total[i][1]:
                    dropi.append(listwhich(selectcol(homogrp,"ptrg"),"==",j)[0])
            else:      
                dropi.append(listwhich(selectcol(homogrp,"hsag"),"==",total[i][0])[0])
                dropi.append(listwhich(selectcol(homogrp,"ptrg"),"==",total[i][1])[0])
        if len(dropi) > 0:
            homogrp["dat"] = listdrop(homogrp["dat"],dropi)
        ### homogrp might be empty after dropping, for example retro-gene
        if len(homogrp["dat"]) == 0:
            homogrp["dat"] = copy.deepcopy(localmatch)
            ### insert localmatch to orthologous group
        else:
            for i in range(0,len(localmatch)):
                if isinstance(localmatch[i][0],list):
                    if listasint(localmatch[i][0],"max")[0] > max(listasint(selectcol(homogrp,"hsag"),"max")):
                        homogrp["dat"].append(localmatch[i])
                    else:
                        inserti = min(listwhich(listasint(selectcol(homogrp,"hsag"),"min"),">",listasint(localmatch[i][0],"max")[0]))
                        homogrp["dat"].insert(inserti,localmatch[i])
                else:
                    if localmatch[i][0] > max(listasint(selectcol(homogrp,"hsag"),"max")):
                        homogrp["dat"].append(localmatch[i])
                    else:
                        inserti = min(listwhich(listasint(selectcol(homogrp,"hsag"),"min"),">",localmatch[i][0]))
                        homogrp["dat"].insert(inserti,localmatch[i])
        ##----------- organize orthologous exon group -------------##
        for i in range(0,len(homogrp["dat"])):
            if isinstance(homogrp["dat"][i][0],list) and len(homogrp["dat"][i][0]) == 1 and len(homogrp["dat"][i][1]) == 1:
                homogrp["dat"][i][0] = homogrp["dat"][i][0][0]
                homogrp["dat"][i][1] = homogrp["dat"][i][1][0]
    ##=====================================================
    #print("optimise ruler alignment based on homologous exon group")
    hsgroup = selectcol(regrhs,"group")
    ttgroup = selectcol(regrtt,"group")
    for hsg in hsgroup:
        if not(hsg in listasint(selectcol(homogrp,"hsag"),deal="all")):
            if hsg < max(listasint(selectcol(homogrp,"hsag"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"hsag"),deal="min"),">",hsg)),[hsg,0])
            else:
                homogrp["dat"].append([hsg,0])
    for ttg in ttgroup:
        if not(ttg in listasint(selectcol(homogrp,"ptrg"),deal="all")):
            if ttg < max(listasint(selectcol(homogrp,"ptrg"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"ptrg"),deal="min"),">",ttg)),[0,ttg])
            else:
                homogrp["dat"].append([0,ttg])
    ##=====================================================##
    ##====== test if lost group during optimise       =====##
    ##=====================================================##
    refhsa = selectcol(regrhs,"group")
    refptr = selectcol(regrtt,"group")
    tsthsa = []
    tstptr = []
    for i in range(0,len(homogrp["dat"])):
        if isinstance(homogrp["dat"][i][0],list):
            tsthsa = tsthsa + homogrp["dat"][i][0]
            tstptr = tstptr + homogrp["dat"][i][1]
        else:
            if homogrp["dat"][i][0] != 0:
                tsthsa.append(homogrp["dat"][i][0])
            if homogrp["dat"][i][1] != 0:
                tstptr.append(homogrp["dat"][i][1])  
    ##----------------------------------------
    diffhsa = setdiff(refhsa,tsthsa)
    diffptr = setdiff(refptr,tstptr)
    if len(diffhsa) != 0 or len(diffptr) != 0:
        raise SyntaxError("homo_colinearity_test: group loss after homo-corlinearity optimise")
    ##----------------------------------------
    if len(homogrp["dat"])>1:
        for i in range(1,len(homogrp["dat"])):
            if homogrp["dat"][i][0] != 0:
                test1 = min(listasint(homogrp["dat"][i][0])) > max(listasint(selectcol(homogrp,"hsag")[0:i],"max"))
            else:
                test1 = True
            if homogrp["dat"][i][1] != 0:
                test2 = min(listasint(homogrp["dat"][i][1])) > max(listasint(selectcol(homogrp,"ptrg")[0:i],"max"))
            else:
                test2 = True           
            if not(test1) or not(test2):
                raise SyntaxError("homo_colinearity_test: corlinearity test fail"+str()+":"+str())
    ##=====================================##
    ##============ return results =========##
    ##=====================================##
    return homogrp

def orgexongroup(genhs,regrhs,gentt,regrtt,refhomo):
    ## organize ruler group information
    store = []
    newgroup = 0
    sti = listwhich(regrhs["coln"],"==","start")[0]
    eni = listwhich(regrhs["coln"],"==","end")[0]
    seqi = listwhich(regrhs["coln"],"==","seq")[0]
    for i in range(0,len(refhomo["dat"])):
        newgroup = newgroup + 1
        if not(isinstance(refhomo["dat"][i][0],list)):
            hsg = refhomo["dat"][i][0]
            ttg = refhomo["dat"][i][1]
            if hsg == 0 and ttg != 0:
                hspos = "None"
                ttstart = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][sti]
                ttend = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][eni]
                ttpos = ":".join([str(ttstart),str(ttend)])
                iden = [0]
                grptype = "0-1"
            elif hsg != 0 and ttg == 0:
                hsstart = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][sti]
                hsend = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][eni]
                hspos = ":".join([str(hsstart),str(hsend)])
                ttpos = "None"
                iden = [0]
                grptype = "1-0" 
            else:
                hsstart = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][sti]
                hsend = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][eni]
                hspos = ":".join([str(hsstart),str(hsend)])
                ##----------------------------------------
                ttstart = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][sti]
                ttend = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][eni]
                ttpos = ":".join([str(ttstart),str(ttend)])
                ##----------------------------------------
                hsseqtmp = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][seqi]
                ttseqtmp = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][seqi]
                ##----------------------------------------
                iden = calIdentity(hsseqtmp, ttseqtmp, method="global")
                grptype = "1-1"
        else:
            hsgtmp = []
            hssttmp = []
            hsentmp = []
            hsseqtmp = ""
            for hsii in refhomo["dat"][i][0]:
                hsseqtmp = hsseqtmp + regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsii)[0]][seqi]
                hssttmp.append(regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsii)[0]][sti])
                hsentmp.append(regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsii)[0]][eni])
                hsgtmp.append(hsii)
            ##----------------------------------------
            ttgtmp = []
            ttsttmp = []
            ttentmp = []
            ttseqtmp = ""
            for ttii in refhomo["dat"][i][1]:
                ttseqtmp = ttseqtmp + regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttii)[0]][seqi]
                ttsttmp.append(regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttii)[0]][sti])
                ttentmp.append(regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttii)[0]][eni])                
                ttgtmp.append(ttii)
            ##----------------------------------------
            hspos = ":".join([str(min(hssttmp)),str(max(hsentmp))])
            ttpos = ":".join([str(min(ttsttmp)),str(max(ttentmp))])
            #### it is important to use global alignment
            #sumalignlo = pairwisealign(hsseqtmp,ttseqtmp,method="global")           
            iden = calIdentity(hsseqtmp, ttseqtmp, method="global")
            if len(hsgtmp) == 1 and len(ttgtmp) > 1:
                grptype = "1-N"
            elif len(hsgtmp) > 1 and len(ttgtmp) == 1:
                grptype = "N-1"
            else:
                grptype = "N-N"
        #matchtype = refhomo["dat"][i][2]
        store.append([newgroup,genhs,hspos,gentt,ttpos,iden[0],grptype])

    return {"dat":store, "coln":["group","hsagen","hsapos","ptrgen","ptrpos","identity","grouptype"]}

##======================================================##
##          detection of orthologous isoforms           ##
##======================================================##

## summarize cds or utr composition
def regComSum(genestd,regidxd,dattype="cds/utr"):
    genename = genestd["coln"]
    genest = genestd["dat"]
    regidxname = regidxd["coln"]
    regidx = regidxd["dat"]

    regnum = len(regidx)
    samnum = len(genest)
    coli = listwhich(genename,"==","ensemblt")[0]
    isonam = ["id","start"]
    for n in list(range(0,len(genest))):
        isonam.append(genest[n][coli])
    
    commatrix = []
    starti = listwhich(regidxname,"==","start")[0]
    endi = listwhich(regidxname,"==","end")[0]
    idi = listwhich(regidxname,"==","id")[0]
    typei = listwhich(genename,"==",dattype)[0]
    for i in list(range(0,regnum)):
        regposabb = ":".join([str(regidx[i][starti]),str(regidx[i][endi])]) # compare the string rather pos num
        regstatus = []
        regstatus.append(regidx[i][idi])                  # the first column stores cds region id
        regstatus.append(int(regidx[i][starti]))  
        for j in list(range(0,samnum)):
            regpostmp = genest[j][typei].split("|")
            regstatmp = 0
            locnum = len(regpostmp)
            for m in list(range(0,locnum)):
                if regpostmp[m] == regposabb:
                    regstatmp = 1
                    break
            regstatus.append(regstatmp)
        commatrix.append(regstatus)
    ## drop void cds region
    commatd = {"dat":commatrix,"coln":isonam}
    countcom = [0 for _ in range(0,len(commatrix))]
    for i in range(0,len(commatrix)):
        countcom[i] = sum(commatrix[i][1:])
    seli = listwhich(countcom,"!=",0)
    regsum = selectlsls(commatd,seli,"all")
    return regsum

## detect conserved isoform at isoform level
def consertiveIso(exoncon,exonttcon,cdshs,cdstt,homoexondat,identhres):
    ##========================================================##
    ##------------- set conservative threshold  --------------##
    ##========================================================##
    ##========================================================##
    ##-------------- estimate exon group list  ---------------##
    ##========================================================##
    sti = listwhich(cdshs["coln"],"==","start")[0]
    eni = listwhich(cdshs["coln"],"==","end")[0]
    grpi = listwhich(homoexondat["coln"],"==","group")[0]
    ##---------- hsa -------------
    hsaexongroup = []
    for hse in exoncon:
        st = int(cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hse)[0]][sti])
        en = int(cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hse)[0]][eni])
        grptmpi = intersect(listwhich(selectcol(homoexondat,"hsast"),"<=",st),listwhich(selectcol(homoexondat,"hsaen"),">=",en))[0]
        hsaexongroup.append(homoexondat["dat"][grptmpi][grpi])
    hsaegrp = unique(hsaexongroup)
    ##---------- hsa -------------
    ptrexongroup = []
    for pte in exonttcon:
        st = int(cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",pte)[0]][sti])
        en = int(cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",pte)[0]][eni])
        grptmpi = intersect(listwhich(selectcol(homoexondat,"ptrst"),"<=",st),listwhich(selectcol(homoexondat,"ptren"),">=",en))[0]
        ptrexongroup.append(homoexondat["dat"][grptmpi][grpi])
    ptregrp = unique(ptrexongroup)
    ##--------- compare ----------
    if hsaegrp != ptregrp:
        passtest = 0
    else:
        passtest = 1

    ##========================================================##
    ##----------- need further sequence test ??? -------------##
    ##========================================================##
    infostr = ""
    if passtest == 1:
        seqi = listwhich(cdshs["coln"],"==","seq")[0]
        idenlist = []
        infotmp = []
        for i in range(0,len(hsaegrp)):
            hsae = selectele(exoncon,listwhich(hsaexongroup,"==",hsaegrp[i]))
            ptre = selectele(exonttcon,listwhich(ptrexongroup,"==",hsaegrp[i]))

            hsaseq = ""
            ptrseq = ""
            for he in hsae:
                hsaseq = hsaseq + cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",he)[0]][seqi]
            for pe in ptre:
                ptrseq = ptrseq + cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",pe)[0]][seqi]

            identmp = calIdentity(hsaseq, ptrseq, method="global")[0]
            idenlist.append(identmp)
            infotmp.append(hsae[0]+"-"+ptre[0]+"|"+str(identmp))    
        
        if len(listwhich(idenlist,"<",identhres)) > 0:
            passtest = 0

        infostr = ";".join(infotmp)
    
    suminfo = {"passtest":passtest, "info":infostr}
    
    return suminfo

## detect conserved isoform at gene level
def orthoiso(species1, species2, genhs, genhsstr, strandhs, cdshs, gentt, genttstr, strandtt, cdstt, conexondat, identhres):
    ##========================================================##
    ##--------------- reorganize conexondat  ----------------##
    ##========================================================##
    homoexon = []
    grpi = listwhich(conexondat["coln"],"==","Group")[0]
    hsapi = listwhich(conexondat["coln"],"==",species1+"Pos")[0]
    ptrpi = listwhich(conexondat["coln"],"==",species2+"Pos")[0]
    ideni = listwhich(conexondat["coln"],"==","Iden")[0]
    for i in range(0,len(conexondat["dat"])):
        ##---------- hsa -------------
        if conexondat["dat"][i][hsapi] == "None":
            hsast = 0
            hsaen = 0
        else:
            hsagrouppos = conexondat["dat"][i][hsapi].split(":")
            hsast = int(hsagrouppos[0])
            hsaen = int(hsagrouppos[1])
        ##---------- ptr -------------
        if conexondat["dat"][i][ptrpi] == "None":
            ptrst = 0
            ptren = 0
        else:
            ptrgrouppos = conexondat["dat"][i][ptrpi].split(":")
            ptrst = int(ptrgrouppos[0])
            ptren = int(ptrgrouppos[1])
        ##---------- organize ------------
        homoexon.append([conexondat["dat"][i][grpi],hsast,hsaen,ptrst,ptren,conexondat["dat"][i][ideni]])
        
    homoexondat = {"dat":homoexon,"coln":["group","hsast","hsaen","ptrst","ptren","iden"]}
    ##========================================================##
    ##-------------     estimate exon group     --------------##
    ##========================================================##
    sti = listwhich(cdshs["coln"],"==","start")[0]
    eni = listwhich(cdshs["coln"],"==","end")[0]   
    ##---------- hsa -------------
    for i in range(0,len(cdshs["dat"])):
        grptmpi = intersect(listwhich(selectcol(homoexondat,"hsast"),"<=",int(cdshs["dat"][i][sti])),listwhich(selectcol(homoexondat,"hsaen"),">=",int(cdshs["dat"][i][eni])))[0]
        cdshs["dat"][i].insert(0,homoexondat["dat"][grptmpi][grpi])
    cdshs["coln"].insert(0,"group")
    ##---------- ptr -------------
    for i in range(0,len(cdstt["dat"])):
        grptmpi = intersect(listwhich(selectcol(homoexondat,"ptrst"),"<=",int(cdstt["dat"][i][sti])),listwhich(selectcol(homoexondat,"ptren"),">=",int(cdstt["dat"][i][eni])))[0]
        cdstt["dat"][i].insert(0,homoexondat["dat"][grptmpi][grpi])
    cdstt["coln"].insert(0,"group")

    ##========================================================##
    ##---------- estimate isoform CDS composition  -----------##
    ##========================================================##
    geninfohs = aslsls(genhsstr, datcolumns=["ensemblt","orf","cds","aaseq"])
    geninfott = aslsls(genttstr, datcolumns=["ensemblt","orf","cds","aaseq"])
    ##============ estimate cds ruler and gap =================
    cdsmtt = regComSum(genestd=geninfott,regidxd=cdstt,dattype="cds")
    if strandtt == 1:
        isott = listsort(cdsmtt,by=["start"],increase = True)
    if strandtt == -1:
        isott = listsort(cdsmtt,by=["start"],increase = False)
    isott = selectlsls(isott,list(range(0,len(isott["dat"]))),selectele(isott["coln"],[0]+list(range(2,len(isott["coln"])))))
    ##--------------------------------
    cdsmhs = regComSum(genestd=geninfohs,regidxd=cdshs,dattype="cds")
    if strandhs == 1:
        isohs = listsort(cdsmhs,by=["start"],increase = True)
    if strandhs == -1:
        isohs = listsort(cdsmhs,by=["start"],increase = False)
    isohs = selectlsls(isohs,list(range(0,len(isohs["dat"]))),selectele(isohs["coln"],[0]+list(range(2,len(isohs["coln"])))))
    
    ##========================================================##
    ##------------ estimate homologous isoform  --------------##
    ##========================================================##
    isoformhs = selectele(isohs["coln"],list(range(1,len(isohs["coln"]))))
    conseriso = []
    conhsgen = []
    conttgen = []
    isohsa = []
    exoniden = []
    for j in range(1,len(isohs["coln"])):
        #print("analysis " + isoformhs[j-1] + ":", end="")
        exoncon = selectele(selectcol(isohs,"id"),listwhich(selectcol(isohs,isohs["coln"][j]),"==",1))
        tariso = "."
        for e in range(1,len(isott["coln"])):
            #print(".",end="")
            exonttcon = selectele(selectcol(isott,"id"),listwhich(selectcol(isott,isott["coln"][e]),"==",1))
            ## test consertive
            conser = consertiveIso(exoncon,exonttcon,cdshs,cdstt,homoexondat,identhres=identhres)
            if conser["passtest"] == 1:
                ##================================================##
                ### test ORF shift using R language!!!!!!!!!!!!!! ##
                ##================================================##
                tariso = isott["coln"][e]
                isohsa.append(isoformhs[j-1])
                conseriso.append(isott["coln"][e])
                conhsgen.append(genhs)
                conttgen.append(gentt)
                exoniden.append(conser["info"])
        
        #if tariso == ".":
            #print(" fail")
        #else:
            #print(" succeed")

    summ = {"dat":transpose([conhsgen,conttgen,isohsa,conseriso,exoniden]),"coln":[species1,species2,species1+"iso",species2+"iso","exonmatch"]}
    return summ

##======================================================##
##                   analysis process                   ##
##======================================================##
##****************************************************************************************************************************##
#orthog = pd.read_table('/Users/jeffma/MJFPhD/HsSpDB/step6_homologs/2homologs_gene_blastp_inparanoid8.dat',header=0,sep='\t')
#blastexon = pd.read_table('/Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/blastn_mapping.tab',header=0,sep='\t')

#orthogene = pd.read_table('/Users/majinfa/MJFPhD/HsSpDB/step6_homologs/ortho_hsa_mmu.txt',header=0,sep='\t')
#blastexon = pd.read_table('/Users/majinfa/MJFPhD/HsSpDB/step7_orthologousExon/2homologous_exons_whole_one_to_one_withinhomogen_hsa_mmu.dat',header=0,sep='\t')

#isocom1 = pd.read_table('/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0hsa/mysql_traninfo.tab',header=0,sep='\t')
#cdscom1 = pd.read_table('/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0hsa/mysql_exon_seqinfo.tab',header=0,sep='\t')

#isocom2 = pd.read_table('/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0ptr/mysql_traninfo.tab',header=0,sep='\t')
#cdscom2 = pd.read_table('/Users/jeffma/MJFPhD/HsSpDB/step5_db_for_Diamond/0ptr/mysql_exon_seqinfo.tab',header=0,sep='\t')

#species1 = "hsa"
#species2 = "ptr"

#identhreshold = 0.8
#coverthreshold = 0.8

#mapscore = 2
#misscore = -2
#gapscore = -1
#minexon = 2

def EGIO(data_dict):

    genhs = data_dict["raw_data"][0].iloc[data_dict["order"]][str(data_dict["raw_data"][12])]

    print("analysis of gene " + genhs + " ")

    strandhsstr = unique(list(data_dict["raw_data"][2].iloc[listwhich(list(data_dict["raw_data"][2]["EnsemblG"]),"==",genhs)]["Strand"]))[0]
    if strandhsstr == "-":
        strandhs = -1
    else:
        strandhs = 1
        
    genhscdstmp1 = data_dict["raw_data"][2].iloc[listwhich(list(data_dict["raw_data"][2]["EnsemblG"]),"==",genhs)][["EnsemblT","Orf","Exoncom"]]
    genhscdsdat = genhscdstmp1.iloc[listwhich(list(genhscdstmp1["Orf"]),"!=",".|.")]
        
    cdshscds = data_dict["raw_data"][3].iloc[listwhich(list(data_dict["raw_data"][3]["EnsemblG"]),"==",genhs)][["ID","Start","End","Seq"]]
    genhscds = translatecdna(genhscdsdat,cdshscds,strandhsstr)

    genhsstr = transform(data=genhscds)
    cdshsstr = transform(data=cdshscds)
    #####==============================================#####
    #####=================== chimpanzee ===============#####
    #####==============================================#####
    gentt = data_dict["raw_data"][0].iloc[data_dict["order"]][str(data_dict["raw_data"][13])]
    strandttstr = unique(list(data_dict["raw_data"][4].iloc[listwhich(list(data_dict["raw_data"][4]["EnsemblG"]),"==",gentt)]["Strand"]))[0]
    if strandttstr == "-":
        strandtt = -1
    else:
        strandtt = 1

    genttcdstmp1 = data_dict["raw_data"][4].iloc[listwhich(list(data_dict["raw_data"][4]["EnsemblG"]),"==",gentt)][["EnsemblT","Orf","Exoncom"]]
    genttcdsdat = genttcdstmp1.iloc[listwhich(list(genttcdstmp1["Orf"]),"!=",".|.")]
    
    cdsttcds = data_dict["raw_data"][5].iloc[listwhich(list(data_dict["raw_data"][5]["EnsemblG"]),"==",gentt)][["ID","Start","End","Seq"]]
    genttcds = translatecdna(genttcdsdat,cdsttcds,strandttstr)
        
    ##-----------------------------------
    genttstr = transform(data=genttcds)
    cdsttstr = transform(data=cdsttcds)
    #####==============================================#####
    #####=================== blastexon  ===============#####
    #####==============================================#####
    if len(data_dict["raw_data"][1]) > 0:
        selecthspti = intersect(listwhich(list(data_dict["raw_data"][1][str(data_dict["raw_data"][12])]),"==",genhs),listwhich(list(data_dict["raw_data"][1][str(data_dict["raw_data"][13])]),"==",gentt))
        homoexontmp = data_dict["raw_data"][1].iloc[selecthspti][[str(data_dict["raw_data"][12]),str(data_dict["raw_data"][13])]]
        ## transform into my basic function fitting structure
        homoexon = []
        if len(homoexontmp) > 0:
            for hei in range(0,len(homoexontmp)):
                homoexon.append([homoexontmp.iloc[hei][str(data_dict["raw_data"][12])],homoexontmp.iloc[hei][str(data_dict["raw_data"][13])]])
        homocdsdat = {"dat":homoexon, "coln":[str(data_dict["raw_data"][12]),str(data_dict["raw_data"][13])]}
    else:
        homocdsdat = {"dat":[], "coln":[str(data_dict["raw_data"][12]),str(data_dict["raw_data"][13])]}

    #####==============================================#####
    #####===================  analysis  ===============#####
    #####==============================================#####
    cdshs = aslsls(cdshsstr, datcolumns=["id","start","end","seq"])
    cdstt = aslsls(cdsttstr, datcolumns=["id","start","end","seq"])
    ##------------------ estimate cds ruler ------------------
    #print("estimate unique exon region:")
    cdsrulerhs = regRuler(regidxd=cdshs,genstr=strandhs)
    cdsrulertt = regRuler(regidxd=cdstt,genstr=strandtt)
    ##------------ estimate CDS homologous group  ------------
    #print("organize homologous exon group:")
    #regrhs=cdsrulerhs;regrtt=cdsrulertt
    refhomo = homo_colinearity_test(data_dict["raw_data"][12], data_dict["raw_data"][13], cdsrulerhs,cdsrulertt,cdshs,cdstt,homocdsdat, identhres=data_dict["raw_data"][9], coverthres=data_dict["raw_data"][10], minexon=data_dict["raw_data"][11], mapscore=data_dict["raw_data"][6], misscore=data_dict["raw_data"][7], gapscore=data_dict["raw_data"][8])
        
    ##=======================================================================##
    ##                         website back end scripts                      ##
    ##=======================================================================##
    #print("summarize exon groups")
    ## genhs=genhs;regrhs=cdsrulerhs;gentt=gentt;regrtt=cdsrulertt;refhomo=refhomo
    sumexon = orgexongroup(genhs,cdsrulerhs,gentt,cdsrulertt,refhomo)
    
    ##=======================================================================##
    orthoexon = pd.DataFrame(sumexon["dat"])
    orthoexon.columns = ["Group",str(data_dict["raw_data"][12])+"EnsemblG",str(data_dict["raw_data"][12])+"Pos",str(data_dict["raw_data"][13])+"EnsemblG",str(data_dict["raw_data"][13])+"Pos","Iden","Type"]
    
    #####==============================================#####
    #####============= read blastexon data ============#####
    #####==============================================#####
    #print("Detect orthologous isoforms")
    homoexon = []
    if len(orthoexon) > 0:
        for hei in range(0,len(orthoexon)):
            homoexon.append([orthoexon.iloc[hei]["Group"],orthoexon.iloc[hei][str(data_dict["raw_data"][12])+"Pos"],orthoexon.iloc[hei][str(data_dict["raw_data"][13])+"Pos"],orthoexon.iloc[hei]["Iden"],orthoexon.iloc[hei]["Type"]])
    conexondat = {"dat":homoexon, "coln":["Group",str(data_dict["raw_data"][12])+"Pos",str(data_dict["raw_data"][13])+"Pos","Iden","Type"]}
    #####==============================================#####
    #####============= detect homo isoofrm ============#####
    #####==============================================#####
    sumiso = orthoiso(data_dict["raw_data"][12], data_dict["raw_data"][13],genhs, genhsstr, strandhs, cdshs, gentt, genttstr, strandtt, cdstt, conexondat,identhres=data_dict["raw_data"][9])
    
    #sumiso = orthoiso(str(species1),str(species2),genhs, genhsstr, strandhs, cdshs, gentt, genttstr, strandtt, cdstt, conexondat,identhres=identhreshold)
        
    #print("------------------------------------")
    #####==============================================#####
    #####============= merge all the data =============#####
    #####==============================================#####

    return {"oe":sumexon, "ot":sumiso}

##****************************************************************************************************************************##

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'detection of orthologous exons and isoforms using EGIO')
    parser.add_argument('--orthog',
        nargs = '?',
        help = 'path of file contains orthologous gene pairs')
    parser.add_argument('--blastn',
        nargs = '?',
        help = 'path of file contains homologous exon pairs')
    parser.add_argument('--species1',
        nargs = '?',
        help = 'query species')
    parser.add_argument('--species2',
        nargs = '?',
        help = 'subject species')
    parser.add_argument('--isocom1',
        nargs = '?',
        default = None,
        help = 'path of file contains isoform information of species1')
    parser.add_argument('--cdscom1',
        nargs = '?',
        default = None,
        help = 'path of file contains cds information of species1')
    parser.add_argument('--isocom2',
        nargs = '?',
        default = None,
        help = 'path of file contains isoform information of species2')
    parser.add_argument('--cdscom2',
        nargs='?',
        default = None,
        help = 'path of file contains cds information of species2')
    parser.add_argument('--identhres',
        nargs = '?',
        default = 0.8,
        type = float,
        help = 'identity threshold of pairwise alignment')
    parser.add_argument('--coverthres',
        nargs = '?',
        default = 0.8,
        type = float,
        help = 'recoprocal coverage threshold of pairwise alignment')
    parser.add_argument('--match',
        nargs = '?',
        default = 2,
        type = int,
        help = 'match score of dynamic programming')
    parser.add_argument('--mismatch',
        nargs = '?',
        default = -2,
        type = int,
        help = 'mismatch score of dynamic programming')
    parser.add_argument('--gap',
        nargs = '?',
        default = -1,
        type = int,
        help = 'gapscore score of dynamic programming')
    parser.add_argument('--minexon',
        nargs = '?',
        default = 2,
        type = int,
        help = 'exons less than minexon size to be optimised after dynamic programming')

    parser.add_argument('--pnum',
        nargs = '?',
        default = 1,
        type = int,
        help = 'cpu core numbers to run EGIO')
    
    args = parser.parse_args()

    #outpath = os.getcwd()

    orthogene = pd.read_table(str(args.orthog),header=0,sep='\t')
    blastexon = pd.read_table(str(args.blastn),header=0,sep='\t')

    isocom1 = pd.read_table(str(args.isocom1),header=0,sep='\t')
    cdscom1 = pd.read_table(str(args.cdscom1),header=0,sep='\t')

    isocom2 = pd.read_table(str(args.isocom2),header=0,sep='\t')
    cdscom2 = pd.read_table(str(args.cdscom2),header=0,sep='\t')

    ##------------------------------

    data_in = (orthogene, blastexon, isocom1, cdscom1, isocom2, cdscom2, args.match, args.mismatch, args.gap, args.identhres, args.coverthres,  args.minexon, args.species1, args.species2)
    
    p = Pool(args.pnum)
    res_l = []
    for gi in range(0,len(orthogene)):
        res = p.apply_async(EGIO,args=({"raw_data":data_in,"order":gi},))
        res_l.append(res)

    p.close()
    p.join()

    exoniden = [["Group",'hsaEnsemblG',"hsaPos","ptrEnsemblG",'ptrPos','Iden',"Type"]]
    isoiden = [["hsa","ptr","isohsa","isoptr","exoniden"]]
    for res in res_l:
        datatmp = res.get()
        exoniden = exoniden + datatmp["oe"]["dat"]
        isoiden = isoiden + datatmp["ot"]["dat"]

    exonidendf = pd.DataFrame(exoniden)
    isoidendf = pd.DataFrame(isoiden)

    exonidendf.to_csv(os.getcwd()+'/ExonGroup_testpro.txt', sep='\t', header=False, index=False)
    isoidendf.to_csv(os.getcwd()+'/OrthoIso_testpro.txt', sep='\t', header=False, index=False)


##****************************************************************************************************************************##




    #runEGIO(args.orthog, args.blastn, args.species1, args.species2, args.isocom1, args.cdscom1, args.isocom2, args.cdscom2, args.identhres,args.coverthres, args.match, args.mismatch, args.gap, args.minexon, args.pnum)

# cd /Users/jeffma/MJFPhD/HsSpDB/step0_GEIO_script/test/
# python __EGIO.py --orthog homolog_gene.dat --blastn blastn_hsa_ptr.dat --species1 hsa --species2 ptr --isocom1 mysql_traninfo_hsa.tab --cdscom1 mysql_exon_hsa.tab --isocom2 mysql_traninfo_ptr.tab --cdscom2 mysql_exon_ptr.tab

#python __EGIO.py --orthog homolog_gene.dat --blastn blastn_mapping.tab --species1 hsa --species2 ptr --isocom1 hsa.tran --cdscom1 hsa.exon --isocom2 ptr.tran --cdscom2 ptr.exon --identhres 0.8 --coverthres 0.8 --match 2 --mismatch -2 --gap -1
