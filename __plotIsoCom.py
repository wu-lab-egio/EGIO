import argparse
import os
import pandas as pd
import sys
import string
import ctypes
import itertools
import copy
import math
from plotnine import *


##======================================================##
align= ctypes.cdll.LoadLibrary(os.getcwd() + '/___pairwisealign.so')
align.pairwisealign.argtypes = [ctypes.c_char_p]
align.pairwisealign.restype = ctypes.c_char_p

##======================================================##
## correct alignment based on cds ruler
def pairwisealign(seqref,seqtst):
    if len(seqref) == 0 and len(seqtst) > 0:
        correctref = "".join(["-" for _ in range(0,len(seqtst))])
        correcttst = seqtst
    elif len(seqref) > 0 and len(seqtst) == 0:
        correctref = seqref
        correcttst = "".join(["-" for _ in range(0,len(seqref))])
    elif len(seqref) > 0 and len(seqtst) > 0:
        hsseqlo = ctypes.c_char_p(bytes(seqref, 'utf-8'))
        ttseqlo = ctypes.c_char_p(bytes(seqtst, 'utf-8'))
        sumalignlotmp = align.pairwisealign(hsseqlo,ttseqlo,"local",-10)
        #sumalignlotmp = alignmentg.blastWhole(hsseqlo,ttseqlo)
        sumalignlo = sumalignlotmp.decode().split("|")
        lohs = sumalignlo[0].replace("B","-")
        lott = sumalignlo[1].replace("B","-")
        ##====== find excluded sequence ============
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
        ##===================================
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
    ##====== find excluded sequence ============
    ##===================================
    listsum = {"corref":correctref,"cortst":correcttst}
    return listsum

## def which as in R
def which(judge):
    backdata = [x for x in range(0,len(judge)) if judge[x]]
    return backdata

## def which as in R
def listwhich(lsdata,deal,value):
    if deal == "==":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] == value]
    if deal == "!=":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] != value]
    if deal == ">":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] > value]
    if deal == ">=":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] >= value]
    if deal == "<":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] < value]
    if deal == "<=":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] <= value]         
    return(judge)
    
## def intersect as in R
def intersect(ls1, ls2):
    store = [x for x in ls1 if x in ls2]
    return store

## type conversion of a certain column
def typeinto(listd,colname,astype=int):
    listdat = listd["dat"]
    listname = listd["coln"]
    coli = listwhich(listname,"==",colname)[0]
    for i in range(0,len(listdat)):
        listdat[i][coli] = astype(listdat[i][coli])
    info = {"dat":listdat,"coln":listname}
    return info

## define unique
def unique(listdat):
    store = []
    for i in range(0,len(listdat)):
        if not(listdat[i] in store):
            store.append(listdat[i])
    store.sort(reverse = False)
    return store

## define get a column
def selectcol(listd,target):
    seli = listd["coln"].index(target)
    if not(target in listd["coln"]):
        raise SyntaxError("selectcol: wrong colname is given, please check")
    store =[listd["dat"][x][seli] for x in range(0,len(listd["dat"]))]
    return store

## define select element
def selectele(lsdata,targeti):
    store = [lsdata[x] for x in targeti]
    return store

## list calculation
def listcal(lsdata,calcu,value):
    if calcu == "+":
        listtmp = [x+value for x in lsdata]
    if calcu == "-":
        listtmp = [x-value for x in lsdata]
    if calcu == "*":
        listtmp = [x*value for x in lsdata]
    if calcu == "/":
        if value == 0:
            raise SyntaxError("listcal: 0 can not be divided")
        listtmp = [x/value for x in lsdata]
    if calcu == "%":
        if value == 0:
            raise SyntaxError("listcal: 0 can not be divided")
        listtmp = [x%value for x in lsdata]
    if calcu == "^":
        listtmp = [x^value for x in lsdata]
    return listtmp

## list list calculation
def listlistcal(ls1,calcu,ls2):
    if len(ls1) != len(ls2):
        raise SyntaxError("listlistcal: the two lists do not have the same length")
    
    if calcu == "+":
        listtmp = [ls1[x]+ls2[x] for x in range(0,len(ls1))]
    if calcu == "-":
        listtmp = [ls1[x]-ls2[x] for x in range(0,len(ls1))]
    if calcu == "*":
        listtmp = [ls1[x]*ls2[x] for x in range(0,len(ls1))]
    if calcu == "/":
        if 0 in ls2:
            raise SyntaxError("listlistcal: 0 can not be divided")
        listtmp = [ls1[x]/ls2[x] for x in range(0,len(ls1))]
    if calcu == "%":
        if 0 in ls2:
            raise SyntaxError("listlistcal: 0 can not be divided")
        listtmp = [ls1[x]%ls2[x] for x in range(0,len(ls1))]
    if calcu == "^":
        listtmp = [ls1[x]^ls2[x] for x in range(0,len(ls1))]
    return listtmp

## list drop based row index [one dimension]
def listdrop(lsdata,dropi):
    store = [lsdata[x] for x in range(0,len(lsdata)) if not(x in dropi)]
    return store

## merge two list by column binding
def transpose(lsls):
    new_list = [[lsls[j][i] for j in range(0,len(lsls))] for i in range(0,len(lsls[0]))]
    return new_list

## sort list list
def listsort(listd,by,increase):
    if len(by) > 2:
        raise SyntaxError("the max element to sort is 2")
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

## def list-list merge
def listmerge(listd1, listd2, by):
    list1 = listd1["dat"]
    lsname1 = listd1["coln"]
    list2 = listd2["dat"]
    lsname2 = listd2["coln"]
    dropi = listwhich(lsname2,"==",by)
    reflist = selectcol(listd1,by)
    tstlist = selectcol(listd2,by)
    sorti = [0 for _ in range(0,len(reflist))]
    for i in range(0,len(reflist)):
        refid = reflist[i]
        sorti[i] = listwhich(tstlist,"==",refid)[0]
    newtst = []
    for i in range(0,len(sorti)):
        newtst.append(list2[sorti[i]])
    list1_t = transpose(list1)
    list2_t = transpose(newtst)
    list2_tt = listdrop(list2_t,dropi)
    lsname2t = listdrop(lsname2,dropi)
    store = list1_t + list2_tt
    listff = transpose(store)
    listffname = lsname1 + lsname2t
    listsum = {"dat":listff,"coln":listffname}
    return listsum

## select listlist
def selectlsls(listd,rowi,coln="all"):
    if coln == "all":
        coln = listd["coln"]
    coli = [listd["coln"].index(x) for x in coln]
    store = [[listd["dat"][y][x] for x in coli] for y in rowi]
    listsum = {"dat":store,"coln":coln}
    return listsum

## drop duplicatecd element
def dropduplicated(listdat):
    if len(listdat) == 0:
        unilist = []
    if len(listdat) > 0:
        unilist = [listdat[0]]
        for i in range(0,len(listdat)):
            if not(listdat[i] in unilist):
                unilist.append(listdat[i])
    return unilist

## calculate a column and return a lsls
def listColCal(listd,targetcol,calcu,value):
    seli = listwhich(listd["coln"],"==",targetcol)[0]
    listdat = listd["dat"]
    listtmp = copy.deepcopy(listdat)
    for i in range(0,len(listdat)):
        if calcu == "+":
            listtmp[i][seli] = float(listdat[i][seli]) + value
        if calcu == "-":
            listtmp[i][seli] = float(listdat[i][seli]) - value
        if calcu == "*":
            listtmp[i][seli] = float(listdat[i][seli]) * value
        if calcu == "/":
            listtmp[i][seli] = float(listdat[i][seli]) / value
        if calcu == "%":
            listtmp[i][seli] = float(listdat[i][seli]) % value
        if calcu == "^":
            listtmp[i][seli] = float(listdat[i][seli]) ** value
    listsum = {"dat":listtmp,"coln":listd["coln"]}
    return listsum

##======================================================##
##                   advanced function                  ##
##======================================================##

def testquery(speciespair, orthologpair, wdpath):
    #wdpath = os.getcwd()
    if wdpath[-1] != "/":
        wdpath = wdpath + "/"

    inpath = wdpath + "extrainfo/"
    ##===============================================
    infiles = os.listdir(inpath)
    wdfiles = os.listdir(wdpath)
    ##===============================================
    species1 = speciespair.split("-")[0]
    species2 = speciespair.split("-")[1]

    gene1 = orthologpair.split("-")[0]
    gene2 = orthologpair.split("-")[1]
    ##===============================================
    tran1 = species1+".tran"
    tran2 = species2+".tran"

    exon1 = species1+".exon"
    exon2 = species2+".exon"

    egfile = "ExonGroup_testpro_"+species1+"_"+species2+".txt"

    if not(tran1 in infiles):
        raise SyntaxError(tran1 + " is not found")
    if not(tran2 in infiles):
        raise SyntaxError(tran2 + " is not found")
    if not(exon1 in infiles):
        raise SyntaxError(exon1 + " is not found")
    if not(exon2 in infiles):
        raise SyntaxError(exon2 + " is not found")
    if not(egfile in wdfiles):
        raise SyntaxError(egfile + " is not found")
    ##===============================================
    egfile = pd.read_table(wdpath+egfile,header=0, sep="\t")

    genelist1 = list(egfile[species1+"EnsemblG"])
    genelist2 = list(egfile[species2+"EnsemblG"])

    genepair = []
    for i in range(0,len(egfile)):
        genepair.append("-".join([genelist1[i],genelist2[i]]))
    
    #if not(orthologpair in genepair):
    #    raise SyntaxError(orthologpair + " is not found")
    ##===============================================
    tranfile1 = pd.read_table(inpath+tran1,header=0, sep="\t")
    tranfile2 = pd.read_table(inpath+tran2,header=0, sep="\t")

    exonfile1 = pd.read_table(inpath+exon1,header=0, sep="\t")
    exonfile2 = pd.read_table(inpath+exon2,header=0, sep="\t")    
    ##===============================================
    tranidx1 = tranfile1[tranfile1.EnsemblG == gene1].index.tolist()
    if len(tranidx1) == 0:
        raise SyntaxError(gene1 + " is not recorded in the transcriptome")

    tranidx2 = tranfile2[tranfile2.EnsemblG == gene2].index.tolist()
    if len(tranidx2) == 0:
        raise SyntaxError(gene2 + " is not recorded in the transcriptome")

    exonidx1 = exonfile1[exonfile1.EnsemblG == gene1].index.tolist()
    exonidx2 = exonfile2[exonfile2.EnsemblG == gene2].index.tolist()
    egidx = listwhich(genepair,"==",orthologpair)
    ##===============================================
    genestr1tmp = []
    for i in tranidx1:
        if tranfile1.iloc[i]["Orf"] != ".|." :
            genestr1tmp.append("#".join([tranfile1.iloc[i]["EnsemblT"],
                    tranfile1.iloc[i]["Orf"],
                    tranfile1.iloc[i]["Strand"],
                    tranfile1.iloc[i]["Exoncom"]]))
    genestr1 = "@".join(genestr1tmp)

    genestr2tmp = []
    for i in tranidx2:
        if tranfile2.iloc[i]["Orf"] != ".|." :
            genestr2tmp.append("#".join([tranfile2.iloc[i]["EnsemblT"],
                    tranfile2.iloc[i]["Orf"],
                    tranfile2.iloc[i]["Strand"],
                    tranfile2.iloc[i]["Exoncom"]]))
    genestr2 = "@".join(genestr2tmp)
    ##------------------------------------------------
    cdsstr1tmp = []
    for i in exonidx1:
        cdsstr1tmp.append("#".join([exonfile1.iloc[i]["ID"],
                    str(exonfile1.iloc[i]["Start"]),
                    str(exonfile1.iloc[i]["End"]),
                    exonfile1.iloc[i]["Seq"]]))
    cdsstr1 = "@".join(cdsstr1tmp)

    cdsstr2tmp = []
    for i in exonidx2:
        cdsstr2tmp.append("#".join([exonfile2.iloc[i]["ID"],
                    str(exonfile2.iloc[i]["Start"]),
                    str(exonfile2.iloc[i]["End"]),
                    exonfile2.iloc[i]["Seq"]]))
    cdsstr2 = "@".join(cdsstr2tmp)
    ##------------------------------------------------
    egstrtmp = []
    for i in egidx:
        egstrtmp.append("#".join([str(egfile.iloc[i]["Group"]),
                    str(egfile.iloc[i][species1+"Pos"]),
                    str(egfile.iloc[i][species2+"Pos"]),
                    str(egfile.iloc[i]["Iden"])]))
    egstr = "@".join(egstrtmp)
    ##===============================================
    store = {"genestr1":genestr1, "genestr2":genestr2, "cdsstr1":cdsstr1, "cdsstr2":cdsstr2, "egstr":egstr}
    
    return store
    ##****************************************************************
    ##****************************************************************
    ##****************************************************************

def extract_CDS_from_Exon(orfpos,exoncom,strtmp):
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
    cdsposstr = [":".join([str(CDS[x][0]),str(CDS[x][1])]) for x in range(0,len(CDS))]

    return "|".join(cdsposstr)

## AsDataFrame: 
def aslsls(datastr,typename="isoinfo/codon/cds/motif/exong"):
    datastr = datastr.replace("^","|") # replace "^" to ";", which ";" is not allowed during php calling py
    rowdata = datastr.split("@")       # seperate row information, which is seperated by "@"
    rowlen = len(rowdata)              # seperate col information, which is seperated by "#"
    datainfo = []                                 # set a list to store seperated information
    for i in list(range(0,rowlen)):               # traversal to seperate each row to seperated columns
        datainfo.append(rowdata[i].split("#"))    # seperate row to seperated columns
    #dfdata.columns.values.tolist()
    if typename == "isoinfo":
        datcolumns=["ensemblt","orf","strand","cds"] # rename columns
        orfi = datcolumns.index("orf")
        cdsi = datcolumns.index("cds")
        stri = datcolumns.index("strand")
        for i in range(0,len(datainfo)):
            datainfo[i][cdsi] = extract_CDS_from_Exon(datainfo[i][orfi],datainfo[i][cdsi] ,datainfo[i][stri] )

    elif typename == "cds":
        datcolumns=["id","start","end","seq"]    # rename columns     
    elif typename == "motif":
        datcolumns=["type","uniid","leftaa","leftp","rigtaa","rigtp","annoid","ensemblt"]
    elif typename == "codon":
        datcolumns=["codon","abbname","sinname"]
    elif typename == "exong":
        datcolumns=["exong","hsapos","tttpos","iden"]
    
    datsum = {"dat":datainfo,"coln":datcolumns}
    return datsum

## summarize cds or utr composition
def regComSum(genestd,regidxd):
    genename = genestd["coln"]
    genest = genestd["dat"]
    regidxname = regidxd["coln"]
    regidx = regidxd["dat"]

    regnum = len(regidx)
    samnum = len(genest)
    coli = listwhich(genename,"==","ensemblt")[0]
    isonam = ["id"]
    for n in list(range(0,len(genest))):
        isonam.append(genest[n][coli])
    
    commatrix = []
    starti = listwhich(regidxname,"==","start")[0]
    endi = listwhich(regidxname,"==","end")[0]
    idi = listwhich(regidxname,"==","id")[0]
    typei = listwhich(genename,"==","cds")[0]
    for i in list(range(0,regnum)):
        gstart = int(regidx[i][starti])
        gend = int(regidx[i][endi])
        regposabb = ":".join([str(regidx[i][starti]),str(regidx[i][endi])]) # compare the string rather pos num
        regstatus = []
        regstatus.append(regidx[i][idi])               # the first column stores cds region id
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

## drop regions that mapping to no isoform
def newinfo(comd,regidxd):
    cdsregname = selectcol(comd,"id")
    seli = []
    idi = listwhich(regidxd["coln"],"==","id")[0]
    for i in range(0,len(regidxd["dat"])):
        if regidxd["dat"][i][idi] in cdsregname:
            seli.append(i)
    newreg = selectlsls(regidxd,seli,"all")
    return newreg

## merge overlapped sequence
def mergeOverlap(datad,grpstr):
    if grpstr == "+":
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
    if grpstr == "-":
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

## get region ruler
def regRuler(regidxd,genstr):
    regidxd = typeinto(regidxd,"start",int)
    regidxd = typeinto(regidxd,"end",int)
    regname = regidxd["coln"]
    if genstr == "+":
        regcomd = listsort(regidxd,by=["start","end"],increase = True)
    if genstr == "-":
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
        gsum.append([gid[i],gst,gen])
    gsumname = ["group","start","end"]
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

    ruler_t = gsum_t + posinfo_t + addinfo_t
    ruler = transpose(ruler_t)
    rulername = gsumd["coln"] + posinfod["coln"] + addinfod["coln"]
    rulerd = {"dat":ruler,"coln":rulername}
    return rulerd

##----------------------------------------
def regAlignTwo(regrhs,regrtt,exong):
    egi = exong["coln"].index("exong")
    hsaegposi = exong["coln"].index("hsapos")
    tttegposi = exong["coln"].index("tttpos")
    ideni = exong["coln"].index("iden")
    
    refseq = ""
    tstref = ""
    for i in range(0,len(exong["dat"])):
        if exong["dat"][i][hsaegposi] != "None":
            idx1 = listwhich(selectcol(regrhs,"start"),"<",int(exong["dat"][i][hsaegposi].split(":")[2]))
            idx2 = listwhich(selectcol(regrhs,"end"),">",int(exong["dat"][i][hsaegposi].split(":")[1]))
            idx = intersect(idx1, idx2)
            hsaseq = ""
            for j in idx:
                hsaseq = hsaseq + selectcol(regrhs,"seq")[j]
        else:
            hsaseq = ""
        
        if exong["dat"][i][tttegposi] != "None":
            idx1 = listwhich(selectcol(regrtt,"start"),"<",int(exong["dat"][i][tttegposi].split(":")[2]))
            idx2 = listwhich(selectcol(regrtt,"end"),">",int(exong["dat"][i][tttegposi].split(":")[1]))
            idx = intersect(idx1, idx2)
            tttseq = ""
            for j in idx:
                tttseq = tttseq + selectcol(regrtt,"seq")[j]
        else:
            tttseq = ""
        
        alignment = pairwisealign(hsaseq,tttseq)
        refseq = refseq + alignment["corref"]
        tstref = tstref + alignment["cortst"]

    listsum = {"corref":refseq,"cortst":tstref}
    return listsum

## get gap data, a common function to estimate gap data
def estiGap(rulermatrix,method):
    if method == "one":
        delinfo = {"dat":[],"coln":["width","start","end","alistart","aliend"]}
        insinfo = {"dat":[],"coln":["width","start","end","alistart","aliend"]}
    if method == "two":
        deltmp = []
        cortst = rulermatrix["cortst"]
        delgap = 0
        count = 0
        for i in range(0,len(cortst)):
            if cortst[i] != "-":
                count = count + 1
            if cortst[i] == "-" and delgap == 0:
                delgap = 1
                delalist = i
                delst = count
            if cortst[i] == "-" and delgap == 1:
                continue
            if cortst[i] != "-" and delgap == 1:
                delgap = 0
                delalien = i
                deltmp.append([delalien-delalist,delst+1,delst+delalien-delalist,delalist+1,delalien])
        ##-------------------------------------------------------------
        instmp = []
        corref = rulermatrix["corref"]
        insgap = 0
        count = 0
        for i in range(0,len(corref)):
            if corref[i] != "-":
                count = count + 1
            if corref[i] == "-" and insgap == 0:
                insgap = 1
                insalist = i
                insst = count
            if corref[i] == "-" and insgap == 1:
                continue
            if corref[i] != "-" and insgap == 1:
                insgap = 0
                insalien = i
                instmp.append([insalien-insalist,insst+1,insst+insalien-insalist,insalist+1,insalien])
        ##-------------------------------------------------------------
        delinfo = {"dat":deltmp,"coln":["width","start","end","alistart","aliend"]}
        insinfo = {"dat":instmp,"coln":["width","start","end","alistart","aliend"]}
        ##=============================================================
    sumgap = {"del":delinfo, "ins":insinfo}
    return sumgap

## re-estimate ruler position
def regRulerGap(regrulerd,rulergapd):
    regruler = regrulerd["dat"]
    rulername = regrulerd["coln"]
    rulergap = rulergapd["dat"]
    #alipos = pd.DataFrame(columns=["gap","accgap","alistart","aliend"])
    if len(rulergap) > 0:
        accg = []
        temg = []
        insti = listwhich(rulername,"==","instart")[0]
        ineni = listwhich(rulername,"==","inend")[0]
        for i in range(0,len(regruler)):
            st = regruler[i][insti]
            en = regruler[i][ineni]
            accgi = listwhich(selectcol(rulergapd,"start"),"<=",st)
            if len(accgi) == 0:
                accgtmp = 0
            else:
                accgtmp = sum(selectele(selectcol(rulergapd,"width"),accgi)) ## wrong with width, nnot defined yet
            j1 = listwhich(selectcol(rulergapd,"start"),">",st)
            j2 = listwhich(selectcol(rulergapd,"start"),"<",en)
            temgi = intersect(j1,j2)
            if len(temgi) == 0:
                temgtmp = 0
            else:
                temgtmp = sum(selectele(selectcol(rulergapd,"width"),temgi))
            accg.append(accgtmp)
            temg.append(temgtmp)
    if len(rulergap) == 0:
        accg = [0 for _ in range(0,len(regruler))]
        temg = [0 for _ in range(0,len(regruler))]
    
    insttmp = selectcol(regrulerd,"instart")
    inentmp = selectcol(regrulerd,"inend")
    alistart = listlistcal(insttmp,"+",accg)
    aliendtmp = listlistcal(inentmp,"+",accg)
    aliend = listlistcal(aliendtmp,"+",temg)
    
    alipos=[]
    alipos.append(temg)
    alipos.append(accg)
    alipos.append(alistart)
    alipos.append(aliend)
    aliposname = ["gap","accgap","alistart","aliend"]

    regruler_t = transpose(regruler)
    regrulerg_t = regruler_t + alipos
    regrulerg = transpose(regrulerg_t)
    reggname =rulername + aliposname

    regrulgd = {"dat":regrulerg,"coln":reggname}
    return regrulgd

## reunite
def reunite(hsr,ttr):
    st_hs = selectcol(hsr,"alistart")
    en_hs = selectcol(hsr,"aliend")
    st_tt = selectcol(ttr,"alistart")
    en_tt = selectcol(ttr,"aliend")
    ## re-organize range section
    store = []
    for i in range(0,len(st_hs)):
        store.append([st_hs[i],en_hs[i]])
    for j in range(0,len(st_tt)):
        store.append([st_tt[j],en_tt[j]])
    unistore = dropduplicated(store)
    unistored = {"dat":unistore,"coln":["alistart","aliend"]}
    liststore = listsort(unistored,by=["alistart","aliend"],increase=True)
    ## get final unique range
    if len(liststore["dat"]) == 1:
        sortedlist = liststore["dat"]    
    if len(liststore["dat"]) > 1:
        dropi = []
        for i in range(1,len(liststore["dat"])):
            if liststore["dat"][i][0] < liststore["dat"][i-1][1]:
                liststore["dat"][i][0] = min(liststore["dat"][i][0],liststore["dat"][i-1][0])
                liststore["dat"][i][1] = max(liststore["dat"][i][1],liststore["dat"][i-1][1])
                dropi.append(i-1)
        if len(dropi) > 0:
            sortedlist = listdrop(liststore["dat"],dropi)
        else:
            sortedlist = liststore["dat"]

    sortedd = {"dat":sortedlist,"coln":["alistart","aliend"]}
    ##  re-arrange group name
    sti = listwhich(hsr["coln"],"==","alistart")[0]
    eni = listwhich(hsr["coln"],"==","aliend")[0]
    gpi = listwhich(hsr["coln"],"==","group")[0]
    for i in range(0,len(hsr["dat"])):
        j1 = listwhich(selectcol(sortedd,"alistart"),"<=",hsr["dat"][i][sti])
        j2 = listwhich(selectcol(sortedd,"aliend"),">=",hsr["dat"][i][eni])
        hsr["dat"][i][gpi] = intersect(j1,j2)[0] + 1

    for i in range(0,len(ttr["dat"])):
        j1 = listwhich(selectcol(sortedd,"alistart"),"<=",ttr["dat"][i][sti])
        j2 = listwhich(selectcol(sortedd,"aliend"),">=",ttr["dat"][i][eni])
        ttr["dat"][i][gpi] = intersect(j1,j2)[0] + 1
    
    infosum = {"hs":hsr,"tt":ttr,"group":sortedd}
    return infosum

## align position to exon
def assignReg(regrulerd,regidxd,rulergapd,genstr):
    regidxd = typeinto(regidxd,"start",int)
    regidxd = typeinto(regidxd,"end",int)
    if genstr == "+":
        regidxd = listsort(regidxd,by=["start","end"],increase = True)
    if genstr == "-":
        regidxd = listsort(regidxd,by=["start","end"],increase = False)
    regruler = regrulerd["dat"]
    rulername = regrulerd["coln"]
    rulergap = rulergapd["dat"]
    regidx = regidxd["dat"]
    regname = regidxd["coln"]  
    
    instart = [0 for _ in range(0,len(regidx))]
    inend = [0 for _ in range(0,len(regidx))]
    sti = listwhich(rulername,"==","start")[0]
    eni = listwhich(rulername,"==","end")[0]
    insti = listwhich(rulername,"==","instart")[0]
    ineni = listwhich(rulername,"==","inend")[0]
    groupi = listwhich(rulername,"==","group")[0]
    regsti = listwhich(regname,"==","start")[0]
    regeni = listwhich(regname,"==","end")[0]
    for i in range(0,len(regruler)):
        prst = regruler[i][sti]
        pren = regruler[i][eni]
        inst = regruler[i][insti]
        inen = regruler[i][ineni]
        j1 = listwhich(selectcol(regidxd,"start"),">=",prst)
        j2 = listwhich(selectcol(regidxd,"end"),"<=",pren)
        selecti = intersect(j1,j2)
        
        for j in range(0,len(selecti)):
            refst = regidx[selecti[j]][regsti]
            refen = regidx[selecti[j]][regeni]
            if genstr == "+":
                instart[selecti[j]] = inst + abs(prst-refst)
                inend[selecti[j]] = inen - abs(pren-refen)
            if genstr == "-":
                instart[selecti[j]] = inst + abs(pren-refen)
                inend[selecti[j]] = inen - abs(prst-refst)
    group = [0 for _ in range(0,len(regidx))]
    for i in range(0,len(regidx)):
        j1 = listwhich(selectcol(regrulerd,"instart"),"<=",instart[i])
        j2 = listwhich(selectcol(regrulerd,"inend"),">=",inend[i])
        selecti = intersect(j1,j2)
        if len(selecti) > 0:
            group[i] = regruler[selecti[0]][groupi]

    if len(rulergap) > 0:
        accg = []
        temg = []
        for i in range(0,len(instart)):
            st = instart[i]
            en = inend[i]
            accgi = listwhich(selectcol(rulergapd,"start"),"<=",st)
            if len(accgi) == 0:
                accgtmp = 0
            else:
                accgtmp = sum(selectele(selectcol(rulergapd,"width"),accgi)) ## wrong with width, nnot defined yet
            j1 = listwhich(selectcol(rulergapd,"start"),">",st)
            j2 = listwhich(selectcol(rulergapd,"start"),"<",en)
            temgi = intersect(j1,j2)
            if len(temgi) == 0:
                temgtmp = 0
            else:
                temgtmp = sum(selectele(selectcol(rulergapd,"width"),temgi))
            accg.append(accgtmp)
            temg.append(temgtmp)
    if len(rulergap) == 0:
        accg = [0 for _ in range(0,len(regidx))]
        temg = [0 for _ in range(0,len(regidx))]
    
    regidxdata = []
    regidxdata.append(group)
    regidxdata.append(selectcol(regidxd,"id"))
    regidxdata.append(selectcol(regidxd,"seq"))
    regidxdata.append(temg)
    regidxdata.append(accg)
    alistart = listlistcal(instart,"+",accg)
    regidxdata.append(alistart)
    aliendtmp = listlistcal(inend,"+",accg)
    aliend = listlistcal(aliendtmp,"+",temg)
    regidxdata.append(aliend)
    regsizetmp = listlistcal(inend,"-",instart)
    regsize = listcal(regsizetmp,"+",1)
    regidxdata.append(regsize)

    newregidx = transpose(regidxdata)
    sumnew = {"dat":newregidx,"coln":["group","id","seq","gap","accgap","alistart","aliend","size"]}
    return sumnew

## define extract canvas information
def extractCanvasTwo(rulerd):
    rulersize = len(rulerd["corref"])
    cavinfo = [["CANVAS",1,rulersize,0,0,"canvas","canvas"]]
    cansum = {"dat":cavinfo,"coln":["type","x","wid","y","gcol","cli","legend"]}
    return cansum

## extract ruler information
def extractRulerTwo(groupd,hsruler,rulerd,ttruler,detailspecies):
    species1 = detailspecies.split("-")[0]
    species2 = detailspecies.split("-")[1]
    ##
    store = []
    ## organize group tag
    yitmp = 1
    gsti = listwhich(groupd["coln"],"==","alistart")[0]
    geni = listwhich(groupd["coln"],"==","aliend")[0]

    hsgrpseq = rulerd["corref"]
    ttgrpseq = rulerd["cortst"]

    for i in range(0,len(groupd["dat"])):
        start = groupd["dat"][i][gsti]
        end = groupd["dat"][i][geni]
        width = end - start # + 1
        grtmp = i + 1
        tmp = ["T",round(width*0.5+start,0),0,yitmp+1,grtmp,"EGIT",str(grtmp)]
        store.append(tmp)
        ##----------------------------------
        grphstmp = hsgrpseq[(start-1):end]
        grptttmp = ttgrpseq[(start-1):end]
        tmp = ["P",start,width,yitmp+0.7,"EGI"+str(grtmp),"EGI","EG"+str(grtmp)]
        store.append(tmp)
        tmp = ["P",start,width,yitmp+1.3,"EGI"+str(grtmp),"EGI","EG"+str(grtmp)]
        store.append(tmp)
        tmp = ["P",start+width,width,yitmp+1.3,"EGI"+str(grtmp),"EGI","EG"+str(grtmp)]
        store.append(tmp)
        tmp = ["P",start+width,width,yitmp+0.7,"EGI"+str(grtmp),"EGI","EG"+str(grtmp)]
        store.append(tmp)
    store.append(["T",int(start),int(width),yitmp+1,grtmp,"RAWT","EGI"])
    ## organize hs
    yitmp = 3
    rsti = listwhich(hsruler["coln"],"==","alistart")[0]
    reni = listwhich(hsruler["coln"],"==","aliend")[0]
    rgri = listwhich(hsruler["coln"],"==","group")[0]
    ridi = listwhich(hsruler["coln"],"==","id")[0]
    for i in range(0,len(hsruler["dat"])):
        start = hsruler["dat"][i][rsti]
        width = hsruler["dat"][i][reni] - hsruler["dat"][i][rsti]# + 1
        grtmp = hsruler["dat"][i][rgri]
        antmp = hsruler["dat"][i][ridi]
        tmp = ["P",int(start),int(width),yitmp+0.7,"UIh"+str(grtmp),"UIh","EI:"+species1+"-"+species2]
        store.append(tmp)
        tmp = ["P",int(start),int(width),yitmp+1.3,"UIh"+str(grtmp),"UIh","EI:"+species1+"-"+species2]
        store.append(tmp)
        tmp = ["P",int(start)+int(width),int(width),yitmp+1.3,"UIh"+str(grtmp),"UIh","EI:"+species1+"-"+species2]
        store.append(tmp)
        tmp = ["P",int(start)+int(width),int(width),yitmp+0.7,"UIh"+str(grtmp),"UIh","EI:"+species1+"-"+species2]
        store.append(tmp)
    store.append(["T",0,0,yitmp+1,grtmp,"RAWT","Exon ideo. "+species1])
    ## organize primate
    yitmp = 4
    for i in range(0,len(ttruler["dat"])):
        start = ttruler["dat"][i][rsti]
        width = ttruler["dat"][i][reni] - ttruler["dat"][i][rsti]# + 1
        grtmp = ttruler["dat"][i][rgri]
        antmp = ttruler["dat"][i][ridi]
        tmp = ["P",int(start),int(width),yitmp+0.7,"UIt"+str(grtmp),"UIt","EI:"+species1+"-"+species2]
        store.append(tmp)
        tmp = ["P",int(start),int(width),yitmp+1.3,"UIt"+str(grtmp),"UIt","EI:"+species1+"-"+species2]
        store.append(tmp)
        tmp = ["P",int(start)+int(width),int(width),yitmp+1.3,"UIt"+str(grtmp),"UIt","EI:"+species1+"-"+species2]
        store.append(tmp)
        tmp = ["P",int(start)+int(width),int(width),yitmp+0.7,"UIt"+str(grtmp),"UIt","EI:"+species1+"-"+species2]
        store.append(tmp)
    store.append(["T",0,0,yitmp+1,grtmp,"RAWT","Exon ideo. "+species2])
    tmp = ["L",0,0,yitmp+1,0,"IL","IL"]
    store.append(tmp)
    storesum = {"dat":dropduplicated(store),"coln":["type","x","wid","y","gcol","cli","legend"]}
    return storesum

## extract isoform information
def extractISOTwo(isoalid,gensumd,species):

    store = []
    orfgroup = []
    isoname = isoalid["coln"][8:]

    isti = listwhich(isoalid["coln"],"==","alistart")[0]
    ieni = listwhich(isoalid["coln"],"==","aliend")[0]
    igri = listwhich(isoalid["coln"],"==","group")[0]
    iidi = listwhich(isoalid["coln"],"==","id")[0]
    plotgrp = 0
    for i in range(0,len(isoname)):
        comstatus = selectcol(isoalid,isoname[i])
        entmp = isoname[i]
        for j in range(0,len(comstatus)):
            plotgrp = plotgrp + 1
            if comstatus[j] == 0:
                continue
            start = isoalid["dat"][j][isti]
            width = isoalid["dat"][j][ieni] - isoalid["dat"][j][isti] #+ 1
            yitmp = i + 1
            grtmp = isoalid["dat"][j][igri]

            tmp = ["P",int(start),int(width),yitmp+0.7,species+str(plotgrp),"EXON"+str(grtmp),"exon_"+species]
            store.append(tmp)
            tmp = ["P",int(start),int(width),yitmp+1.3,species+str(plotgrp),"EXON"+str(grtmp),"exon_"+species]
            store.append(tmp)
            tmp = ["P",int(start)+int(width),int(width),yitmp+1.3,species+str(plotgrp),"EXON"+str(grtmp),"exon_"+species]
            store.append(tmp)
            tmp = ["P",int(start)+int(width),int(width),yitmp+0.7,species+str(plotgrp),"EXON"+str(grtmp),"exon_"+species]
            store.append(tmp)
            
            texttmp = ["T",0,0,yitmp+1,0,"RAWT",entmp]

            store.append(texttmp)
            #store.append(tmp)
    
    tmp = ["L",0,0,yitmp+1,0,"IL","IL"]
    store.append(tmp)

    storesum = {"dat":dropduplicated(store),"coln":["type","x","wid","y","gcol","cli","legend"]}
    return storesum

## summary sumdata
def summaryplotTwo(canvasd, rulerd,isodathsd,isodatttd):
    y1 = max(selectcol(rulerd,"y"))
    y2 = max(selectcol(isodathsd,"y"))
    y3 = max(selectcol(isodatttd,"y"))

    ytmp = y1 + y2 + y3 + 10
    #ytmp = round(selectcol(canvasd,"y")[0]+1,0)
    
    if len(rulerd["dat"]) > 0:
        rulerd = listColCal(rulerd,"y","-",ytmp)
        rulerd = listColCal(rulerd,"y","*",-1)
        ytmp = min(selectcol(rulerd,"y"))    
    if len(isodathsd["dat"]) > 0:
        isodathsd = listColCal(isodathsd,"y","-",ytmp)
        isodathsd = listColCal(isodathsd,"y","*",-1)
        ytmp = min(selectcol(isodathsd,"y"))
    if len(isodatttd["dat"]) > 0:
        isodatttd = listColCal(isodatttd,"y","-",ytmp)
        isodatttd = listColCal(isodatttd,"y","*",-1)
        ytmp = min(selectcol(isodatttd,"y"))
    
    yi = listwhich(canvasd["coln"],"==","y")[0]
    canvasd["dat"][0][yi] = ytmp

    sumdat = rulerd["dat"] + isodathsd["dat"] + isodatttd["dat"]
    
    listsum = {"dat":sumdat,"coln":rulerd["coln"]}
    return listsum

##======================================================##
##                   start analysis here                ##
##======================================================##

## main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract exon mappings from reciprocal BLASTN')
    parser.add_argument('--speciespair',
        nargs='?',
        help='species pair, eg. species1-species2')
    parser.add_argument('--orthologpair',
        nargs='?',
        help='ensembl gene pair, eg. ensembl_gene1-ensembl_gene2')
    parser.add_argument('--outpath',
        nargs='?',
        help='ensembl gene pair, eg. ensembl_gene1-ensembl_gene2')

    args = parser.parse_args()
    ##-------------------------------------------
    wdpath = os.getcwd()
    codonstr = "AAA#Lys#K@AAT#Asn#N@AAG#Lys#K@AAC#Asn#N@ATA#Ile#I@ATT#Ile#I@ATG#Met#M@ATC#Ile#I@AGA#Arg#R@AGT#Ser#S@AGG#Arg#R@AGC#Ser#S@ACA#Thr#T@ACT#Thr#T@ACG#Thr#T@ACC#Thr#T@TAA#Ter#*@TAT#Tyr#Y@TAG#Ter#*@TAC#Tyr#Y@TTA#Leu#L@TTT#Phe#F@TTG#Leu#L@TTC#Phe#F@TGA#Ter#*@TGT#Cys#C@TGG#Trp#W@TGC#Cys#C@TCA#Ser#S@TCT#Ser#S@TCG#Ser#S@TCC#Ser#S@GAA#Glu#E@GAT#Asp#D@GAG#Glu#E@GAC#Asp#D@GTA#Val#V@GTT#Val#V@GTG#Val#V@GTC#Val#V@GGA#Gly#G@GGT#Gly#G@GGG#Gly#G@GGC#Gly#G@GCA#Ala#A@GCT#Ala#A@GCG#Ala#A@GCC#Ala#A@CAA#Gln#Q@CAT#His#H@CAG#Gln#Q@CAC#His#H@CTA#Leu#L@CTT#Leu#L@CTG#Leu#L@CTC#Leu#L@CGA#Arg#R@CGT#Arg#R@CGG#Arg#R@CGC#Arg#R@CCA#Pro#P@CCT#Pro#P@CCG#Pro#P@CCC#Pro#P"
    ##------------ prepare data -----------------
    speciespair = args.speciespair
    orthologpair = args.orthologpair
    
    species1, species2 = speciespair.split("-")
    test = testquery(speciespair, orthologpair, wdpath)
    ##-------------------------------------------
    geninfohs = aslsls(test['genestr1'], typename="isoinfo")
    cdshs = aslsls(test['cdsstr1'], typename="cds")
    geninfott = aslsls(test['genestr2'], typename="isoinfo")
    cdstt = aslsls(test['cdsstr2'], typename="cds")

    strandhs = selectcol(geninfohs,"strand")[0]
    strandtt = selectcol(geninfott,"strand")[0]

    codon = aslsls(codonstr, typename="codon")
    exongroup = aslsls(test['egstr'], typename="exong")
    ##============ estimate cds ruler and gap =================
    cdsmtt = regComSum(genestd=geninfott,regidxd=cdstt)
    cdsinfott = newinfo(comd=cdsmtt,regidxd=cdstt)
    cdsmhs = regComSum(genestd=geninfohs,regidxd=cdshs)
    cdsinfohs = newinfo(comd=cdsmhs,regidxd=cdshs)
    ##------------ estimate cds and utr ruler -----------------
    cdsrulerhs = regRuler(regidxd=cdshs,genstr=strandhs)
    cdsrulertt = regRuler(regidxd=cdstt,genstr=strandtt)
    ##=-=-=-=-=-=-=-=-=-=-= optimize R1 =-=-=-=-=-=-=-=-=-=-=
    #regrhs=cdsrulerhs;regrtt=cdsrulertt;exong = exongroup
    rulermtwo = regAlignTwo(regrhs=cdsrulerhs,regrtt=cdsrulertt,exong = exongroup)
    ##------------- estimate cds and utr gap ------------------
    cdsgap = estiGap(rulermatrix=rulermtwo,method="two")
    cdsdel = cdsgap["del"]
    cdsins = cdsgap["ins"]
    cdshsr = regRulerGap(regrulerd=cdsrulerhs,rulergapd=cdsgap["ins"])
    cdsttr = regRulerGap(regrulerd=cdsrulertt,rulergapd=cdsgap["del"])
    ##------------------ reunite -----------------------------
    unitedinfo = reunite(hsr=cdshsr,ttr=cdsttr)
    hsr = unitedinfo["hs"]
    ttr = unitedinfo["tt"]
    groupinfo = unitedinfo["group"]
    ##========== estimate aligned exon position ===============
    alicdshs = assignReg(regrulerd=hsr,regidxd=cdshs,rulergapd=cdsins,genstr=strandhs)
    alicdstt = assignReg(regrulerd=ttr,regidxd=cdstt,rulergapd=cdsdel,genstr=strandtt)

    isoinfohs = listmerge(listd1=alicdshs, listd2=cdsmhs, by="id")
    isoinfott = listmerge(listd1=alicdstt, listd2=cdsmtt, by="id")
    ##=============== prepare plot data =======================
    canff = extractCanvasTwo(rulerd = rulermtwo)
    rulerff = extractRulerTwo(groupd=groupinfo,rulerd = rulermtwo,hsruler=hsr,ttruler=ttr,detailspecies = speciespair)

    isoffhs = extractISOTwo(isoalid=isoinfohs,gensumd=geninfohs,species=species1)
    isofftt = extractISOTwo(isoalid=isoinfott,gensumd=geninfott,species=species2)
    ##=============== summarize plot data =====================
    sumff = summaryplotTwo(canvasd=canff, rulerd=rulerff,isodathsd=isoffhs,isodatttd=isofftt)
    ##======================================================##
    ##                    end analysis here                 ##
    ##======================================================##

    pyggplot = pd.DataFrame(sumff["dat"])
    pyggplot.columns = sumff["coln"]

    egiexon = pyggplot[pyggplot.type=="P"]
    isoname = pyggplot[pyggplot.cli=="RAWT"]
    egitext = pyggplot[pyggplot.cli=="EGIT"]

    canvas_y = max(pyggplot["y"])*0.15
    if canvas_y > 8:
        canvas_y = 8

    p = ggplot(egiexon,aes(x="x",y="y",group="gcol",fill="legend")) + geom_polygon(colour ="black")# + ylim((0,25))
    p = p + scale_y_continuous(breaks = list(isoname["y"]), labels = list(isoname["legend"]))
    p = p + annotate("text",label=list(egitext["legend"]), x=list(egitext["x"]), y=list(egitext["y"]), size = 6, colour = "black")
    p = p + xlab(" vs ".join(orthologpair.split("-")))
    p = p + theme(legend_position="none")
    p = p + theme(axis_text_x = element_blank(),
			   axis_text_y = element_text(color="black"),
               axis_ticks= element_blank(),
               axis_title_y= element_blank(),
			   axis_title_x= element_text(color="black"),
               figure_size=(6.5, canvas_y)
			   )

    if args.outpath:
        if args.outpath[-1] == "/":
            p.save(args.outpath + orthologpair+".pdf", width=6.5, height=canvas_y)
        else:
            p.save(args.outpath + "/" + orthologpair+".pdf", width=6.5, height=canvas_y)
    else:
        p.save(orthologpair+".pdf", width=6.5, height=canvas_y)


# cd /Users/jeffma/EGIO-main
# python __plotIsoCom.py --speciespair hsa-mmu --orthologpair ENSG00000002330-ENSMUSG00000024959 --outpath /Users/jeffma/EGIO-main/plot
# python __plotIsoCom.py --speciespair hsa-mmu --orthologpair ENSG00000099889-ENSMUSG00000000325 --outpath /Users/jeffma/EGIO-main/plot
# python __plotIsoCom.py --speciespair hsa-mmu --orthologpair ENSG00000099968-ENSMUSG00000009112 --outpath /Users/jeffma/EGIO-main/plot



## speciespair = "hsa-mmu"
## orthologpair = "ENSG00000002330-ENSMUSG00000024959"