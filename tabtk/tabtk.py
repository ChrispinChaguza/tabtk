#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
import multiprocessing
import argparse
import itertools
from progress.bar import ChargingBar
import re
import subprocess
from pathlib import Path
import numpy as np
import shutil
import pandas as pd
import networkx as nx
import math
import datetime
import statistics
import glob
from .__init__ import __version__

def pairwiseSimilarity(firstRecord,secondRecord):
    recordMatch = 0
    recordMisMatch = 0

    firstRecordName = list(firstRecord.keys())[0]
    secondRecordName = list(secondRecord.keys())[0]

    for r,s in zip(firstRecord[firstRecordName],secondRecord[secondRecordName]):
        if not (r == 0 and s == 0):
            if r == s:
                if (r in [0,1] and s in [0,1]):
                    recordMatch = recordMatch + 1
                else:
                    print(f"--->>>Found characters other than 0 and 1 when comparing items {firstRecordName} and {secondRecordName}... exiting...")
                    sys.exit()
            else:
                if (r in [0,1] and s in [0,1]):
                    recordMisMatch = recordMisMatch + 1
                else:
                    print(f"--->>>Found characters other than 0 and 1 when comparing items {firstRecordName} and {secondRecordName}... exiting...")
                    sys.exit()
        else:
            if (r in [0,1] and s in [0,1]):
                pass
            else:
                print(f"--->>>Found characters other than 0 and 1 when comparing items {firstRecordName} and {secondRecordName}... exiting...")
                sys.exit()

    percentSimilarity = recordMatch/(recordMatch + recordMisMatch)*100

    similarityData = [firstRecordName,secondRecordName,percentSimilarity]

    #print('.', end='', flush=True)

    return(similarityData)


def alignmentSimilarity(ij):
    seqMatch=0
    seqMisMatch=0

    for r,s in zip(str(i.seq).replace("-","N").replace("X","N").upper(),str(j.seq).replace("-","N").replace("X","N").upper()):
        if r==s and not (r=="N" or s=="N"):
            seqMatch=seqMatch+1
        else:
            if r=="N" or s=="N":
                continue
            else:
                seqMisMatch=seqMisMatch+1

    seqPID=round(seqMatch/(seqMatch+seqMisMatch)*100,4) if (seqMatch+seqMisMatch)!=0 else 0

    seqLen=round(len(i.seq),4)

    seqCov=round((seqMatch+seqMisMatch)/seqLen*100,4) if seqLen!=0 else 0

    retValue = f"{i.id}\t{j.id}\t{seqPID}"

    return(retValue)


def mashSimilarity(cmdValues):
    if not shutil.which("mash"):
        print("Install MASH and try again... exiting...")
        sys.exit()

    else:
        pass

    tmpOut = f"mash.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"
    tmpOutR = f"mash.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

    print("Calculating percent similarity... may take longer... be patient...")

    mashSketchCmd = f"mash sketch -p {cmdValues["threads"]} -k {cmdValues["length"]} -o {tmpOut} -l {cmdValues["data"]} -s {cmdValues["size"]}"
    mashDistCmd = f"mash dist -p {cmdValues["threads"]} -s {cmdValues["size"]} {f"{tmpOut}.msh"} {f"{tmpOut}.msh"} > {tmpOutR}"

    subprocess.call(mashSketchCmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    subprocess.call(mashDistCmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    mashResults = []

    with open(tmpOutR,"r") as mashFile:
        for mashOutput in mashFile:
            item = mashOutput.strip().replace("/","\t").split("\t")
            mashResults.append(f"{Path(item[0]).stem}\t{Path(item[1]).stem}\t{float(item[4])/float(item[5])*100}")

    if os.path.exists(tmpOut):
        os.remove(tmpOut)
    else:
        pass

    if os.path.exists(tmpOutR):
        os.remove(tmpOutR)
    else:
        pass

    if os.path.exists(f"{tmpOut}.msh"):
        os.remove(f"{tmpOut}.msh")
    else:
        pass

    return(mashResults)    

def blastSimilarity(cmdValues):
    if not shutil.which("blastn"):
        print("Install BLAST and try again... exiting...")
        sys.exit()

    else:
        pass

    tmpFileNames = []
    fileNames = []

    tmpCmd="wc -l "+str(cmdValues["data"])+" | awk '{print $1}'"
    dataFileLines=subprocess.Popen(tmpCmd,shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    dataLines = int(dataFileLines.communicate()[0])
    dataFileLines.kill()

    bar = ChargingBar("Loading data", max=int(dataLines))
    for fileName in open(cmdValues["data"],"r"):
        fileNames.append(str(fileName).strip())

        bar.next()
    bar.finish()

    for i in fileNames:
        tmpFileNames.append(f"{Path(i).stem}.tmpfna")

        with open(f"{Path(i).stem}.tmpfna","w") as fh:
            for j in SeqIO.parse(i,"fasta"):
                fh.write(f">{Path(os.path.basename(i)).stem}\n{j.seq}\n") 

    tmpOut = f"blast.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

    Cmd1 = f"cat {" ".join(tmpFileNames)} > {tmpOut}"
    Cmd2 = f"makeblastdb -in {tmpOut} -dbtype nucl"

    subprocess.call(Cmd1,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    subprocess.call(Cmd2,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    print("Calculating percent similarity... may take longer... be patient...")

    tmpOut1 = f"blast.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"
    tmpOut2 = f"blast.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

    blastCmd1 = f"echo \'qseqid sseqid length pident qcovs\' | tr \' \' \'\\t\' > {tmpOut2}"
    blastCmd2 = f"blastn -db {tmpOut} -query {tmpOut} -num_threads {cmdValues["threads"]} -outfmt \'6 qseqid sseqid length pident qcovs\' -evalue 0.001 >> {tmpOut2}"

    subprocess.call(blastCmd1,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    subprocess.call(blastCmd2,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    blastData = pd.read_csv(tmpOut2,delimiter="\t",dtype={'pident': float,'qcovs': float})
    blastData["tmp"] = blastData.apply(lambda x: (x["qcovs"] * x["pident"])/100,axis=1)

    blastData["qseqidR"] = blastData.apply(lambda x: x["qseqid"] if x["qseqid"] > x["sseqid"] else x["sseqid"],axis=1)
    blastData["sseqidR"] = blastData.apply(lambda x: x["qseqid"] if x["qseqid"] < x["sseqid"] else x["sseqid"],axis=1)

    blastData["qseqid"] = blastData["qseqidR"]
    blastData["sseqid"] = blastData["sseqidR"]

    blastData["cov"] = blastData[["qseqid","sseqid","tmp"]].groupby(["qseqid","sseqid"])["tmp"].transform(lambda x: statistics.mean(x))
    blastDataR = blastData.groupby(["qseqid","sseqid"]).first().reset_index()

    blastData = blastDataR.drop(columns=["tmp","length","pident","qseqidR","sseqidR"])
    blastDataFinal = blastData[["sseqid","qseqid","cov"]].rename(columns={"sseqid":"seq1","qseqid":"seq2","cov":"pid"}).values.tolist()

    bar = ChargingBar("Cleaning temporary files", max=len(fileNames))
    for fileName in tmpFileNames:
        if os.path.exists(fileName):
            os.remove(fileName)
        else:
            pass
            
        bar.next() 
    bar.finish() 

    if os.path.exists(tmpOut1):
        os.remove(tmpOut1)
    else:
        pass

    if os.path.exists(tmpOut2):
        os.remove(tmpOut2)
    else:
        pass

    for file in glob.glob(f"{tmpOut}*"):
        if os.path.exists(file):
            os.remove(file)
        else:
            pass

    return(blastDataFinal)


def fastaniSimilarity(cmdValues):
    if not shutil.which("fastANI"):
        print("Install fastANI and try again... exiting...")
        sys.exit()

    else:
        pass

    tmpOut = f"mash.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

    fastaniCmd = f"fastANI --rl {cmdValues["data"]} --ql {cmdValues["data"]} -k {cmdValues["length"]} --matrix --fragLen {cmdValues["fragment"]} --minFraction {cmdValues["fraction"]} --output {tmpOut} -t {cmdValues["threads"]}"

    print("Calculating percent similarity... may take longer... be patient...")

    subprocess.call(fastaniCmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    fastaniResults = []

    with open(tmpOut,"r") as fastaniFile:
        for itemPos,fastaniOutput in enumerate(fastaniFile):
            if itemPos == 0:
                pass
            else:
                item = fastaniOutput.strip().split("\t")
                fastaniResults.append(f"{Path(item[0]).stem}\t{Path(item[1]).stem}\t{float(item[2])}")

    if os.path.exists(tmpOut):
        os.remove(tmpOut)
    else:
        pass

    if os.path.exists(f"{tmpOut}.matrix"):
        os.remove(f"{tmpOut}.matrix")
    else:
        pass

    return(fastaniResults)


def similarity(cmdValues):

    if sum([1 if cmdValues["aln"] else 0, 1 if cmdValues["mash"] else 0, 1 if cmdValues["fastani"] else 0]) > 1:
        print("Specify none or only one of these options --aln/-a, --mash/-m, and --fastani/-n... exiting...")
        sys.exit()
    else:
        pass

    if (not cmdValues["blast"]) and (not cmdValues["aln"]) and (not cmdValues["mash"]) and (not cmdValues["fastani"]):
        dataMatrixDF = loadDataTable(cmdValues)

        itemPairs = itertools.permutations([item for pos,item in enumerate(dataMatrixDF.keys()) if pos>0], r=2)
        tmpPairs = itertools.permutations([item for pos,item in enumerate(dataMatrixDF.keys()) if pos>0], r=2)

        print("Calculating pairwise similarity... may take longer... be patient...")

        with multiprocessing.Pool(processes=cmdValues["threads"]) as pool:
            args = [({firstItem: dataMatrixDF[firstItem]}, {secondItem: dataMatrixDF[secondItem]}) for firstItem,secondItem in itemPairs]
            similarityResults = pool.starmap(pairwiseSimilarity, args)

        with open(cmdValues["output"],"w") as simData:
            simData.write("seq1\tseq2\tpid\n")

            bar = ChargingBar('Writing similarity values', max=len(similarityResults))
            for itemPair in similarityResults:
                simData.write(f"{itemPair[0]}\t{itemPair[1]}\t{itemPair[2]}\n")

                bar.next()
            bar.finish()

    elif (not cmdValues["blast"]) and (not cmdValues["aln"]) and (cmdValues["mash"]) and (not cmdValues["fastani"]):
        similarityResults = mashSimilarity(cmdValues)

        with open(cmdValues["output"],"w") as simData:
            simData.write("seq1\tseq2\tpid\n")

            bar = ChargingBar('Writing similarity values', max=len(similarityResults))
            for itemPair in similarityResults:
                simData.write(f"{itemPair}\n")

                bar.next()
            bar.finish()

    elif (not cmdValues["blast"]) and (not cmdValues["aln"]) and (not cmdValues["mash"]) and (cmdValues["fastani"]):
        similarityResults = fastaniSimilarity(cmdValues)

        with open(cmdValues["output"],"w") as simData:
            simData.write("seq1\tseq2\tpid\n")

            bar = ChargingBar('Writing similarity values', max=len(similarityResults))
            for itemPair in similarityResults:
                simData.write(f"{itemPair}\n")

                bar.next()
            bar.finish()

    elif (not cmdValues["blast"]) and (cmdValues["aln"]) and (not cmdValues["mash"]) and (not cmdValues["fastani"]):
        alignment={}

        for i in SeqIO.parse(cmdValues['data'],"fasta"):
            alignment[i.id]=i

        fhandle=open(str(cmdValues['outputFile']),"w")
        fhandle.write("seq1\tseq2\tpid\n")

        tmpList = sorted(map(sorted, combinations(set([i for i in alignment.keys()]), 2)))

        with multiprocessing.Pool(processes=cmdValues['threads']) as pool:
            args=[(alignment[m],alignment[n]) for m,n in tmpList]
            results=pool.starmap(alignmentSimilarity, args)

    elif (cmdValues["blast"]) and (not cmdValues["aln"]) and (not cmdValues["mash"]) and (not cmdValues["fastani"]):
        similarityResults = blastSimilarity(cmdValues)

        with open(cmdValues["output"],"w") as simData:
            simData.write("seq1\tseq2\tpid\n")

            bar = ChargingBar('Writing similarity values', max=len(similarityResults))
            for itemPair in similarityResults:
                simData.write(f"{"\t".join([str(r) for r in itemPair])}\n")

                bar.next()
            bar.finish()

    else:
        pass


def transpose(cmdValues):
    if not shutil.which("datamash"):
        print("Install 'datamash' software and try again... exiting...")
        sys.exit()

    else:
        pass

    bar = ChargingBar("Transposing table", max=1)
    datamashCmd = f"datamash transpose < {cmdValues["data"]} > {f"{Path(cmdValues["output"]).stem}.transposed.tsv"}"

    subprocess.call(datamashCmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    bar.next()
    bar.finish()


def lineages1Levels(cmdValues):
    np.random.rand(1)

    print(cmdValues["thresholds"][0])
    snpDistDF = pd.read_csv(cmdValues["data"],delimiter="\t",header=0)

    Graph = nx.Graph()
    Graph.add_nodes_from(snpDistDF.seq2.to_list())
    Graph.add_nodes_from(snpDistDF.seq1.to_list())

    snpDistDF1 = snpDistDF[snpDistDF.pid >= cmdValues["thresholds"][0]]
    seqPair = [(snpDistDF1.seq1[i],snpDistDF1.seq2[i]) for i in snpDistDF1.index]

    Graph.add_edges_from(seqPair)
    networkClusters = tuple(nx.connected_components(Graph))

    fhandle = open(f"{cmdValues['output']}","w")

    for i,j in enumerate(networkClusters):
        for s in j:
            FirstLevelLin = str(i+1)
            seqLabel = s
            clusterPrefix = cmdValues["prefix"]
            separator = cmdValues["separator"]

            fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}\n")

    fhandle.close();


def lineages2Levels(cmdValues):
    np.random.rand(1)

    snpDistDF = pd.read_csv(cmdValues["data"],delimiter="\t",header=0)

    Graph = nx.Graph()
    Graph.add_nodes_from(snpDistDF.seq2.to_list())
    Graph.add_nodes_from(snpDistDF.seq1.to_list())

    snpDistDF1 = snpDistDF[snpDistDF.pid >= cmdValues["thresholds"][0]]
    seqPair = [(snpDistDF1.seq1[i],snpDistDF1.seq2[i]) for i in snpDistDF1.index]

    Graph.add_edges_from(seqPair)
    networkClusters = tuple(nx.connected_components(Graph))

    fhandle = open(f"{cmdValues['output']}","w")

    for i,j in enumerate(networkClusters):
        snpDistDF2 = snpDistDF[snpDistDF.seq1.isin(j) & snpDistDF.seq2.isin(j)]
        snpDistDF3 = snpDistDF2[snpDistDF2.pid >= cmdValues["thresholds"][1]]

        Graph1 = nx.Graph()
        Graph1.add_nodes_from(j)

        seqPair1 = [(snpDistDF3.seq1[z],snpDistDF3.seq2[z]) for z in snpDistDF3.index]

        Graph1.add_edges_from(seqPair1)
        networkClusters1 = tuple(nx.connected_components(Graph1))

        for r,k in enumerate(networkClusters1):
            for s in k:
                FirstLevelLin = str(i+1)
                seqLabel = s
                SecondLevelLin = str(r+1)
                clusterPrefix = cmdValues["prefix"]
                separator = cmdValues["separator"]

                if len(networkClusters1)==1:
                   fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}\n")
                else:
                    fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}.{SecondLevelLin}\n")

    fhandle.close();


def lineages3Levels(cmdValues):
    np.random.rand(1)

    snpDistDF = pd.read_csv(cmdValues["data"],delimiter="\t",header=0)

    Graph = nx.Graph()
    Graph.add_nodes_from(snpDistDF.seq2.to_list())
    Graph.add_nodes_from(snpDistDF.seq1.to_list())

    snpDistDF1 = snpDistDF[snpDistDF.pid >= cmdValues["thresholds"][0]]
    seqPair = [(snpDistDF1.seq1[i],snpDistDF1.seq2[i]) for i in snpDistDF1.index]

    Graph.add_edges_from(seqPair)
    networkClusters = tuple(nx.connected_components(Graph))

    fhandle = open(f"{cmdValues['output']}","w")

    for i,j in enumerate(networkClusters):
        snpDistDF2 = snpDistDF[snpDistDF.seq1.isin(j) & snpDistDF.seq2.isin(j)]
        snpDistDF3 = snpDistDF2[snpDistDF2.pid >= cmdValues["thresholds"][1]]

        Graph1 = nx.Graph()
        Graph1.add_nodes_from(j)

        seqPair1 = [(snpDistDF3.seq1[z],snpDistDF3.seq2[z]) for z in snpDistDF3.index]

        Graph1.add_edges_from(seqPair1)
        networkClusters1 = tuple(nx.connected_components(Graph1))

        for p,q in enumerate(networkClusters1):
            snpDistDF2R = snpDistDF[snpDistDF.seq1.isin(q) & snpDistDF.seq2.isin(q)]
            snpDistDF3R = snpDistDF2R[snpDistDF2R.pid >= cmdValues["thresholds"][2]]

            Graph2 = nx.Graph()
            Graph2.add_nodes_from(q)

            seqPair1R = [(snpDistDF3R.seq1[z],snpDistDF3R.seq2[z]) for z in snpDistDF3R.index]

            Graph2.add_edges_from(seqPair1R)
            networkClusters2 = tuple(nx.connected_components(Graph2))

            for r,k in enumerate(networkClusters2):
                for s in k:
                    FirstLevelLin = str(i+1)
                    SecondLevelLin = str(p+1)
                    ThirdLevelLin = str(r+1)
                    clusterPrefix = cmdValues["prefix"]
                    seqLabel = s
                    separator = cmdValues["separator"]

                    if len(networkClusters1)==1:
                        fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}\n")
                    else:
                        if len(networkClusters2)==1:
                            fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}.{SecondLevelLin}\n")
                        else:
                            fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}.{SecondLevelLin}.{ThirdLevelLin}\n")
                        
    fhandle.close();


def lineages4Levels(cmdValues):
    np.random.rand(1)

    snpDistDF = pd.read_csv(cmdValues["data"],delimiter="\t",header=0)

    Graph = nx.Graph()
    Graph.add_nodes_from(snpDistDF.seq2.to_list())
    Graph.add_nodes_from(snpDistDF.seq1.to_list())

    snpDistDF1 = snpDistDF[snpDistDF.pid >= cmdValues["thresholds"][0]]
    seqPair = [(snpDistDF1.seq1[i],snpDistDF1.seq2[i]) for i in snpDistDF1.index]

    Graph.add_edges_from(seqPair)
    networkClusters = tuple(nx.connected_components(Graph))

    fhandle = open(f"{cmdValues['output']}","w")

    for i,j in enumerate(networkClusters):
        snpDistDF2 = snpDistDF[snpDistDF.seq1.isin(j) & snpDistDF.seq2.isin(j)]
        snpDistDF3 = snpDistDF2[snpDistDF2.pid >= cmdValues["thresholds"][1]]

        Graph1 = nx.Graph()
        Graph1.add_nodes_from(j)

        seqPair1 = [(snpDistDF3.seq1[z],snpDistDF3.seq2[z]) for z in snpDistDF3.index]

        Graph1.add_edges_from(seqPair1)
        networkClusters1 = tuple(nx.connected_components(Graph1))

        for p,q in enumerate(networkClusters1):
            snpDistDF2R = snpDistDF[snpDistDF.seq1.isin(q) & snpDistDF.seq2.isin(q)]
            snpDistDF3R = snpDistDF2R[snpDistDF2R.pid >= cmdValues["thresholds"][2]]

            Graph2 = nx.Graph()
            Graph2.add_nodes_from(q)

            seqPair1R = [(snpDistDF3R.seq1[z],snpDistDF3R.seq2[z]) for z in snpDistDF3R.index]

            Graph2.add_edges_from(seqPair1R)
            networkClusters2 = tuple(nx.connected_components(Graph2))

            for x,y in enumerate(networkClusters2):

                snpDistDF2RR = snpDistDF[snpDistDF.seq1.isin(y) & snpDistDF.seq2.isin(y)]
                snpDistDF3RR = snpDistDF2RR[snpDistDF2RR.pid >= cmdValues["thresholds"][3]]

                Graph3 = nx.Graph()
                Graph3.add_nodes_from(y)

                seqPair1RR = [(snpDistDF3RR.seq1[z],snpDistDF3RR.seq2[z]) for z in snpDistDF3RR.index]

                Graph3.add_edges_from(seqPair1RR)
                networkClusters3 = tuple(nx.connected_components(Graph3))

                for r,k in enumerate(networkClusters3):

                    for s in k:
                        FirstLevelLin = str(i+1)
                        SecondLevelLin = str(p+1)
                        ThirdLevelLin = str(x+1)
                        FourthLevelLin = str(r+1)
                        clusterPrefix = cmdValues["prefix"]
                        seqLabel = s
                        separator = cmdValues["separator"]

                        if len(networkClusters1)==1:
                            fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}\n")
                        else:
                            if len(networkClusters2)==1:
                                fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}.{SecondLevelLin}\n")
                            else:
                                if len(networkClusters3)==1:
                                    fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}.{SecondLevelLin}.{ThirdLevelLin}\n")
                                else:
                                    fhandle.write(f"{seqLabel}\t{clusterPrefix}{separator}{FirstLevelLin}.{SecondLevelLin}.{ThirdLevelLin}.{FourthLevelLin}\n")

    fhandle.close();


def cluster(cmdValues):
    for item in cmdValues['thresholds']:
        if float(item)<0 or float(item)>100:
            print("Clustering thresholds should be between 0 and 100")
            sys.exit()

        else:
            pass

    for itemPos,item in enumerate(cmdValues['thresholds']):
        if itemPos == 0:
            pass
        else:
            if cmdValues['thresholds'][itemPos-1]<cmdValues['thresholds'][itemPos]:
                pass
            else:
                print("Clustering thresholds should be provided in increasing order... exiting...")
                sys.exit()

    if len(cmdValues['thresholds'])==1:
        lineages1Levels(cmdValues)
    elif len(cmdValues['thresholds'])==2:
        lineages2Levels(cmdValues)
    elif len(cmdValues['thresholds'])==3:
        lineages3Levels(cmdValues)
    elif len(cmdValues['thresholds'])==4:
        lineages4Levels(cmdValues)
    else:
        print("Maximum number of nested should be 4... exiting...")
        sys.exit()


def order(cmdValues):
    with open(cmdValues["output"],"w") as outputFile:
        dataMatrixDF = loadDataTable(cmdValues)

        itemNames = [re.split("\t|,|;",str(item).strip())[0] for item in open(cmdValues["rownames"],"r")]

        bar = ChargingBar('Saving ordered table', max=len(itemNames))
        for itemNum,itemRowName in enumerate(itemNames):
            if itemNum == 0:
                headerName = list(dataMatrixDF.keys())[0]
                outputFile.write(f"{headerName}\t{"\t".join([str(i) for i in dataMatrixDF[headerName]])}\n")
            else:
                pass

            if itemRowName in list(dataMatrixDF.keys()):
                outputFile.write(f"{itemRowName}\t{"\t".join([str(i) for i in dataMatrixDF[itemRowName]])}\n")

            else:
                outputFile.write(f"{itemRowName}\n")

            bar.next()
        bar.finish()

def sample(cmdValues):
    np.random.seed(1)

    dataMatrixDF = loadDataTable(cmdValues)

    rows = list(dataMatrixDF.keys())
    columns = np.arange(len(dataMatrixDF[rows[0]])).tolist()
    
    selectedColumns = np.random.choice(columns,int(np.ceil(float(cmdValues["colfreq"])*len(columns))),replace=False).tolist()
    selectedRows = [rows[0]]
    selectedRows.extend(np.random.choice(rows[1:],int(np.ceil(float(cmdValues["rowfreq"])*len(rows[1:]))),replace=False).tolist())

    with open(cmdValues["output"],"w") as outputFile:

        bar = ChargingBar('Saving table', max=len(selectedRows)*len(selectedColumns))
        for i,j in enumerate(selectedRows):
            tmpStr = []

            for r,s in enumerate(selectedColumns):
                if r == 0:
                    tmpStr.append(str(j))
                    tmpStr.append(str(dataMatrixDF[j][r]))
                else:
                    tmpStr.append(str(dataMatrixDF[j][r]))

                bar.next()
            outputFile.write(f"{"\t".join(tmpStr)}\n")
        
        bar.finish()


def rename(cmdValues):
    with open(cmdValues["output"],"w") as outputFile:
        dataMatrixDF = loadDataTable(cmdValues)

        columnIndexValues = []

        if not os.path.exists(cmdValues["newdata"]):
            print(f"File {cmdValues["newdata"]} does not exist... exiting...")
            sys.exit()
        else:
            pass

        if not isinstance(cmdValues["columns"],str):
            columnIndexValues = [int(i) for i in [re.split("\t|,|;",str(item).strip()) for item in cmdValues["columns"]][0]]

        else:
            columnIndexValues = [int(i) for i in cmdValues["columns"]]

        newRowNameValuesDF = {}

        tmpCmd="wc -l "+str(cmdValues["newdata"])+" | awk '{print $1}'"
        dataFileLines=subprocess.Popen(tmpCmd,shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        dataLines = int(dataFileLines.communicate()[0])
        dataFileLines.kill()

        bar = ChargingBar('Loading file with new names', max=int(dataLines))
        with open(cmdValues["newdata"],"r") as newData:
            for rowData in newData:
                newRowNameValues = []

                if len(rowData.strip().split("\t")) > 0:
                    if max(columnIndexValues) < 0 or max(columnIndexValues) > len(rowData.strip().split("\t")):
                        print(f"Values for --column/-c should be >0 and <number of columns in the new data table... exiting...")
                        sys.exit()
                    else:   
                        pass
                else:
                    pass

                for columnIndex in columnIndexValues:
                    newRowNameValues.append(rowData.strip().split("\t")[columnIndex-1])

                newRowNameValuesDF[rowData.strip().split("\t")[0]] = "_".join(newRowNameValues)

            bar.next()
        bar.finish()

        bar = ChargingBar('Saving renamed table', max=len(dataMatrixDF.keys()))
        for rowIndex,rowName in enumerate(dataMatrixDF.keys()):

            if rowIndex == 0:
                outputFile.write(f"{rowName}\t{"\t".join([str(i) for i in dataMatrixDF[rowName]])}\n")
            else:
                if rowName in newRowNameValuesDF.keys():
                    outputFile.write(f"{newRowNameValuesDF[rowName]}\t{"\t".join([str(i) for i in dataMatrixDF[rowName]])}\n")

                else:
                    outputFile.write(f"{rowName}\t{"\t".join([str(i) for i in dataMatrixDF[rowName]])}\n")

            bar.next()
        bar.finish()


def filter(cmdValues):
    if cmdValues["minCutoff"]>=0 and cmdValues["maxCutoff"]<=100:
        if cmdValues["minCutoff"] > cmdValues["maxCutoff"]:
            print(f"--->>>Value of --min/-a should be less than or equal to value of --max/-b... exiting...")
            sys.exit()
        else:
            pass
        
        dataMatrixDF = loadDataTable(cmdValues)

        tmpMatrixData = {}

        with open(cmdValues["output"],"w") as outputFile:
            columnValues = []

            bar = ChargingBar('Filtering table', max=len(dataMatrixDF[list(dataMatrixDF.keys())[0]]))
            for itemPos,itemValue in enumerate(dataMatrixDF[list(dataMatrixDF.keys())[0]]):
                itemMatrixColValues = []

                for tmpPos,itemMatrixName in enumerate(dataMatrixDF.keys()):
                    if tmpPos != 0:
                        itemMatrixColValues.append(dataMatrixDF[itemMatrixName][itemPos])
                    else:
                        pass
                        
                rowMatrixSum = sum(itemMatrixColValues)

                if rowMatrixSum >= (cmdValues["minCutoff"]/100)*len(itemMatrixColValues) and \
                    rowMatrixSum <= (cmdValues["maxCutoff"]/100)*len(itemMatrixColValues):

                    columnValues.append(itemPos)
                else:
                    pass

                bar.next()
            bar.finish()

            bar = ChargingBar('Saving filtered table', max=len(list(dataMatrixDF.keys())))
            for itemMatrixName in dataMatrixDF.keys():
                tmpColData = []

                for colPosition in columnValues:
                    tmpColData.append(str(dataMatrixDF[itemMatrixName][colPosition]))

                bar.next()
                outputFile.write(f"{str(itemMatrixName)}\t{"\t".join(tmpColData)}\n")
            bar.finish()

    else:
        if cmdValues["minCutoff"] > cmdValues["maxCutoff"]:
            print(f"--->>>Value of --min/-a should be less than or equal to value of --max/-b... exiting...")
            sys.exit()
        else:
            print(f"--->>>Value of --min/-a snd --max/-b values should be between 0 and 100... exiting...")
            sys.exit()


def makeMapPED(i,j,k,l):
    if i == 0:
        pass
    else:
        tmpStr = ""

        for m,n in enumerate(k):
            if m == 0:
                tmpStr = f"{str(j)}\t{str(j)}\t0\t0\t1\t{str(l)}"

            else:
                tmpN = 0
                if(int(n) == 0):
                    tmpN = 2
                else:
                    tmpN = 1

                tmpStr = f"{tmpStr}\t{str(tmpN)} {str(tmpN)}"

        print('.', end='', flush=True)

        return(tmpStr)


def gwas(cmdValues):
    checkFlags = [1 for flag in [cmdValues["pyseer"],cmdValues["plink"],cmdValues["gemma"],cmdValues["fastlmm"]] if flag == True]

    if sum(checkFlags) > 1:
        print(f"Only one these flags should be specified: --pyseer/-s, --plink/-l, --gemma/g, and --fastlmm/-m... exiting...")
        sys.exit()
    else:
        if sum(checkFlags) == 0:
            cmdValues["plink"] = True

        else:
            pass

    dataMatrixDF = loadDataTable(cmdValues)

    phenoData={}

    if os.path.exists(cmdValues["phenotype"]):
        phenoDataTmp = [str(item).strip() for item in open(cmdValues["phenotype"],"r")]

        phenoData = {}

        bar = ChargingBar('Processing phenotype data', max=len(phenoDataTmp))
        for m in phenoDataTmp:
            tmpR = re.split("\t|,|;",str(m).strip()) 
            phenoData[tmpR[0]]=tmpR[1]

            bar.next()
        bar.finish()

    else:
        print(f"Phenotype file {cmdValues["phenotype"]} not found... exiting...")
        sys.exit()

    if cmdValues["plink"]:
        with open(f"{Path(cmdValues["output"]).stem}.plink.map","w") as mapFile:
            with open(f"{Path(cmdValues["output"]).stem}.plink.ped","w") as pedFile:

                tmpItems = [(i,j) for i,j in enumerate(dataMatrixDF.keys()) if i>0]

                with multiprocessing.Pool(processes=cmdValues["threads"]) as pool:
                    args = [(i,j,dataMatrixDF[j],phenoData[j]) for i,j in tmpItems]
                    resultsMapPED = pool.starmap(makeMapPED, args)

                bar = ChargingBar("Saving PLINK PED file", max=len(list(dataMatrixDF.keys())))
                for eachResult in resultsMapPED:
                    pedFile.write(f"{eachResult}\n")

                    bar.next()
                bar.finish()

                bar = ChargingBar("Saving PLINK MAP file", max=len(dataMatrixDF[list(dataMatrixDF.keys())[0]]))
                for k,l in enumerate(dataMatrixDF[list(dataMatrixDF.keys())[0]]):
                    if k != 0:
                        mapFile.write(f"26\t{str(l)}\t0\t{str(k)}\n")

                    else:
                        pass

                    bar.next()
                bar.finish()

    elif cmdValues["fastlmm"]:
        with open(f"{Path(cmdValues["output"]).stem}.fastlmm.map","w") as mapFile:
            with open(f"{Path(cmdValues["output"]).stem}.fastlmm.ped","w") as pedFile:

                tmpItems = [(i,j) for i,j in enumerate(dataMatrixDF.keys()) if i>0]

                with multiprocessing.Pool(processes=cmdValues["threads"]) as pool:
                    args = [(i,j,dataMatrixDF[j],phenoData[j]) for i,j in tmpItems]
                    resultsMapPED = pool.starmap(makeMapPED, args)

                bar = ChargingBar("Saving FASTLMM PED file", max=len(list(dataMatrixDF.keys())))
                for eachResult in resultsMapPED:
                    pedFile.write(f"{eachResult}\n")

                    bar.next()
                bar.finish()


                bar = ChargingBar("Saving FASTLMM MAP file", max=len(dataMatrixDF[list(dataMatrixDF.keys())[0]]))
                for k,l in enumerate(dataMatrixDF[list(dataMatrixDF.keys())[0]]):
                    if k != 0:
                        mapFile.write(f"26\t{str(l)}\t0\t{str(k)}\n")

                    else:
                        pass

                    bar.next()
                bar.finish()

    elif cmdValues["gemma"]:
        with open(f"{Path(cmdValues["output"]).stem}.plink.map","w") as mapFile:
            with open(f"{Path(cmdValues["output"]).stem}.plink.ped","w") as pedFile:

                tmpItems = [(i,j) for i,j in enumerate(dataMatrixDF.keys()) if i>0]

                with multiprocessing.Pool(processes=cmdValues["threads"]) as pool:
                    args = [(i,j,dataMatrixDF[j],phenoData[j]) for i,j in tmpItems]
                    resultsMapPED = pool.starmap(makeMapPED, args)

                bar = ChargingBar("Saving PLINK PED file", max=len(list(dataMatrixDF.keys())))
                for eachResult in resultsMapPED:
                    pedFile.write(f"{eachResult}\n")

                    bar.next()
                bar.finish()

                bar = ChargingBar("Saving PLINK MAP file", max=len(dataMatrixDF[list(dataMatrixDF.keys())[0]]))
                for k,l in enumerate(dataMatrixDF[list(dataMatrixDF.keys())[0]]):
                    if k != 0:
                        mapFile.write(f"26\t{str(l)}\t0\t{str(k)}\n")

                    else:
                        pass

                    bar.next()
                bar.finish()

        if not shutil.which("plink"):
            print("Install 'plink' software dependency and try again... exiting...")
            sys.exit()
        else:
            pass

        bar = ChargingBar("Generating GEMMA binary PED files", max=1)
        plinkCmd = f"plink -file {f"{Path(cmdValues["output"]).stem}.plink"} -make-bed -out {f"{Path(cmdValues["output"]).stem}.gemma"} --threads {cmdValues["threads"]}"

        subprocess.call(plinkCmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
        bar.next()
        bar.finish()

        os.remove(f"{Path(cmdValues["output"]).stem}.plink.map")
        os.remove(f"{Path(cmdValues["output"]).stem}.plink.ped")
        os.remove(f"{Path(cmdValues["output"]).stem}.gemma.log")

    else:
        with open(f"{Path(cmdValues["output"]).stem}.pyseer.tsv","w") as pyseerTab:
            with open(f"{Path(cmdValues["output"]).stem}.pyseer.pheno.tsv","w") as pyseerPheno:

                transpose(cmdValues)

                transposedFileName = f"{Path(cmdValues["output"]).stem}.transposed.tsv"
                pyseerFileName = f"{Path(cmdValues["output"]).stem}.pyseer.tsv"

                shutil.move(transposedFileName,pyseerFileName)

                bar = ChargingBar("Generating PYSEER files", max=len(phenoData.keys()))
                for itemPos,item in enumerate(phenoData.keys()):
                    if itemPos == 0:
                        pyseerPheno.write("samples\tphenotype\n")
                    else:
                        pyseerPheno.write(f"{item}\t{phenoData[item]}\n")

                    bar.next()
                bar.finish()


def queryKmers(seqNum,seqFile,kmerFile,cmdValues):
    seqFileName = Path(os.path.basename(seqFile)).stem
    tmpOutFileA = f"bifrost.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

    tmpCmd1 = f"Bifrost build -t 3 -k 31 -i -d -r {seqFile} -o {tmpOutFileA}.graph" 
    subprocess.call(tmpCmd1,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    if cmdValues["approximate"]:
        tmpCmd2 = f"Bifrost query -a -t 5 -e 1 -q {kmerFile} -g {tmpOutFileA}.graph.gfa.gz -o {tmpOutFileA}"
        subprocess.call(tmpCmd2,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    else:
        tmpCmd2 = f"Bifrost query -t 5 -e 1 -q {kmerFile} -g {tmpOutFileA}.graph.gfa.gz -o {tmpOutFileA}" 
        subprocess.call(tmpCmd2,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    tmpCmd3 = f"echo {seqFileName} > {tmpOutFileA}.out.tsv"
    subprocess.call(tmpCmd3,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    tmpCmd4 = f"grep -v query_name {tmpOutFileA}.tsv | sort -k1 -n > {tmpOutFileA}.out.tmp.tsv"
    subprocess.call(tmpCmd4,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    tmpCmd5 = "awk '{print $2}' "+f"{tmpOutFileA}.out.tmp.tsv >> {tmpOutFileA}.out.tsv"
    subprocess.call(tmpCmd5,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    tmpCmd6 = f"echo KmerID > {tmpOutFileA}.head.tsv"
    subprocess.call(tmpCmd6,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    tmpCmd7 = "awk '{print $1}' "+f"{tmpOutFileA}.out.tmp.tsv | sed \"s:query_name:KmerID:g\" >> {tmpOutFileA}.head.tsv";
    subprocess.call(tmpCmd7,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    tmpCmd8 = f"datamash transpose < {tmpOutFileA}.head.tsv > {tmpOutFileA}.header.tsv"
    subprocess.call(tmpCmd8,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    tmpCmd9 = f"datamash transpose < {tmpOutFileA}.out.tsv > {tmpOutFileA}.tmp.tsv"
    subprocess.call(tmpCmd9,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

    kmerData = {seqFileName: {"kmerFile": [str(i).strip() for i in open(f"{tmpOutFileA}.tmp.tsv","r")][0],"kmerHeaderFile": [str(j).strip() for j in open(f"{tmpOutFileA}.header.tsv","r")][0] }}

    if os.path.exists(f"{tmpOutFileA}.graph.bfi"):
        os.remove(f"{tmpOutFileA}.graph.bfi")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.graph.gfa.gz"):
        os.remove(f"{tmpOutFileA}.graph.gfa.gz")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.out.tsv"):
        os.remove(f"{tmpOutFileA}.out.tsv")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.tsv"):
        os.remove(f"{tmpOutFileA}.tsv")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.out.tmp.tsv"):
        os.remove(f"{tmpOutFileA}.out.tmp.tsv")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.out.tsv"):
        os.remove(f"{tmpOutFileA}.out.tsv")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.head.tsv"):
        os.remove(f"{tmpOutFileA}.head.tsv")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.header.tsv"):
        os.remove(f"{tmpOutFileA}.header.tsv")
    else:
        pass

    if os.path.exists(f"{tmpOutFileA}.tmp.tsv"):
        os.remove(f"{tmpOutFileA}.tmp.tsv")
    else:
        pass

    return(kmerData)


def table(cmdValues):
    if not os.path.exists(cmdValues["data"]):
        print(f"File {cmdValues["data"]} not found... exiting...")
        sys.exit()
    else:
        pass

    if not shutil.which("awk"):
        print("Install awk and try again... exiting...")
        sys.exit()
    else:
        pass

    if not shutil.which("dsk") and not shutil.which("dsk2ascii"):
        print("Install DSK and try again... exiting...")
        sys.exit()
    else:
        pass

    if not shutil.which("Bifrost"):
        print("Install Bifrost and try again... exiting...")
        sys.exit()
    else:
        pass

    if not shutil.which("datamash"):
        print("Install datamash and try again... exiting...")
        sys.exit()
    else:
        pass

    kmersFile = ""

    if cmdValues["kmers"] != "":
        kmersFile = cmdValues["kmers"]
    else:
        kmersFile = f"{Path(cmdValues["output"]).stem}.fasta"

    with open(cmdValues["data"],"r") as seqFile:

        tmpCmd="wc -l "+str(cmdValues["data"])+" | awk '{print $1}'"
        dataFileLines = subprocess.Popen(tmpCmd,shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        dataLines = int(dataFileLines.communicate()[0])
        dataFileLines.kill()

        bar = ChargingBar("Generating k-mers (DSK)", max=int(dataLines))
        for fileNameValue in seqFile:
            fileName = fileNameValue.strip()

            if os.path.exists(fileName):
                if os.path.isfile(fileName) and os.path.getsize(fileName) > 0:
                    pass
                else:
                    print(f"File {fileName} in {cmdValues["data"]} is empty... exiting...")
                    sys.exit()
            else:
                print(f"File {fileName} in {cmdValues["data"]} does not exist... exiting...")
                sys.exit()

            bar.next()
        bar.finish()

    if cmdValues["dsk"] and not cmdValues["bifrost"]:

        if cmdValues["kmers"] == "":
            bar = ChargingBar("Generating k-mers (DSK)", max=2)
            dskCmdTmpFile = f"KMERS.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

            dsk1Cmd = f"dsk -kmer-size {cmdValues["length"]} -abundance-min {cmdValues["minabundance"]} -max-memory {cmdValues["maxmemory"]} -file {cmdValues["data"]} -out {dskCmdTmpFile} -nb-cores {cmdValues["threads"]}"
            subprocess.call(dsk1Cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
            bar.next()      

            tmpKmersFileR = f"KMERS.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

            dsk2Cmd = f"dsk2ascii -file {dskCmdTmpFile}.h5 -out {tmpKmersFileR} -nb-cores {cmdValues["threads"]}"
            subprocess.call(dsk2Cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
            bar.next()

            awkCmd = "awk '{print \">\"NR\"\\n\"$1}' "+str(tmpKmersFileR)+" > "+str(kmersFile)
            subprocess.call(awkCmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
            bar.next()
            bar.finish()

            if os.path.exists(f"{dskCmdTmpFile}.h5"):
                os.remove(f"{dskCmdTmpFile}.h5")
            else:
                pass

            if os.path.exists(f"{dskCmdTmpFile}"):
                os.remove(f"{dskCmdTmpFile}")
            else:
                pass

            if  os.path.exists(f"{dskCmdTmpFile}.txt"):
                os.remove(f"{dskCmdTmpFile}.txt")
            else:
                pass

            if  os.path.exists(f"{tmpKmersFileR}.txt"):
                os.remove(f"{tmpKmersFileR}.txt")
            else:
                pass

            if  os.path.exists(f"{tmpKmersFileR}"):
                os.remove(f"{tmpKmersFileR}")
            else:
                pass

        else:
            pass

    else:
        if cmdValues["kmers"] != "":
            tmpKmersFile = f"KMERS.{str(datetime.datetime.now()).replace(" ","").replace("-",".").replace(":",".")}"

            bar = ChargingBar("Generating unitigs (Bifrost)", max=1)
            bifrostCmd = f"Bifrost build -t {cmdValues["threads"]} -k {cmdValues["length"]} -i -d -r {cmdValues["data"]} -o {tmpKmersFile} -v -a -f -n"

            subprocess.call(bifrostCmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
            bar.next()
            bar.finish()
        
            with open(kmersFile,"w") as kmerData:
                for seqNum,seq in enumerate(SeqIO.parse(f"{tmpKmersFile}.fasta","fasta")):
                    kmerData.write(f">{(seqNum+1)}\n{seq.seq}\n")
   
                if os.path.exists(f"{tmpKmersFile}.bfi"):
                    os.remove(f"{tmpKmersFile}.bfi")
                else:
                    pass

                if os.path.exists(f"{tmpKmersFile}.fasta"):
                    os.remove(f"{tmpKmersFile}.fasta")
                else:
                    pass

        else:
            pass

    tmpCmd="wc -l "+str(cmdValues["data"])+" | awk '{print $1}'"
    dataFileLines = subprocess.Popen(tmpCmd,shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    dataLines = int(dataFileLines.communicate()[0])
    dataFileLines.kill()

    dataMatrixDF = []

    bar = ChargingBar("Processing sequence files", max=int(dataLines))
    for seqFile in open(cmdValues["data"],"r"):
        dataMatrixDF.append(str(seqFile).strip().split("\t")[0])

        bar.next()
    bar.finish()

    seqFileDF = [(pos,item) for pos,item in enumerate(dataMatrixDF)]

    print("Query k-mers against sequences... may take longer... be patient...")

    with multiprocessing.Pool(processes=cmdValues["threads"]) as pool:
        args = [(seqFileNum,seqFileName,kmersFile,cmdValues) for seqFileNum,seqFileName in seqFileDF]
        kmerResults = pool.starmap(queryKmers, args)

    with open(f"{Path(cmdValues["output"]).stem}.table.tsv","w") as kmerData:
        kmerFiles = []
        kmerHeaderFiles = []

        bar = ChargingBar("Writing k-mers table", max=len(kmerResults))
        for eachResultPos,eachResult in enumerate(kmerResults):
            seqName = list(eachResult.keys())[0]

            kmerFiles.append(eachResult[seqName])
            kmerHeaderFiles.append(eachResult[seqName])

            if eachResultPos == 0:
                kmerData.write(f"{eachResult[seqName]["kmerHeaderFile"]}\n")
                kmerData.write(f"{eachResult[seqName]["kmerFile"]}\n")
            else:
                kmerData.write(f"{eachResult[seqName]["kmerFile"]}\n")

            bar.next()
        bar.finish()


def loadDataTable(cmdValues):
    if not (os.path.exists(cmdValues["data"]) and os.path.isfile(cmdValues["data"])):
        print(f"--->>>Input file {cmdValues["data"]} does not exist... exiting...")
        sys.exit()

    else:
        dataMatrixDF = {}

        sepFormat = "\t" if cmdValues["format"] == "tsv" else ","

        if not (os.path.exists(cmdValues["data"]) and os.path.isfile(cmdValues["data"])):
            print(f"--->>>Input file {cmdValues["data"]} does not exist... exiting...")
            sys.exit()
        else:
            with open(cmdValues["data"],"r") as dataFile:

                if not shutil.which("awk"):
                    print("Install awk and try again... exiting...")
                    sys.exit()
                else:
                    pass

                tmpCmd="wc -l "+str(cmdValues["data"])+" | awk '{print $1}'"
                dataFileLines = subprocess.Popen(tmpCmd,shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
                dataLines = int(dataFileLines.communicate()[0])
                dataFileLines.kill()

                bar = ChargingBar('Loading data table', max=int(dataLines))
                for dataPos,dataRow in enumerate(dataFile):
                    if len(dataRow.strip().split(sepFormat)) == 0:
                        pass
                    else:
                        if dataPos == 0:
                            dataRowValues = dataRow.strip().split(sepFormat)
                            dataRowValuesName = dataRowValues[0]
                            dataRowValuesCols = [item for item in dataRowValues[1:]]
                            dataMatrixDF[dataRowValuesName] = dataRowValuesCols
                        else:
                            dataRowValues = dataRow.strip().split(sepFormat)
                            dataRowValuesName = dataRowValues[0]
                            dataRowValuesCols = [int(item) for item in dataRowValues[1:]]
                            dataMatrixDF[dataRowValuesName] = dataRowValuesCols

                            if set(dataRowValuesCols) == {0} or set(dataRowValuesCols) == {1} or set(dataRowValuesCols) == {0,1}:
                                pass
                            else:
                                print(f"--->>>Row {dataPos} in {cmdValues["data"]} contains unexpected other values besides 0 and 1... exiting...")
                                sys.exit()

                    bar.next()
                bar.finish()

        return(dataMatrixDF)


def citation():
    print(f"Chrispin Chaguza. 2026. matrixToSimilarity: calculating percent similarity between pairs of items (rows) \
given tabular text file. https://github.com/ChrispinChaguza/matrixToSimilarity")


def printVersion():
    print(f"tabkit {__version__}")


def tabtkMain():

    optionsCmd = argparse.ArgumentParser(sys.argv[0],
                 usage=argparse.SUPPRESS,
                 description='tabtk: A toolkit for sequence similarity analysis based on tabular data',
                 prefix_chars='-',
                 add_help=True,
                 epilog='Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2026')

    optionsCommand = optionsCmd.add_subparsers(dest="command")

    similarityOptions = optionsCommand.add_parser("similarity")

    similarityOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file (input should be file of fasta files when used with --mash/-m and --fastani/-n options or alignment file when used with --aln/-a option)')
    similarityOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    similarityOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.similarity.tsv",
                        help='Output percent similarity')
    similarityOptions.add_argument('--length','-l',action='store',required=False,nargs=1,
                        metavar='length',dest='length',default=31,type=int,
                        help='K-mer length for --mash/-m or --fastani/-n (default: 31)')
    similarityOptions.add_argument('--size','-s',action='store',required=False,nargs=1,
                        metavar='size',dest='size',default=1000000,type=int,
                        help='Sketch size for MASH; use with --mash (default: 1000000)')
    similarityOptions.add_argument('--mash','-m',action='store_true',default=False,
                        dest='mash',help='Calculate similarity using MASH')
    similarityOptions.add_argument('--fastani','-n',action='store_true',default=False,
                        dest='fastani',help='Calculate similarity using fastANI')
    similarityOptions.add_argument('--blast','-b',action='store_true',default=False,
                        dest='blast',help='Calculate similarity using BLAST')
    similarityOptions.add_argument('--aln','-a',action='store_true',default=False,
                        dest='aln',help='Calculate similarity based on a given sequence alignment') 
    similarityOptions.add_argument('--fraction','-c',action='store',required=False,nargs=1,
                        metavar='fraction',dest='fraction',default=0.70,type=float,
                        help='Minimum overlap fraction with shorter sequence for ANI calculation; use with --fastani (default: 0.70)')
    similarityOptions.add_argument('--fragment','-r',action='store',required=False,nargs=1,
                        metavar='fragment',dest='fragment',default=500,type=int,
                        help='Fragment length for calculating ANI; use with --fastani (default: 500)')
    similarityOptions.add_argument('--threads','-t',action='store',required=False,nargs=1,
                        metavar='threads',dest='threads',default=5,type=int,
                        help='Number of threads (default: 5)')


    orderOptions = optionsCommand.add_parser("order")

    orderOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file')
    orderOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    orderOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.ordered.tsv",
                        help='Output ordered data file')
    orderOptions.add_argument('--alphabetic','-a',action='store_true',required=False,
                        dest='alphabetic',default=False,
                        help='Show similarity values as percentage instead of fraction')
    orderOptions.add_argument('--rownames','-r',action='store',required=False,nargs=1,
                        metavar='rownames',dest='rownames',
                        help='Order by item names in a given text file (one name per row)')
    orderOptions.add_argument('--threads','-t',action='store',required=False,nargs=1,
                        metavar='threads',dest='threads',default=5,type=int,
                        help='Number of threads (default: 5)')

    renameOptions = optionsCommand.add_parser("rename")

    renameOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file')
    renameOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    renameOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.renamed.tsv",
                        help='Output ordered data file')
    renameOptions.add_argument('--columns','-c',action='store',required=True,nargs="*",
                        metavar='columns',dest='columns',
                        help='Selected column numbers (e.g., --columns 1 2 3 4)')
    renameOptions.add_argument('--newdata','-n',action='store',required=True,nargs=1,
                        metavar='newdata',dest='newdata',
                        help='Order by item names in a given text file (one name per row)')


    filterOptions = optionsCommand.add_parser("filter")

    filterOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file')
    filterOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    filterOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.filter.tsv",
                        help='Output percent similarity')
    filterOptions.add_argument('--min','-a',action='store',required=False,nargs=1,
                        metavar='minCutoff',dest='minCutoff',default=0,type=int,
                        help='Select with frequency ≥minCutoff percent (overrides --columns/-c and --rows/-r)')
    filterOptions.add_argument('--max','-b',action='store',required=False,nargs=1,
                        metavar='maxCutoff',dest='maxCutoff',default=100,type=int,
                        help='Select with frequency ≤maxCutoff percent (overrides --columns/-c and --rows/-r)')
    filterOptions.add_argument('--columns','-c',action='store',required=False,nargs=1,
                        metavar='columns',dest='columns',
                        help='Select only columns present in a given text file (one column value per row)')
    filterOptions.add_argument('--rows','-r',action='store',required=False,nargs=1,
                        metavar='rows',dest='rows',
                        help='Select only rows present in a given text file (one column value per row)')
    filterOptions.add_argument('--threads','-t',action='store',required=False,nargs=1,
                        metavar='threads',dest='threads',default=5,type=int,
                        help='Number of threads (default: 5)')


    sampleOptions = optionsCommand.add_parser("sample")

    sampleOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file')
    sampleOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    sampleOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.sample.tsv",
                        help='Output file')
    sampleOptions.add_argument('--cf','-c',action='store',required=False,nargs=1,
                        metavar='colfreq',dest='colfreq',default=0,type=float,
                        help='Randomly select ≥colfreq percent of the columns')
    sampleOptions.add_argument('--rf','-r',action='store',required=False,nargs=1,
                        metavar='rowfreq',dest='rowfreq',default=0,type=float,
                        help='Randomly select ≥rowfreq percent of the rows')


    transposeOptions = optionsCommand.add_parser("transpose")

    transposeOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file')
    transposeOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    transposeOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.filter.tsv",
                        help='Output percent similarity')


    clusterOptions = optionsCommand.add_parser("cluster")

    clusterOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file')
    clusterOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    clusterOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.clusters.tsv",
                        help='Output clusters')
    clusterOptions.add_argument('--prefix','-p',action='store',required=False,nargs=1,
                        metavar='prefix',dest='prefix',default="CLS",
                        help='Cluster label prefix')
    clusterOptions.add_argument('--sep','-s',action='store',required=False,nargs=1,
                        metavar='separator',dest='separator',default="-",
                        help='Lineage label separator')
    clusterOptions.add_argument('--thresholds','-t',action='store',required=True,nargs="*",
                        metavar='thresholds',dest='thresholds',
                        help='Cluster or lineage thresholds')


    gwasOptions = optionsCommand.add_parser("gwas")

    gwasOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input tabular text file')
    gwasOptions.add_argument('--format','-f',action='store',required=False,nargs=1,
                        metavar='format',dest='format',default="tsv",
                        choices=['tsv','csv'],
                        help='Data format (default: tsv)')
    gwasOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="out.gwas.tsv",
                        help='Output percent similarity')
    gwasOptions.add_argument('--phenotype','-p',action='store',required=False,nargs=1,
                        metavar='phenotype',dest='phenotype',
                        help='Input phenotype text file (row name,phenotype); no header')
    gwasOptions.add_argument('--pyseer','-s',action='store_true',required=False,
                        dest='pyseer',default=False,
                        help='Output presence and absence data for Pyseer')
    gwasOptions.add_argument('--plink','-l',action='store_true',required=False,
                        dest='plink',default=False,
                        help='Output presence and absence data for Pyseer')
    gwasOptions.add_argument('--fastlmm','-m',action='store_true',required=False,
                        dest='fastlmm',default=False,
                        help='Output presence and absence data for Pyseer')
    gwasOptions.add_argument('--gemma','-g',action='store_true',required=False,
                        dest='gemma',default=False,
                        help='Output presence and absence data for GEMMA')
    gwasOptions.add_argument('--threads','-t',action='store',required=False,nargs=1,
                        metavar='threads',dest='threads',default=5,type=int,
                        help='Number of threads (default: 5)')


    kmerMatrixOptions = optionsCommand.add_parser("table")

    kmerMatrixOptions.add_argument('--data','-d',action='store',required=True,nargs=1,
                        metavar='data',dest='data',
                        help='Input text files containing fasta sequence file names')
    kmerMatrixOptions.add_argument('--output','-o',action='store',required=False,nargs=1,
                        metavar='output',dest='output',default="db.kmers.fasta",
                        help='Output percent similarity')
    kmerMatrixOptions.add_argument('--length','-l',action='store',required=False,nargs=1,
                        metavar='length',dest='length',default=31,type=int,
                        help='K-mer length (default: 31)')
    kmerMatrixOptions.add_argument('--min','-m',action='store',required=False,nargs=1,
                        metavar='minabundance',dest='minabundance',default=1,type=int,
                        help='Minimum k-mer abundance (default: 1)')
    kmerMatrixOptions.add_argument('--maxmem','-x',action='store',required=False,nargs=1,
                        metavar='maxmemory',dest='maxmemory',default=31,type=int,
                        help='Maximum memory (default: 5000)')
    kmerMatrixOptions.add_argument('--dsk','-s',action='store_true',default=True,
                        dest='dsk',help='DSK to identify k-mers')
    kmerMatrixOptions.add_argument('--bifrost','-b',action='store_true',default=False,
                        dest='bifrost',help='Bifrost to identify variable length k-mers or unitigs')
    kmerMatrixOptions.add_argument('--kmers','-k',action='store',required=False,nargs=1,
                        metavar='kmers',dest='kmers',default="",
                        help='Only query input sequence against a given k-mer database')
    kmerMatrixOptions.add_argument('--approx','-a',action='store_true',required=False,
                        dest='approximate',default=False,
                        help='Search inexact k-mers (default: False)')
    kmerMatrixOptions.add_argument('--threads','-t',action='store',required=False,nargs=1,
                        metavar='threads',dest='threads',default=5,type=int,
                        help='Number of threads (default: 5)')


    versionOptions = optionsCommand.add_parser("version")

    versionOptions.add_argument('--version','-v',action='store_false',default=True,
                        dest='version',help='Show software version')


    citationOptions = optionsCommand.add_parser("citation")

    citationOptions.add_argument('--citation','-c',action='store_false',default=True,
                        dest='citation',help='Show software citation information')


    options = optionsCmd.parse_args(args=None if sys.argv[1:] else ['--help'])

    cmdValues = {}

    if options.command == "similarity":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'length': int(options.length[0:][0]) if isinstance(options.length,list) else int(options.length),
                 'fragment': int(options.fragment[0:][0]) if isinstance(options.fragment,list) else int(options.fragment),
                 'fraction': int(options.fraction[0:][0]) if isinstance(options.fraction,list) else int(options.fraction),
                 'size': int(options.size[0:][0]) if isinstance(options.size,list) else int(options.size),
                 'mash': options.mash,
                 'fastani': options.fastani,
                 'blast': options.blast,
                 'aln': options.aln,
                 'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else int(options.threads)
                 }
        
        similarity(cmdValues)

    elif options.command == "filter":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'minCutoff': int(options.minCutoff[0:][0]) if isinstance(options.minCutoff,list) else int(options.minCutoff),
                 'maxCutoff': int(options.maxCutoff[0:][0]) if isinstance(options.maxCutoff,list) else int(options.maxCutoff),
                 'rows': options.rows[0:][0] if isinstance(options.rows,list) else options.rows,
                 'columns': options.columns[0:][0] if isinstance(options.columns,list) else options.columns,
                 'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else int(options.threads)
                 }

        filter(cmdValues)

    elif options.command == "sample":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'colfreq': float(options.colfreq[0:][0]) if isinstance(options.colfreq,list) else float(options.colfreq),
                 'rowfreq': float(options.rowfreq[0:][0]) if isinstance(options.rowfreq,list) else float(options.rowfreq)
                 }

        sample(cmdValues)

    elif options.command == "order":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'rownames': options.rownames[0:][0] if isinstance(options.rownames,list) else options.rownames,
                 'alphabetic': options.alphabetic,
                 'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else int(options.threads)
                 }

        order(cmdValues)

    elif options.command == "rename":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'newdata': options.newdata[0:][0] if isinstance(options.newdata,list) else options.newdata,
                 'columns': options.columns[0:] if isinstance(options.columns,list) else options.columns
                 }

        rename(cmdValues)

    elif options.command == "cluster":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'prefix': options.prefix[0:][0] if isinstance(options.prefix,list) else options.prefix,
                 'separator': options.separator[0:][0] if isinstance(options.separator,list) else options.separator,
                 'thresholds': [float(i) for i in options.thresholds[0:]] if isinstance(options.thresholds,list) else float(options.thresholds)
                 }

        cluster(cmdValues)

    elif options.command == "transpose":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output
                 }

        transpose(cmdValues)

    elif options.command == "gwas":
        cmdValues = {'data': options.data[0:][0],
                 'format': options.format[0:][0] if isinstance(options.format,list) else options.format,
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'phenotype': options.phenotype[0:][0] if isinstance(options.phenotype,list) else options.phenotype,
                 'pyseer': options.pyseer,
                 'fastlmm': options.fastlmm,
                 'plink': options.plink,
                 'gemma': options.gemma,
                 'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else int(options.threads)
                 }

        gwas(cmdValues)

    elif options.command == "table":
        cmdValues = {'data': options.data[0:][0],
                 'output': options.output[0:][0] if isinstance(options.output,list) else options.output,
                 'dsk': options.dsk,
                 'bifrost': options.bifrost,
                 'minabundance': int(options.minabundance[0:][0]) if isinstance(options.minabundance,list) else int(options.minabundance),
                 'maxmemory': int(options.maxmemory[0:][0]) if isinstance(options.maxmemory,list) else int(options.maxmemory),
                 'length': int(options.length[0:][0]) if isinstance(options.length,list) else int(options.length),
                 'kmers': options.kmers[0:][0] if isinstance(options.kmers,list) else options.kmers,
                 'approximate': options.approximate,
                 'threads': int(options.threads[0:][0]) if isinstance(options.threads,list) else int(options.threads)
                 }

        table(cmdValues)

    elif options.command == "version":
        printVersion()

    elif options.command == "citation":
        citation()

    else:
        sys.exit()


if __name__=="__main__":
    tabtkMain()
