import os
import sys
import time
import math
import numpy as np
import operator
import collections
from sets import Set

from matplotlib import pyplot
from matplotlib.lines import Line2D
from matplotlib.pyplot import savefig
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import matplotlib.markers as marks

funcNames = ['Method             ','Kulldorff','EMS','    EBP','  Run Time']

def get_scores(fileName,methodName):
    # noise Level 
    data = []
    with open(fileName) as f:
        for eachLine in f.readlines():
            items = [item.rstrip().lstrip() for item in eachLine.split(',')]
            statFuncName = items[1]
            fileName = items[0]
            funcValue = float(items[2])
            runTime = float(items[3])
            pre = float(items[4])
            rec = float(items[5])
            fmeasure = float(items[6])
            fValues = [float(item) for item in items[7].split(' ')]
            data.append([statFuncName, funcValue, runTime, pre, rec, fmeasure, fValues])
    kulldorff_aver = np.mean([item[1] for item in data if item[0] == 'Kulldorff'])
    kulldorff_aver = "{0:.2f}".format(kulldorff_aver)
    EMS_aver = np.mean([item[1] for item in data if item[0] == 'EMS'])
    EMS_aver = "{0:.2f}".format(EMS_aver)
    EBP_aver = np.mean([item[1] for item in data if item[0] == 'EBP'])
    EBP_aver = "{0:.2f}".format(EBP_aver)
    runTime = np.mean([item[2] for item in data])
    runTime = "{0:.2f}".format(runTime)
    return [methodName,kulldorff_aver.zfill(7),EMS_aver.zfill(6),EBP_aver.zfill(6),runTime.zfill(7)]

def debug():
    graph_IHTP = get_scores('graph_GHTP_CitHepPh_Result.txt','GraphGHTP          ')
    eventTree = get_scores('eventTree_CitHepPh_Result.txt','EventTree          ')
    nphgs = get_scores('nphgs_CitHepPh_Result.txt',    'NPHGS              ')
    graphLapLacian = get_scores('graphLaplacian_CitHepPh_Result.txt',    'GraphLapLacian     ')
    depthFirst = get_scores('depthFirstScan_CitHepPh_Result.txt',    'DepthFirstScan     ')
    print funcNames
    print graph_IHTP
    print eventTree
    print depthFirst
    print nphgs
    print graphLapLacian
    
def get_table():
    graph_IHTP = get_scores('graph_GHTP_CitHepPh_Result.txt','GraphGHTP          ')
    nphgs = get_scores('nphgs_CitHepPh_Result.txt','NPHGS              ')
    eventTree = get_scores('eventTree_CitHepPh_Result.txt','EventTree          ')
    edgeLasso = get_scores('edgeLasso_CitHepPh_Result.txt','EdgeLasso          ')
    depthFirstGraphScan = get_scores('depthFirstGraphScan_CitHepPh_Result.txt','DepthFirstGraphScan')
    averFusedLasso_0_25 = get_scores('genFusedLasso_CitHepPh_Result_0.25_.txt','GenFusedLasso      ')
    print funcNames
    print graph_IHTP
    print averFusedLasso_0_25
    print edgeLasso
    print LTSS
    print eventTree
    print additiveGraphScan
    print depthFirstGraphScan
    print nphgs
    
if __name__ == '__main__':
    debug()
    #print '------------------------------------------'
    #get_table()