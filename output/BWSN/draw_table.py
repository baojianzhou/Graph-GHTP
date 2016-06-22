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

funcNames = ['Method             ','Kulldorff','EMS','EBP','Run Time']

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
            noiseLevel = fileName.split('_')[5].split('.txt')[0]
            data.append([statFuncName, funcValue, runTime, pre, rec, fmeasure, fValues, noiseLevel])
    kulldorff_aver = np.mean([item[1] for item in data if item[0] == 'Kulldorff'])
    kulldorff_aver = "{0:.2f}".format(kulldorff_aver)
    EMS_aver = np.mean([item[1] for item in data if item[0] == 'EMS'])
    print [item[1] for item in data if item[0] == 'EMS']
    EMS_aver = "{0:.2f}".format(EMS_aver)
    EBP_aver = np.mean([item[1] for item in data if item[0] == 'EBP'])
    EBP_aver = "{0:.2f}".format(EBP_aver)
    runTime = np.mean([item[2] for item in data])
    runTime = "{0:.2f}".format(runTime)
    return [methodName,kulldorff_aver.zfill(6),EMS_aver.zfill(6),EBP_aver.zfill(6),runTime.zfill(7)]

def debug():
    eventTree = get_scores('eventTree_BWSN_Result.txt','EventTree')
    graph_IHTP = get_scores('graph_GHTP_BWSN_Result.txt','GraphGHTP')
    print eventTree
    print graph_IHTP 
    
def get_table():
    graph_IHTP = get_scores('graph_GHTP_BWSN_Result.txt','GraphGHTP          ')
    nphgs = get_scores('nphgs_BWSN_Result.txt','NPHGS              ')
    eventTree = get_scores('eventTree_BWSN_Result.txt','EventTree          ')
    depthFirstGraphScan = get_scores('depthFirstGraphScan_BWSN_Result.txt','DepthFirstGraphScan')
    graphLapLacian = get_scores('graphLaplacian_BWSN_Result.txt','GraphLapLacian     ')
    graphIHT = get_scores('graph_IHT_BWSN_Result.txt','GraphIHT           ')
    print funcNames
    print graph_IHTP
    print eventTree
    print nphgs
    print depthFirstGraphScan
    print graphLapLacian
    print graphIHT
    
if __name__ == '__main__':
    #debug()
    print '------------------------------------------------------------'
    get_table()