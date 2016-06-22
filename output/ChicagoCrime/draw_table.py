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

def get_scores(fileName,methodName,date):
    # noise Level 
    data = []
    with open(fileName) as f:
        for eachLine in f.readlines():
            items = [item.rstrip().lstrip() for item in eachLine.split(',')]
            statFuncName = items[1]
            fileName = items[0]
            if fileName.split('_')[0] == date or date == '':
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
    return [methodName,kulldorff_aver.zfill(6),EMS_aver.zfill(6),EBP_aver.zfill(6),runTime.zfill(7)]

def debug():
    graph_IHTP = get_scores('graph_GHTP_CrimeOfChicago_Result.txt','GraphGHTP          ','')
    eventTree = get_scores('eventTree_CrimeOfChicago_Result.txt','EventTree          ','')
    depthFirst = get_scores('depthFirstScan_CrimeOfChicago_Result.txt','DepthFirst         ','')
    nphgs = get_scores('nphgs_CrimeOfChicago_Result.txt','NPHGS              ','')
    print ' Method                        Kulldorff, EMS, EBP    '
    print graph_IHTP
    print eventTree
    print depthFirst
    print nphgs

def debug_average():
    dates = ["2013-07-01", "2013-07-02", "2013-07-03", "2013-07-04",
            "2013-07-05", "2013-07-06", "2013-07-07", "2013-07-08", "2013-07-09", "2013-07-10", "2013-07-11",
            "2013-07-12", "2013-07-13", "2013-07-14", "2013-07-15", "2013-07-16", "2013-07-17", "2013-07-18",
            "2013-07-19", "2013-07-20", "2013-07-21", "2013-07-22", "2013-07-23", "2013-07-24", "2013-07-25",
            "2013-07-26", "2013-07-27", "2013-07-28", "2013-07-29", "2013-07-30", "2013-07-31",]
    eventTree = []
    graph_IHTP = []
    diff = []
    eventTree.append('EventTree          ')
    graph_IHTP.append('GraphGHTP          ')
    diff = ['Difference         ']
    for date in dates:
        graph_IHTP_each = get_scores('graph_GHTP_Traffic_Result.txt','GraphGHTP          ',date)
        eventTree_each = get_scores('eventTree_Traffic_Result.txt','EventTree          ',date)
        eventTree.append(eventTree_each[2])
        graph_IHTP.append(graph_IHTP_each[2])
        diff.append(float(graph_IHTP_each[2]) - float(eventTree_each[2]))
    plt.plot(eventTree[1:])
    plt.plot(graph_IHTP[1:])
    plt.show()
    print eventTree
    print graph_IHTP 
    print diff
    
def get_table():
    graph_IHTP = get_scores('graph_GHTP_CitHepPh_Result.txt','GraphGHTP          ')
    nphgs = get_scores('nphgs_CitHepPh_Result.txt','NPHGS              ')
    eventTree = get_scores('eventTree_CitHepPh_Result.txt','EventTree          ')
    edgeLasso = get_scores('edgeLasso_CitHepPh_Result.txt','EdgeLasso          ')
    LTSS = get_scores('ltss_CitHepPh_Result.txt','LTSS               ')
    depthFirstGraphScan = get_scores('depthFirstGraphScan_CitHepPh_Result.txt','DepthFirstGraphScan')
    additiveGraphScan = get_scores('additiveGraphScan_CitHepPh_Result.txt','AdditiveGraphScan  ')
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
    #debug_average()
    #print '------------------------------------------'
    #get_table()