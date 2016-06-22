import numpy as np
import os
import sys
import time
import math
import collections
import operator
from sets import Set
from matplotlib import pyplot
import matplotlib.pyplot as plt
import matplotlib.markers as marks
from matplotlib.lines import Line2D
from matplotlib.pyplot import savefig

def get_pre_rec_fmeasure(fileName, statFuncID):
    # noise Level 
    x = [ 2, 4, 6, 8, 10]
    data = []
    
    with open(fileName) as f:
        for eachLine in f.readlines():
            items = [item.rstrip().lstrip() for item in eachLine.split(',')]
            fileName = items[0]
            statFuncName = items[1]
            funcValue = float(items[2])
            runTime = float(items[3])
            pre = float(items[4])
            rec = float(items[5])
            fmeasure = float(items[6])
            fValues = [float(item) for item in items[7].split(' ')]
            noiseLevel = fileName.split('_')[5].split('.txt')[0]
            data.append([statFuncName, funcValue, runTime, pre, rec, fmeasure, fValues, noiseLevel])
    precision = [0.0 for item in range(5)]
    precision[0] = np.mean([item[3] for item in data if item[7] == '2' and item[0] == statFuncID])
    precision[1] = np.mean([item[3] for item in data if item[7] == '4' and item[0] == statFuncID])
    precision[2] = np.mean([item[3] for item in data if item[7] == '6' and item[0] == statFuncID])
    precision[3] = np.mean([item[3] for item in data if item[7] == '8' and item[0] == statFuncID])
    print [item[3] for item in data if item[7] == '10' and item[0] == statFuncID]
    precision[4] = np.mean([item[3] for item in data if item[7] == '10' and item[0] == statFuncID])
    recall = [0.0 for item in range(5)]
    recall[0] = np.mean([item[4] for item in data if item[7] == '2' and item[0] == statFuncID])
    recall[1] = np.mean([item[4] for item in data if item[7] == '4' and item[0] == statFuncID])
    recall[2] = np.mean([item[4] for item in data if item[7] == '6' and item[0] == statFuncID])
    recall[3] = np.mean([item[4] for item in data if item[7] == '8' and item[0] == statFuncID])
    recall[4] = np.mean([item[4] for item in data if item[7] == '10' and item[0] == statFuncID])
    fmeasure = [0.0 for item in range(5)]
    fmeasure[0] = np.mean([item[5] for item in data if item[7] == '2' and item[0] == statFuncID])
    fmeasure[1] = np.mean([item[5] for item in data if item[7] == '4' and item[0] == statFuncID])
    fmeasure[2] = np.mean([item[5] for item in data if item[7] == '6' and item[0] == statFuncID])
    fmeasure[3] = np.mean([item[5] for item in data if item[7] == '8' and item[0] == statFuncID])
    fmeasure[4] = np.mean([item[5] for item in data if item[7] == '10' and item[0] == statFuncID])
    print 'precision: ', precision
    print 'recall: ', recall
    print 'fmeasure: ', fmeasure
    return x, precision, recall, fmeasure
    
def draw_precision_recall_fmeasure():
    # input FileList 
    c = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
    markers = ['ro-', 'bs-', 'g^-', 'c+-',  'yD-','m+-', 'k^-']
    statFuncID = 'EBP'
    
    ########################## precision ##################################
    plt.figure(1)
    fig = plt.plot()
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('graph_GHTP_BWSN_Result.txt', statFuncID)
    marker_graph_GHTP, = plt.plot(x, pre, markers[0], markersize=10, label='Graph-GHTP',
             markerfacecolor='white', markeredgecolor=c[0], markeredgewidth=3.0)
    lin_graph_GHTP, = plt.plot(x, pre, linestyle='-', color=c[0], linewidth=2.0, label='Graph-GHTP')
    legend = plt.legend([(marker_graph_GHTP, lin_graph_GHTP)], ['Graph-GHTP'], loc='lower right', shadow=True, prop={'size':14})
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('eventTree_BWSN_Result.txt', statFuncID)
    marker_eventTree, = plt.plot(x, pre, markers[1], markersize=10, label='EventTree',
             markerfacecolor='white', markeredgecolor=c[1], markeredgewidth=3.0)
    lin_eventTree, = plt.plot(x, pre, linestyle='-', color=c[1], linewidth=2.0, label='EventTree')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('nphgs_BWSN_Result.txt', statFuncID)
    marker_nphgs, = plt.plot(x, pre, markers[2], markersize=10, label='NPHGS',
             markerfacecolor='white', markeredgecolor=c[2], markeredgewidth=3.0)
    lin_nphgs, = plt.plot(x, pre, linestyle='-', color=c[2], linewidth=2.0, label='NPHGS')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('depthFirstGraphScan_BWSN_Result.txt', statFuncID)
    marker_depthFirstGraphScan, = plt.plot(x, pre, markers[3], markersize=10, label='DepthFirstGraphScan',
             markerfacecolor='white', markeredgecolor=c[3], markeredgewidth=3.0)
    lin_depthFirstGraphScan, = plt.plot(x, pre, linestyle='-', color=c[3], linewidth=2.0, label='DepthFirstGraphScan')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('graphLaplacian_BWSN_Result.txt', statFuncID)
    marker_graphLaplacian, = plt.plot(x, pre, markers[4], markersize=10, label='GraphLaplacian',
             markerfacecolor='white', markeredgecolor=c[4], markeredgewidth=3.0)
    lin_graphLaplacian, = plt.plot(x, pre, linestyle='-', color=c[4], linewidth=2.0, label='GraphLaplacian')
    #--------------------------------------------------------------------------------------------------------
    
    legend = plt.legend([(marker_graph_GHTP, lin_graph_GHTP),
                         (marker_eventTree, lin_eventTree),
                         (marker_nphgs, lin_nphgs),
                         (marker_depthFirstGraphScan, lin_depthFirstGraphScan),
                         (marker_graphLaplacian, lin_graphLaplacian), ],
                         ['Graph-GHTP', 'EventTree', 'NPHGS', 'DepthFirstGraphScan', 'GraphLaplacian'],
                          loc='lower left', shadow=True, prop={'size':14})
    
    
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(3)
    plt.xlabel(r'Noise Level', fontsize=20)
    plt.ylabel(r'precision', fontsize=20)
    plt.ylim([0.0, 1.05])
    plt.xlim([1.0, 11.0])
    # plt.savefig(os.path.join('../pics/','BWSN_precision_'+statFuncID+'.png'))
    # plt.close()
    # plt.show()
    
    ########################## recall ##################################
    plt.figure(2)
    fig = plt.plot()
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('graph_GHTP_BWSN_Result.txt', statFuncID)
    marker_graph_GHTP, = plt.plot(x, rec, markers[0], markersize=10, label='Graph-GHTP',
             markerfacecolor='white', markeredgecolor=c[0], markeredgewidth=3.0)
    lin_graph_GHTP, = plt.plot(x, rec, linestyle='-', color=c[0], linewidth=2.0, label='Graph-GHTP')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('eventTree_BWSN_Result.txt', statFuncID)
    marker_eventTree, = plt.plot(x, rec, markers[1], markersize=10, label='EventTree',
             markerfacecolor='white', markeredgecolor=c[1], markeredgewidth=3.0)
    lin_eventTree, = plt.plot(x, rec, linestyle='-', color=c[1], linewidth=2.0, label='EventTree')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('nphgs_BWSN_Result.txt', statFuncID)
    marker_nphgs, = plt.plot(x, rec, markers[2], markersize=10, label='NPHGS',
             markerfacecolor='white', markeredgecolor=c[2], markeredgewidth=3.0)
    lin_nphgs, = plt.plot(x, rec, linestyle='-', color=c[2], linewidth=2.0, label='NPHGS')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('depthFirstGraphScan_BWSN_Result.txt', statFuncID)
    marker_depthFirstGraphScan, = plt.plot(x, rec, markers[3], markersize=10, label='DepthFirstGraphScan',
             markerfacecolor='white', markeredgecolor=c[3], markeredgewidth=3.0)
    lin_depthFirstGraphScan, = plt.plot(x, rec, linestyle='-', color=c[3], linewidth=2.0, label='DepthFirstGraphScan')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('graphLaplacian_BWSN_Result.txt', statFuncID)
    marker_graphLaplacian, = plt.plot(x, rec, markers[4], markersize=10, label='GraphLaplacian',
             markerfacecolor='white', markeredgecolor=c[4], markeredgewidth=3.0)
    lin_graphLaplacian, = plt.plot(x, rec, linestyle='-', color=c[4], linewidth=2.0, label='GraphLaplacian')
    #--------------------------------------------------------------------------------------------------------
    legend = plt.legend([(marker_graph_GHTP, lin_graph_GHTP),
                         (marker_eventTree, lin_eventTree),
                         (marker_nphgs, lin_nphgs),
                         (marker_depthFirstGraphScan, lin_depthFirstGraphScan),
                         (marker_graphLaplacian, lin_graphLaplacian), ],
                         ['Graph-GHTP', 'EventTree', 'NPHGS', 'DepthFirstGraphScan', 'GraphLaplacian'],
                         loc='lower left', shadow=True, prop={'size':14})
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(3)
    plt.xlabel(r'Noise Level', fontsize=20)
    plt.ylabel(r'recall', fontsize=20)
    plt.ylim([0.0, 1.05])
    plt.xlim([1.0, 11.0])
    # plt.savefig(os.path.join('../pics/','BWSN_recall_'+statFuncID+'.png'))
    # plt.close()
    # plt.show()
    
    ########################## fmeasure ##################################
    plt.figure(3)
    fig = plt.plot()
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('graph_GHTP_BWSN_Result.txt', statFuncID)
    marker_graph_GHTP, = plt.plot(x, fmeasure, markers[0], markersize=10, label='Graph-GHTP',
             markerfacecolor='white', markeredgecolor=c[0], markeredgewidth=3.0)
    lin_graph_GHTP, = plt.plot(x, fmeasure, linestyle='-', color=c[0], linewidth=2.0, label='Graph-GHTP')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('eventTree_BWSN_Result.txt', statFuncID)
    marker_eventTree, = plt.plot(x, fmeasure, markers[1], markersize=10, label='EventTree',
             markerfacecolor='white', markeredgecolor=c[1], markeredgewidth=3.0)
    lin_eventTree, = plt.plot(x, fmeasure, linestyle='-', color=c[1], linewidth=2.0, label='EventTree')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('nphgs_BWSN_Result.txt', statFuncID)
    marker_nphgs, = plt.plot(x, fmeasure, markers[2], markersize=10, label='NPHGS',
             markerfacecolor='white', markeredgecolor=c[2], markeredgewidth=3.0)
    lin_nphgs, = plt.plot(x, fmeasure, linestyle='-', color=c[2], linewidth=2.0, label='NPHGS')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('depthFirstGraphScan_BWSN_Result.txt', statFuncID)
    marker_depthFirstGraphScan, = plt.plot(x, fmeasure, markers[3], markersize=10, label='DepthFirstGraphScan',
             markerfacecolor='white', markeredgecolor=c[3], markeredgewidth=3.0)
    lin_depthFirstGraphScan, = plt.plot(x, fmeasure, linestyle='-', color=c[3], linewidth=2.0, label='DepthFirstGraphScan')
    #--------------------------------------------------------------------------------------------------------
    x, pre, rec, fmeasure = get_pre_rec_fmeasure('graphLaplacian_BWSN_Result.txt', statFuncID)
    marker_graphLaplacian, = plt.plot(x, fmeasure, markers[4], markersize=10, label='GraphLaplacian',
             markerfacecolor='white', markeredgecolor=c[4], markeredgewidth=3.0)
    lin_graphLaplacian, = plt.plot(x, fmeasure, linestyle='-', color=c[4], linewidth=2.0, label='GraphLaplacian')
    #--------------------------------------------------------------------------------------------------------
    legend = plt.legend([(marker_graph_GHTP, lin_graph_GHTP),
                         (marker_eventTree, lin_eventTree),
                         (marker_nphgs, lin_nphgs),
                         (marker_depthFirstGraphScan, lin_depthFirstGraphScan),
                         (marker_graphLaplacian, lin_graphLaplacian), ],
                         ['Graph-GHTP', 'EventTree', 'NPHGS', 'DepthFirstGraphScan', 'GraphLaplacian'],
                         loc='lower left', shadow=True, prop={'size':14})
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(3)
    plt.xlabel(r'Noise Level', fontsize=20)
    plt.ylabel(r'fmeasure', fontsize=20)
    plt.ylim([0.0, 1.05])
    plt.xlim([1.0, 11.0])
    # plt.savefig(os.path.join('../pics/','BWSN_fmeasure_'+statFuncID+'.png'))
    # plt.close()
    plt.show()
    
if __name__ == '__main__':
    draw_precision_recall_fmeasure()
