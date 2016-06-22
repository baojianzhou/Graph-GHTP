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

def get_funcValues(fileName,statFuncID):
    # noise Level 
    iterNum = 10
    data = []
    with open(fileName) as f:
        for eachLine in f.readlines():
            items = [item.rstrip().lstrip() for item in eachLine.split(',')]
            statFuncName = items[1]
            if statFuncName == statFuncID:
                fileName = items[0]
                funcValue = float(items[2])
                runTime = float(items[3])
                pre = float(items[4])
                rec = float(items[5])
                fmeasure = float(items[6])
                fValues = [float(item) for item in items[7].split(' ')]
                while len(fValues) < iterNum:
                    fValues.append(fValues[len(fValues)-1])
                data.append([statFuncName, funcValue, runTime, pre, rec, fmeasure, fValues])
    averFValues = []
    for i in range(iterNum):
        total_FValue_i = []
        for item in data:
            total_FValue_i.append(item[6][i])
        averFValues.append(np.mean(total_FValue_i))
    print averFValues
    return [item for item in averFValues]
    
def drawAverReduced():
    funcid = 'EMS'
    averOurMethod_IHT = get_funcValues('graph_IHT_CitHepPh_Result.txt', funcid)
    averOurMethod_IHTP = get_funcValues('graph_GHTP_CitHepPh_Result_Updated.txt', funcid)
    averFusedLasso_0_25 = get_funcValues('genFusedLasso_CitHepPh_Result_0.25_.txt',funcid)
    averFusedLasso_1_0 = get_funcValues('genFusedLasso_CitHepPh_Result_1.0_.txt',funcid)
    averFusedLasso_5_0 = get_funcValues('genFusedLasso_CitHepPh_Result_5.0_.txt',funcid)
    averFusedLasso_10_0 = get_funcValues('genFusedLasso_CitHepPh_Result_10.0_.txt',funcid)
    size = len(averFusedLasso_0_25)
    x_ = range(size)
    markers = ['r>-','bs-','go-','cv-','yD-','m+-','k^-']
    #--------------------------------------------------------------------------------------------------------
    marker_fusedLasso_0_25, = plt.plot(x_, averFusedLasso_0_25, markers[1], markersize=10, label= 'Fused Lasso(0.25)',
             markerfacecolor = 'white',markeredgecolor='k',markeredgewidth=3.0)
    lin_fusedLasso_0_25, = plt.plot(x_, averFusedLasso_0_25, linestyle = '-',color = 'k',linewidth = 2.0, label= 'Fused Lasso(0.25)')
    #--------------------------------------------------------------------------------------------------------
    marker_fusedLasso_1_0, = plt.plot(x_, averFusedLasso_1_0, markers[2], markersize=10, label= 'Fused Lasso(1.0)',
             markerfacecolor = 'white',markeredgecolor='green',markeredgewidth=3.0)
    lin_fusedLasso_1_0, = plt.plot(x_, averFusedLasso_1_0, linestyle = '-',color = 'green',linewidth = 2.0, label= 'Fused Lasso(1.0)')
    #--------------------------------------------------------------------------------------------------------
    marker_fusedLasso_5_0, = plt.plot(x_, averFusedLasso_5_0, markers[3], markersize=10, label= 'Fused Lasso(5.0)',
             markerfacecolor = 'white',markeredgecolor='yellow',markeredgewidth=3.0)
    lin_fusedLasso_5_0, = plt.plot(x_, averFusedLasso_5_0, linestyle = '-',color = 'yellow',linewidth = 2.0, label= 'Fused Lasso(5.0)')
    #--------------------------------------------------------------------------------------------------------
    marker_fusedLasso_10_0, = plt.plot(x_, averFusedLasso_10_0, markers[5], markersize=10, label= 'Fused Lasso(10.0)',
             markerfacecolor = 'white',markeredgecolor='m',markeredgewidth=3.0)
    lin_fusedLasso_10_0, = plt.plot(x_, averFusedLasso_10_0, linestyle = '-',color = 'm',linewidth = 2.0, label= 'Fused Lasso(10.0)')
    #--------------------------------------------------------------------------------------------------------
    marker_ourMethod_IHT, = plt.plot(x_, averOurMethod_IHT, markers[0], markersize=10, label= 'IHT',
             markerfacecolor = 'white',markeredgecolor='b',markeredgewidth=3.0)
    lin_ourMethod_IHT, = plt.plot(x_, averOurMethod_IHT, linestyle = '-',color = 'b',linewidth = 2.0, label= 'IHT')
    #--------------------------------------------------------------------------------------------------------
    marker_ourMethod_IHTP, = plt.plot(x_, averOurMethod_IHTP, markers[4], markersize=10, label= 'IHTP',
             markerfacecolor = 'white',markeredgecolor='r',markeredgewidth=3.0)
    lin_ourMethod_IHTP, = plt.plot(x_, averOurMethod_IHTP, linestyle = '-',color = 'r',linewidth = 2.0, label= 'IHTP')
    #--------------------------------------------------------------------------------------------------------
    legend = plt.legend( [(marker_fusedLasso_0_25,lin_fusedLasso_0_25),
                          (marker_fusedLasso_1_0,lin_fusedLasso_1_0),
                          (marker_fusedLasso_5_0,lin_fusedLasso_5_0),
                          (marker_fusedLasso_10_0,lin_fusedLasso_10_0),
                          (marker_ourMethod_IHT,lin_ourMethod_IHT),
                         (marker_ourMethod_IHTP,lin_ourMethod_IHTP) ],
                        ['GenFusedLasso($\lambda=0.25$)',
                         'GenFusedLasso($\lambda=1.0$)',
                         'GenFusedLasso($\lambda=5.0$)',
                         'GenFusedLasso($\lambda=10.0$)',
                         'Graph-IHT','Graph-GHTP'], 
                        loc='lower right', shadow=True,prop={'size':14})
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(3)  # the legend line width
    plt.xlabel(r'iteration', fontsize=20)
    plt.ylabel(r''+funcid+' score', fontsize=20)
    plt.title('CitHepPh', fontsize=20)
    #plt.ylim([0.0,2000])
    #plt.savefig(os.path.join('./pics/','waterNetwork_lambdas'+funcid+'.png'))
    plt.show()
    #plt.close()

if __name__ == '__main__':
    drawAverReduced()