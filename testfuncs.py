import numpy as np
import scipy as sp
import os,sys
import datetime
import math
import json
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

def dif_eqs2(t, f, v, Total_Cross_Section, He, Cross_Sections_dict, A):

    df0dt = -Total_Cross_Section * v * He * f[0] 
    df1dt = Cross_Sections_dict['1s']*v*He*f[0] + A['2s']['1s']*f[2] + A['2p']['1s']*f[3] + A['3p']['1s']*f[5] + A['3s']['1s']*f[4] + A['3d']['1s']*f[6] + A['4s']['1s']*f[7]
    df2dt = Cross_Sections_dict['2s']*v*He*f[0] + A['3p']['2s']*f[5] + A['4p']['2s']*f[8] - A['2s']['1s']*f[2]
    df3dt = Cross_Sections_dict['2p']*v*He*f[0] + A['4s']['2p']*f[7] + A['4d']['2p']*f[9] + A['4f']['2p']*f[10] + A['3s']['2p']*f[4] + A['3d']['2p']*f[6] - A['2p']['1s']*f[3]
    df4dt = Cross_Sections_dict['3s']*v*He*f[0] + A['4p']['3s']*f[8] - (A['3s']['2p']+A['3s']['1s'])*f[4]
    df5dt = Cross_Sections_dict['3p']*v*He*f[0] + A['4d']['3p']*f[9] + A['4s']['3p']*f[7] - (A['3p']['2s'] + A['3p']['1s'])*f[5]
    df6dt = Cross_Sections_dict['3d']*v*He*f[0] + A['4p']['3d']*f[8] - (A['3d']['2p']+A['3d']['1s'])*f[6]    
    df7dt = Cross_Sections_dict['4s']*v*He*f[0] - (A['4s']['3p']+A['4s']['2p']+A['4s']['1s'])*f[7]
    df8dt = Cross_Sections_dict['4p']*v*He*f[0] - (A['4p']['3s']+A['4p']['2s']+A['4p']['3d'])*f[8]
    df9dt = Cross_Sections_dict['4d']*v*He*f[0] - (A['4d']['3p']+A['4d']['2p'])*f[9]
    df10dt= Cross_Sections_dict['4f']*v*He*f[0] - (A['4f']['2p'])*f[10]
    dfdt = [df0dt, df1dt, df2dt, df3dt, df4dt, df5dt, df6dt, df7dt, df8dt, df9dt, df10dt]

    return dfdt

def selectFile(directory):
    try:
        Files = os.listdir(directory)
        for i, File in enumerate(Files):
            print('[{0}] {1}'.format(i, File))
        Filenumber = int(input())
        print(Files[Filenumber])
        Filepath = directory + '/' + Files[Filenumber]
    except:
        # 入力したファイル名が間違っているか，ファイルが存在しない場合，プログラムを終了する．
        print('入力したファイルは存在しません．')
        sys.exit(1)

    return Filepath

def make2Ddict(names):
    dict_2d = {}

    for x in names:
        for y in names:
            if x in dict_2d:
                dict_2d[x][y] = 0.0
            else:
                dict_2d[x] = {y:0.0}

    return dict_2d

def convert_energy(col_Energy_kVu):

    col_Energy_cms = math.sqrt(col_Energy_kVu*10**10*9.64854)*10**2
    
    return col_Energy_cms

def plot_populations(xaxiss, populations, orbits):
    
    if len(xaxiss) != 1:
        selectXaxis = input('横軸を選んでください(time(t) or position(x))')
        if selectXaxis == 'time' or selectXaxis == 't':
            xaxisarray = xaxiss[0]
        else:
            xaxisarray = xaxiss[1]
    else:
        xaxisarray = xaxiss[0]

    isPlotBeforeCollision = input('衝突前のイオンのポピュレーションをプロットしますか?(y/n) > ')

    fig = plt.figure(figsize=(16,9))
    popfig = fig.add_subplot(1,1,1)
    if isPlotBeforeCollision == 'Y' or isPlotBeforeCollision == 'y':
        popfig.plot(xaxisarray, populations[:,0], label='O7+')
    for i in np.arange(1,11):
        popfig.plot(xaxisarray[::100], populations[::100,i], label='{0}'.format(orbits[i-1]))

    if selectXaxis == 'time' or selectXaxis == 't':
        popfig.set_xlabel('Time[s]', fontsize=25)
    else:
        popfig.set_xlabel('Position[cm]', fontsize=25)
    
    popfig.set_ylabel('Population', fontsize=25)
    for tick in popfig.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in popfig.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    popfig.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    popfig.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    popfig.xaxis.offsetText.set_fontsize(15)
    popfig.yaxis.offsetText.set_fontsize(15)
    popfig.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    popfig.xaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))

    #popfig.legend(loc='best')
    plt.show()

    return fig

def outputjson(parameter_dic, output_filename='lastparameter.json'):
    '''
        プログラムの実行に使ったパラメータ群をjsonファイルに出力
    '''
    output_file = open(output_filename, 'w')
    json.dump(parameter_dic, output_file)

def inputjson(input_filename='lastparameters.json'):
    '''
        前回のプログラム実行に使ったパラメータ群を辞書に格納
    '''
    f = open(input_filename)
    param_dic = json.load(f)

    return param_dic
