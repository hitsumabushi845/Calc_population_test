import numpy as np
import os,sys
import math
import gc
import json
import datetime
from scipy.integrate import odeint
from scipy.linalg import block_diag
from testfuncs import *
import matplotlib.pyplot as plt

def testcalc():
    
    orbits = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f']

    if os.path.exists('lastparameters.json'):
        do_you_read_jsonfile = input('前回の設定を使用しますか?[y]/n > ')
        if do_you_read_jsonfile == 'y' or do_you_read_jsonfile == 'Y':
            param_dic = inputjson()
        else:
            param_dic = {}
    else:
        param_dic = {}

    if 'ECCS_Filename' in param_dic:
        ECCSs = np.loadtxt(param_dic['ECCS_Filename'], delimiter=',')
    else:
        print('計算したい電子捕獲断面積ファイルの番号を選んでください．')
        ECCSFilePath = selectFile('ECCS')
        ECCSs = np.loadtxt(ECCSFilePath, delimiter=',')
        param_dic['ECCS_Filename'] = ECCSFilePath 

    if 'Col_Energy' in param_dic:
        Collision_Speed = convert_energy(ECCSs[param_dic['Col_Energy'],0])
        Total_Cross_Section = ECCSs[param_dic['Col_Energy'],1] * 1.e-16
        Cross_Sections_list = ECCSs[param_dic['Col_Energy'],2:]
    else:
        print('計算したい衝突エネルギーの番号を選んでください．')
        for (i, Show_Col_E) in enumerate(ECCSs[:,0]):
            print('[{0}] {1}'.format(i,Show_Col_E))

        #衝突エネルギー[keV/u]を決める
        Collision_Energy = int(input())
        print('{0} keV/u'.format(ECCSs[Collision_Energy,0]))
        default_output_filename = '{0}{1}keVu'.format(default_output_filename,ECCSs[Collision_Energy,0])
        param_dic['Col_Energy'] = Collision_Energy

        #衝突エネルギーから衝突速度計算
        Collision_Speed = convert_energy(ECCSs[Collision_Energy,0])
        print('{0} cm/s'.format(Collision_Speed))

        #衝突エネルギーに対応する電子捕獲断面積をCross_Sections_listに格納，Total_Cross_Sectionに全断面積を格納
        Total_Cross_Section = ECCSs[Collision_Energy,1] * 1.e-16
        Cross_Sections_list = ECCSs[Collision_Energy,2:]

    #電子軌道の名前をkeyに持つ辞書に電子捕獲断面積を格納
    #電子捕獲断面積の単位が10^-16 cm^2なので，cm^2に修正
    Cross_Sections_dict = {}
    for i, orbit in enumerate(orbits):
        Cross_Sections_dict[orbit] = Cross_Sections_list[i] * 1.e-16
    #print(Cross_Sections_dict)
    del Cross_Sections_list
    gc.collect()
    
    if 'AC_Filename' in param_dic:
        ACs = np.loadtxt(param_dic['AC_Filename'], delimiter=',')
    else:
        print('計算したいA係数ファイルの番号を選んでください．')
        ACFilePath = selectFile('AC')
        ACs = np.loadtxt(ACFilePath, delimiter=',')
        param_dic['AC_Filename'] = ACFilePath

    # A係数を格納する2次元辞書
    AC_dict = make2Ddict(orbits)

    for (num1, i) in enumerate(orbits):
        for (num2, f) in enumerate(orbits):
            AC_dict[i][f] = ACs[num1][num2]
    del ACs
    gc.collect()

    if 'particle_number' in param_dic:
        He = param_dic['particle_number']
    else:
        #初期条件設定(衝突前なのでprimary ionが100%
        He = float(input('Heの粒子数を入力してください(単位cm^-3) > '))
        param_dic['particle_number'] = He
    f0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 
    #計算を行う時間幅，ステップ幅，t0を決める
    tend = 2.5e-7
    dt = 1.e-14
    t0 = 0.0

    r = ode(dif_eqs_forode).set_integrator('vode',method='bdf')
    r.set_initial_value(f0,t0).set_f_params(Collision_Speed,Total_Cross_Section,He,Cross_Sections_dict,AC_dict)

    t1 = []
    res = []
    print('Enter 1.e-14 calculation')
    counter = 1
    while r.successful() and r.t < tend:
        if counter % 10000 == 0:
            print('counter = {0}'.format(counter))
        tmp = r.integrate(r.t+dt)
        t1.append(r.t)
        res.append(tmp)
        counter = counter + 1

    t1 = np.array(t1)
    sol1 = np.array(res)

    x1 = t1 * Collision_Speed

    if 'Horizontal_axis' in param_dic:
        if param_dic['Horizontal_axis'] == 't':
            fig = plot_populations([t], sol, orbits)
        else:
            fig = plot_populations([x], sol, orbits)
    else:
        fig, param_dic['Horizontal_axis'] = plot_populations([t,x], sol, orbits)

    plt.clf()

    if do_you_read_jsonfile != 'y' or do_you_read_jsonfile != 'Y':
        f = open('lastparameters.json', 'w')
        json.dump(param_dic, f)

if __name__ == '__main__':

    print('This is the main function.')
    
    testcalc()
