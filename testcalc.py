import numpy as np
import os,sys
import math
import datetime
from scipy.integrate import odeint
from scipy.linalg import block_diag
from testfuncs import *
import matplotlib.pyplot as plt

def testcalc():
    
    orbits = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f']
    d = datetime.datetime.today()
    default_output_filename = d.strftime('%y%m%d_%H%M_')

    print('電子捕獲断面積ファイルを選んでください．')
    ECCSFilePath = selectFile('ECCS')
    ECCSs = np.loadtxt(ECCSFilePath, delimiter=',')

    print('どの衝突エネルギーで計算を行いますか．')
    for (i, Show_Col_E) in enumerate(ECCSs[:,0]):
        print(str(i) + ' ' + str(Show_Col_E))

    #衝突エネルギー[keV/u]を決める
    Collision_Energy = int(input())
    print('{0} keV/u'.format(ECCSs[Collision_Energy,0]))
    default_output_filename = '{0}{1}keVu'.format(default_output_filename,ECCSs[Collision_Energy,0])

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
    
    print('A係数ファイルを選んでください．')
    ACFilePath = selectFile('AC')
    ACs = np.loadtxt(ACFilePath, delimiter=',')

    # A係数を格納する2次元辞書
    AC_dict = make2Ddict(orbits)

    for (num1, i) in enumerate(orbits):
        for (num2, f) in enumerate(orbits):
            AC_dict[i][f] = ACs[num1][num2]

    #print(AC_dict)

    #初期条件設定(衝突前なのでprimary ionが100%
    He = float(input('Heの粒子数を入力してください(単位cm^-3) > '))
    f0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    f0_3 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    default_output_filename = '{0}_{1}percmcube'.format(default_output_filename, He)
    print(default_output_filename)
 
    #計算を行うtime regionを決める
    t = np.linspace(0,2.5e-7,100000)

    #時間とイオンビームの速度から位置を決める
    x = t * Collision_Speed

    #微分方程式をodeintに解かせる
    sol, infodict = odeint(dif_eqs2, f0, t, args=(Collision_Speed,Total_Cross_Section,He,Cross_Sections_dict,AC_dict),full_output=True, printmessg=True)    
    #print('len(sol) = {0}'.format(len(sol)))

    #for key, val in infodict.items():
        #print('infodict[\'{0}\'] = {1}'.format(key, val))


    fig = plot_populations([t,x], sol, orbits)

    saveyn = input('グラフを保存しますか?([y]/n) > ')
    if saveyn == 'y' or saveyn == 'Y':
        figoutput(default_output_filename, fig)

if __name__ == '__main__':

    print('This is the main function.')
    
    testcalc()
