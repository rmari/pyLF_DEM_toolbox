#!/usr/bin/env python

import sys
import numpy as np
import Simulation_LFDEM as lfdem


def get_pvols(d, radius1, radius2):
    if d==2:
        pvol1 = np.pi*radius1**2
        pvol2 = np.pi*radius2**2
    else:
        if d==3:
            pvol1 = (4/3)*np.pi*radius1**3
            pvol2 = (4/3)*np.pi*radius2**3
        else:
            print("d = ", d, "???")

    return pvol1, pvol2

def get_volume(N,vf,vf_ratio,pvol1,pvol2):
    vf1 = vf*vf_ratio
    vf2 = vf-vf1

    volume = N/(vf1/pvol1+vf1/pvol2)

    return volume

def get_N1N2(N, vf1, volume, pvol1):
    N1 = int(np.around(vf1*volume/pvol1))
    N2 = N - N1
    return N1, N2

def get_BoxSize(d,volume, lylx, lzlx):
    if d==2:
        lx = (volume/lzlx)**(1/2)
        ly = 0
        lz  = lzlx*lx
    if d==3:
        lx = (volume/lylx/lzlx)**(1/3)
        ly = lylx*lx
        lz = lzlx*lx
    return lx,ly,lz

def initialRandom(N, N1, radius1, radius2, lx, ly, lz):
    positions = lfdem.Vec3dVector()

    for i in range(N):
        pos = lfdem.vec3d()
        pos.x = lx*np.random.rand()
        pos.y = ly*np.random.rand()
        pos.z = lz*np.random.rand()
        positions.push_back(pos)

    radii = radius1*np.ones(N)
    radii[N1:] = radius2

    radii = lfdem.DoubleVector(radii)

    return positions, radii

def printOutConf(fname,positions, radii, params):
    with open(fname,"w") as outf:
        header = "# np1 np2 vf lx ly lz vf1 vf2 disp\n"+"# "
        for p in ['N1', 'N2', 'vf', 'lx', 'ly', 'lz', 'vf1', 'vf2']:
            header +=" "+str(params[p])
        header += " 0\n"

        outf.write(header)
        for i in range(len(positions)):
            outf.write(str(positions[i].x)+" "+str(positions[i].y)+" "+str(positions[i].z)+" "+str(radii[i])+"\n")


if len(sys.argv) != 3 or len(sys.argv) != 5:
    print(sys.argv[0]+" N vf\n")
    exit(1)

radius1 = 1
radius2 = 1.4

vf_ratio = 0.5

ly_over_lx = 4
lz_over_lx = 1

conf_params = {}

conf_params['d']=3
conf_params['N']= int(sys.argv[1])
conf_params['vf'] = float(sys.argv[2])
if len(sys.argv) == 5:
    ly_over_lx = float(sys.argv[3])
    lz_over_lx = float(sys.argv[4])
else:
    ly_over_lx = 1
    lz_over_lx = 1

conf_params['vf1'] = conf_params['vf']*vf_ratio
conf_params['vf2'] = conf_params['vf'] - conf_params['vf1']



simu = lfdem.Simulation()
sys = simu.getSys()

pvol1, pvol2 =  get_pvols(conf_params['d'], radius1, radius2)
volume = get_volume(conf_params['N'],conf_params['vf'],vf_ratio,pvol1,pvol2)
conf_params['N1'], conf_params['N2'] =  get_N1N2(conf_params['N'],conf_params['vf1'], volume, pvol1)
conf_params['lx'], conf_params['ly'], conf_params['lz'] = get_BoxSize(conf_params['d'],volume, ly_over_lx, lz_over_lx)

positions, radii = initialRandom(conf_params['N'], conf_params['N1'], radius1, radius2,conf_params['lx'], conf_params['ly'], conf_params['lz'])

simu.setDefaultParameters()
simu.p.np_fixed = 0
sys.zero_shear = True
simu.p.kn = 1
simu.p.friction_model = 0
simu.p.integration_method = 0
simu.p.disp_max = 5e-3
simu.p.lubrication_model = 0
# simu.p.contact_relaxation_time = 1e-4

sys.setConfiguration(positions, radii,conf_params['lx'], conf_params['ly'], conf_params['lz'])
sys.setupSystem("rate")

sys.checkNewInteraction()
sys.updateInteractions()

# print(conf_params, volume, pvol1, pvol2)
sys.analyzeState()
print(sys.contact_nb,sys.min_reduced_gap,flush=True)

while sys.contact_nb/conf_params['N']>0.05 and sys.min_reduced_gap<-0.01:
    sys.timeEvolution("time", sys.get_time()+2)
    sys.analyzeState()
    print(sys.contact_nb,sys.min_reduced_gap,flush=True)

fname = "D"+str(conf_params['d'])+"N"+str(conf_params['N'])+"VF"+str(conf_params['vf'])
fname += "Bidi"+str(radius2/radius1)
fname += str(vf_ratio)
if conf_params['d']==3:
    fname += "LyLx"+str(ly_over_lx)
fname += "LzLx"+str(lz_over_lx)
fname += ".dat"

printOutConf(fname,sys.position, radii, conf_params)
