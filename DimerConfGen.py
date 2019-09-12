import sys
import numpy as np
import pyLFDEM as lfdem
import LFDEM_confgen as lfcg
import argparse

def dimervol(d, radius, sep):
    hcap = radius*(1 - sep)
    if d == 2:
        theta = 2*np.arccos(1-hcap/radius)
        areaseg = radius**2*(theta-np.sin(theta))/2
        areap = np.pi*radius**2
        return 2*(areap - areaseg)
    elif d == 3:
        volcap = np.pi*(hcap**2)*(3*radius - hcap)/3
        volp = (4./3.)*np.pi*radius**3
        return 2*(volp - volcap)
    else:
        raise RuntimeError("d = "+str(d)+" unhandled...")


def get_volume(N, vf, vf_ratio, pvol1, pvol2):
    vf1 = vf*vf_ratio
    volume = (N/vf)/((1-vf_ratio)/pvol2 + vf_ratio/pvol1)
    print("v ", volume, ((1-vf_ratio)/pvol2 + vf_ratio/pvol1))
    return volume


def get_N1N2(N, vf1, volume, pvol1):
    N1 = int(np.around(vf1*volume/pvol1))
    N2 = N - N1
    return N1, N2


def get_BoxSize(d, volume, lylx, lzlx):
    if d == 2:
        lx = (volume/lzlx)**(1/2)
        ly = 0
        lz = lzlx*lx
    if d == 3:
        lx = (volume/lylx/lzlx)**(1/3)
        ly = lylx*lx
        lz = lzlx*lx
    return lx, ly, lz

def randUnit3d():
    n = np.random.randn(3)
    return n/np.linalg.norm(n)

def randUnit2d():
    n = np.random.randn(3)
    n[1] = 0
    return n/np.linalg.norm(n)

def periodize(p, box):
    for i in range(3):
        while p[i] > box[i]:
            p[i] -= box[i]
        while p[i] > box[i]:
            p[i] -= box[i]
    return p

def initialRandom(p):
    positions = lfdem.Vec3dVector()
    box = np.array([p['lx'], p['ly'], p['lz']])
    for i in range(p['N']):
        
        pos = box*np.random.rand(3)
        positions.push_back(lfdem.vec3d(*pos))
        if p['d']==2:
            n = p['sep']*randUnit2d()
        if p['d']==3:
            n = p['sep']*randUnit3d()
        positions.push_back(lfdem.vec3d(*(periodize(pos+n, box))))
    radii = p['rad1']*np.ones(2*p['N'])
    radii[2*p['N1']:] = p['rad2']
    radii = lfdem.DoubleVector(radii)

    return positions, radii


def expandConfParams(**args):
    """
        Get a full set of parameters from independent parameters.

        Possible args:
        * 'N' (required)
        * 'volume-fraction' (required)
        * 'dimension' [3]
        * 'radius_ratio' [1]: bidispersity size ratio: one species will have
                              radius 1, the other radius 'radius-ratio'
        * 'volume_fraction_ratio' [0.5]: vf ratio of particles with radius!=1
                                         over particles with radius=1
        * 'gradient_flow_ratio' [1]: size ratio gradient direction
                                     over flow direction
        * 'vorticity_flow_ratio' [1]: size ratio vorticity direction
                                     over flow direction

    """
    conf_params = {}
    conf_params['N'] = args['N']
    conf_params['vf'] = args['volume_fraction']
    conf_params['d'] = args.pop('dimension', 3)
    conf_params['vf_ratio'] = args.pop('volume_fraction_ratio', 0.5)
    conf_params['rad1'] = 1
    conf_params['rad2'] = args.pop('radius_ratio', 1)
    ly_over_lx = args.pop('vorticity_flow_ratio', 1)
    lz_over_lx = args.pop('gradient_flow_ratio', 1)
    conf_params['vf1'] = conf_params['vf']*conf_params['vf_ratio']
    conf_params['vf2'] = conf_params['vf'] - conf_params['vf1']
    conf_params['sep'] = args['separation']

    volume = get_volume(conf_params['N'],
                        conf_params['vf'],
                        conf_params['vf_ratio'],
                        dimervol(conf_params['d'], conf_params['rad1'], conf_params['sep']),
                        dimervol(conf_params['d'], conf_params['rad2'], conf_params['sep']))

    conf_params['N1'], conf_params['N2'] = get_N1N2(conf_params['N'],
                                                    conf_params['vf1'],
                                                    volume,
                                                    dimervol(conf_params['d'], conf_params['rad1'], conf_params['sep']))
    conf_params['lx'], conf_params['ly'], conf_params['lz'] = \
        get_BoxSize(conf_params['d'], volume, ly_over_lx, lz_over_lx)
    conf_params['walls'] = False
    return conf_params

def dimerConf(conf_params):
    conf = lfdem.base_shear_configuration()
    conf.position, conf.radius = initialRandom(conf_params)

    conf.volume_or_area_fraction = conf_params['vf']
    conf.lx, conf.ly, conf.lz = conf_params['lx'], conf_params['ly'], conf_params['lz']

    dimers = lfdem.UnloadedDimerStateVector()
    for i in range(conf_params['N']):
        ds = lfdem.UnloadedDimerState()
        ds.p0 = 2*i
        ds.p1 = 2*i+1
        ds.relaxed_length = 2*conf.radius[2*i]*conf_params['sep']
        dimers.push_back(ds)
    return conf, dimers


def setupSimu(simu, init_conf, dimers):
    # rate = lfdem.DoubleDimQty()
    # rate.unit = lfdem.Unit_hydro
    # rate.dimension = lfdem.Dimension_Force
    # rate.value = 0

    simu.setupFlow()

    PFactory = lfdem.ParameterSetFactory(lfdem.Unit_kn)
    sys = simu.getSys()

    sys.p = PFactory.getParameterSet()

    sys.p.integration_method = 0
    sys.p.disp_max = 5e-3
    sys.p.contact.relaxation_time_tan = 1e-4
    sys.p.contact.friction_model = 0
    sys.p.contact.kn = 1

    sys.mobile_fixed = False
    sys.shear_type = False

    simu.assertParameterCompatibility()
    sys.setupConfiguration(init_conf, lfdem.ControlVariable_rate)
    sys.addDimers(dimers)

    print("setup [ok]")



def runRelax(simu, conf, dimers, conf_params,
                 stop_params={'contact_ratio': 0.05,
                              'min_gap': -0.01}):
    """
        Generate an initial set of positions/radii that is relaxed up
        to small remaining overlaps.

        conf_params in a dict of parameters that must contain:
        'd', 'N', 'N1', 'rad1', 'rad2', 'lx', 'ly', 'lz'.

    """
    for i in range(2*conf_params['N']):
        conf.radius[i] *= 1.02

    setupSimu(simu, conf, dimers)

    simusys = simu.getSys()
    
    contact_nb = lfdem.countNumberOfContact(simusys)
    print("Generating", flush=True, end='')
    while contact_nb[0]/conf_params['N'] > stop_params['contact_ratio']:
        simusys.timeEvolution(simusys.get_time()+2, -1)
        contact_nb = lfdem.countNumberOfContact(simusys)
        # print(".", flush=True, end='')
        print(contact_nb[0]/conf_params['N'], lfdem.evaluateMinGap(simusys))

    for i in range(2*conf_params['N']):
        simusys.conf.radius[i] /= 1.02
    return simu

def printOutConf(fname, positions, radii, params):
    with open(fname, "w") as outf:
        header = "# np np1 np2 vf lx ly lz vf1 vf2 disp\n"+"# "

        for p in ['N', 'N1', 'N2', 'vf', 'lx', 'ly', 'lz', 'vf1', 'vf2']:
            header += " "+str(params[p])
        header += " 0\n"

        outf.write(header)
        np_mobile = 2*params['N']
        if "np_fixed" in params:
            np_mobile -= params['np_fixed']
        for i in range(np_mobile):
            outf.write(str(positions[i].x) + " " +
                       str(positions[i].y) + " " +
                       str(positions[i].z) + " " +
                       str(radii[i]) + "\n")

def printOutDimers(fname, dimers, params):
    with open(fname, "w") as outf:
        for i in range(params['N']):
            outf.write(str(dimers[i].p0) + " " +
                       str(dimers[i].p1) + " " +
                       str(dimers[i].relaxed_length) + "\n")


def generateConf(simu, conf_params,
                 stop_params={'contact_ratio': 0.05,
                              'min_gap': -0.01}):
    """
        Generate an initial set of positions/radii that is relaxed up
        to small remaining overlaps.

        conf_params in a dict of parameters that must contain:
        'd', 'N', 'N1', 'rad1', 'rad2', 'lx', 'ly', 'lz'.

    """
    conf, dimers = dimerConf(conf_params)

    simu = runRelax(simu, conf, dimers, conf_params, stop_params)

    simusys = simu.getSys()
    return np.array(simusys.conf.position), np.array(simusys.conf.radius), dimers


def argParse(in_args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-N',
                        type=int, required=True)
    parser.add_argument('-vf', '--volume-fraction',
                        type=float, required=True)
    parser.add_argument('-vfr', '--volume-fraction-ratio',
                        type=float, default=0.5)
    parser.add_argument('-zx', '--gradient-flow-ratio',
                        type=float, default=1)
    parser.add_argument('-yx', '--vorticity-flow-ratio',
                        type=float, default=1)
    parser.add_argument('-d', '--dimension',
                        type=int, default=3)
    parser.add_argument('-rr', '--radius-ratio',
                        type=float, default=1.4)
    parser.add_argument('-o', '--output')
    parser.add_argument('-s', '--separation',
                        type=float, required=True)

    return vars(parser.parse_args(in_args))

def genConfName(conf_params):
    fname = "D"+str(conf_params['d'])
    fname += "N"+str(conf_params['N'])
    fname += "VF"+str(conf_params['vf'])

    if conf_params['rad2']/conf_params['rad1'] != 1:
        fname += "Bidi"+str(conf_params['rad2']/conf_params['rad1'])
        fname += str(conf_params['vf_ratio'])
    else:
        fname += "Mono"

    if conf_params['d'] == 3:
        fname += "LyLx"+str(conf_params['ly']/conf_params['lx'])
    fname += "LzLx"+str(conf_params['lz']/conf_params['lx'])
    if conf_params['walls']:
        fname += "_walls"
    fname += ".dat"

    return fname

if __name__ == '__main__':
    args = argParse(sys.argv[1:])
    conf_params = expandConfParams(**args)

    conf_name = genConfName(conf_params)
    dimer_name = "Dimer"+conf_name
    simu = lfdem.Simulation()
    print(" LF_DEM version : ", simu.gitVersion())

    conf, rad, dimers = generateConf(simu, conf_params,
                        {'contact_ratio': 8,
                        'min_gap': -2})


    print("\n==========")
    print("Conf : ", conf_name)
    printOutConf(conf_name, conf, rad, conf_params)
    printOutDimers(dimer_name, dimers, conf_params)