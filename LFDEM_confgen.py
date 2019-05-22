#!/usr/bin/env python3

import sys
import numpy as np
import argparse
import pyLFDEM as lfdem


def pvol(d, radius):
    if d == 2:
        return np.pi*radius**2
    elif d == 3:
        return (4/3)*np.pi*radius**3
    else:
        raise RuntimeError("d = "+str(d)+" unhandled...")


def get_volume(N, vf, vf_ratio, pvol1, pvol2):
    vf1 = vf*vf_ratio
    volume = N/(vf1/pvol1+vf1/pvol2)
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


def printOutConf(fname, positions, radii, params, fixed_velocities=None):
    with open(fname, "w") as outf:
        if params['walls']:
            header = "# np np1 np2 vf lx ly lz vf1 vf2 disp np_fixed format\n"+"# "
        else:
            header = "# np np1 np2 vf lx ly lz vf1 vf2 disp\n"+"# "

        for p in ['N', 'N1', 'N2', 'vf', 'lx', 'ly', 'lz', 'vf1', 'vf2']:
            header += " "+str(params[p])
        if params['walls']:
            header += " 0 "+ str(params['np_fixed'])+ " 6\n"
        else:
            header += " 0\n"

        outf.write(header)
        np_mobile = params['N']
        if "np_fixed" in params:
            np_mobile -= params['np_fixed']
        for i in range(np_mobile):
            outf.write(str(positions[i].x) + " " +
                       str(positions[i].y) + " " +
                       str(positions[i].z) + " " +
                       str(radii[i]) + "\n")
        for i in range(np_mobile, params['N']):
            outf.write(str(positions[i].x) + " " +
                       str(positions[i].y) + " " +
                       str(positions[i].z) + " " +
                       str(radii[i]) + " " +
                       str(fixed_velocities[i-np_mobile].x) + " " +
                       str(fixed_velocities[i-np_mobile].y) + " " +
                       str(fixed_velocities[i-np_mobile].z) + "\n")


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

    volume = get_volume(conf_params['N'],
                        conf_params['vf'],
                        conf_params['vf_ratio'],
                        pvol(conf_params['d'], conf_params['rad1']),
                        pvol(conf_params['d'], conf_params['rad2']))

    conf_params['N1'], conf_params['N2'] = get_N1N2(conf_params['N'],
                                                    conf_params['vf1'],
                                                    volume,
                                                    pvol(conf_params['d'], 1))
    conf_params['lx'], conf_params['ly'], conf_params['lz'] = \
        get_BoxSize(conf_params['d'], volume, ly_over_lx, lz_over_lx)

    conf_params['walls'] = args['walls']
    return conf_params


def setupSimu(simu, init_conf):
    rate = lfdem.DoubleDimQty()
    rate.unit = lfdem.Unit_hydro
    rate.dimension = lfdem.Dimension_Force
    rate.value = 0

    simu.setupFlow()

    PFactory = lfdem.ParameterSetFactory(lfdem.Unit_kn)
    sys = simu.getSys()

    sys.p = PFactory.getParameterSet()

    sys.p.flow_type = "shear"
    sys.p.kn = 1
    sys.p.friction_model = 0
    sys.p.integration_method = 0
    sys.p.disp_max = 5e-3
    sys.p.lubrication_model = "none"
    sys.p.contact_relaxation_time_tan = 1e-4

    sys.zero_shear=True
    sys.mobile_fixed = True

    simu.assertParameterCompatibility()
    sys.setupConfiguration(init_conf, lfdem.ControlVariable_rate)
    sys.setVelocityDifference()
    print("setup [ok]")



def generateConf(simu, conf_params,
                 stop_params={'contact_ratio': 0.05,
                              'min_gap': -0.01}):
    """
        Generate an initial set of positions/radii that is relaxed up
        to small remaining overlaps.

        conf_params in a dict of parameters that must contain:
        'd', 'N', 'N1', 'rad1', 'rad2', 'lx', 'ly', 'lz'.

    """
    if not conf_params['walls']:
        conf = lfdem.base_shear_configuration()
        conf.position, conf.radius = initialRandom(conf_params['N'],
                                                 conf_params['N1'],
                                                 conf_params['rad1'],
                                                 conf_params['rad2'],
                                                 conf_params['lx'],
                                                 conf_params['ly'],
                                                 conf_params['lz'])

        conf.volume_or_area_fraction = conf_params['vf']
        conf.lx, conf.ly, conf.lz = conf_params['lx'], conf_params['ly'], conf_params['lz']
    else:
        conf_position = np.random.rand(conf_params['N'], 3)
        conf_position[:, 0] *= conf_params['lx']
        conf_position[:, 1] *= conf_params['ly']
        conf_position[:, 2] *= conf_params['lz']

        conf_radius = conf_params['rad1']*np.ones(conf_params['N'])
        conf_radius[conf_params['N1']:] *= conf_params['rad2']

        conf = lfdem.fixed_velo_configuration()

        conf.volume_or_area_fraction = conf_params['vf']
        conf.lx, conf.ly = conf_params['lx'], conf_params['ly']

        wall_part_rad = 1

        wall_part_nb_x = int(conf.lx/wall_part_rad) + 1
        wall_part_nb_y = int(conf.ly/wall_part_rad) + 1
        x_spacing = conf.lx/wall_part_nb_x
        y_spacing = conf.ly/wall_part_nb_y

        wall_part_nb = wall_part_nb_x*wall_part_nb_y
        wall_part_pos = np.zeros((wall_part_nb, 3))
        wall_part_pos[:, 0] = np.tile(np.linspace(0,
                                                  conf.lx,
                                                  wall_part_nb_x,
                                                  endpoint=False),
                                      wall_part_nb_y)\
                                      + 0.5*x_spacing
        wall_part_pos[:, 1] = np.repeat(np.linspace(0,
                                                    conf.ly,
                                                    wall_part_nb_y,
                                                    endpoint=False),
                                        wall_part_nb_x)\
                                        + 0.5*y_spacing

        #up wall
        wall_part_pos[:, 2] = conf_params['lz'] + wall_part_rad
        conf_position = np.row_stack((conf_position, wall_part_pos))
        conf_radius = np.append(conf_radius, wall_part_rad*np.ones(wall_part_nb))
        #down wall
        wall_part_pos[:, 2] = -wall_part_rad
        conf_position = np.row_stack((conf_position, wall_part_pos))
        conf_radius = np.append(conf_radius, wall_part_rad*np.ones(wall_part_nb))

        conf_fixed_velocities = np.zeros((2*wall_part_nb, 3))
        conf_fixed_velocities[:wall_part_nb, 0] = 1
        conf_fixed_velocities[wall_part_nb:, 0] = -1

        # trick LF_DEM PBCs along z by extending the system along z
        conf_params['lz'] += 4*wall_part_rad
        conf_position[:, 2] += 2*wall_part_rad

        conf.radius = lfdem.DoubleVector(conf_radius)
        conf.position = lfdem.Vec3dVector()
        for p in conf_position:
            pos = lfdem.vec3d()
            pos.x = p[0]
            pos.y = p[1]
            pos.z = p[2]
            conf.position.push_back(pos)

        conf.fixed_velocities = lfdem.Vec3dVector()
        for v in conf_fixed_velocities:
            vel = lfdem.vec3d()
            vel.x = v[0]
            vel.y = v[1]
            vel.z = v[2]
            conf.fixed_velocities.push_back(vel)


        conf.lz = conf_params['lz']
        conf_params['np_fixed'] = conf.fixed_velocities.size()
        conf_params['N'] += conf_params['np_fixed']

    setupSimu(simu, conf)

    sys = simu.getSys()
    if conf_params['walls']:
        simu.p.simulation_mode = 31
    sys.checkNewInteraction()
    sys.updateInteractions()
    contact_nb = lfdem.countNumberOfContact(sys)
    print("Generating", flush=True, end='')
    while contact_nb[0]/conf_params['N'] > stop_params['contact_ratio']\
            and lfdem.evaluateMinGap(sys) < stop_params['min_gap']:
        sys.timeEvolution(sys.get_time()+2, -1)
        contact_nb = lfdem.countNumberOfContact(sys)
        # print(".", flush=True, end='')
        print(contact_nb[0]/conf_params['N'], lfdem.evaluateMinGap(sys))

    if not conf_params['walls']:
        return np.array(sys.position), np.array(sys.radius)
    else:
        return np.array(sys.position), np.array(sys.radius), np.array(sys.fixed_velocities)



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
    parser.add_argument('-w', '--walls', action='store_true')
    parser.add_argument('-s', '--soft', action='store_true')


    args = vars(parser.parse_args(sys.argv[1:]))

    conf_params = expandConfParams(**args)

    conf_name = genConfName(conf_params)

    simu = lfdem.Simulation()
    print(" LF_DEM version : ", simu.gitVersion())

    if args['soft']:
        conf = generateConf(simu, conf_params,
                            {'contact_ratio': 8,
                            'min_gap': -2})
    else:
        conf = generateConf(simu, conf_params)


    print("\n==========")
    print("Conf : ", conf_name)
    if conf_params['walls']:
        printOutConf(conf_name, conf[0], conf[1], conf_params, fixed_velocities=conf[2])
    else:
        printOutConf(conf_name, conf[0], conf[1], conf_params)
