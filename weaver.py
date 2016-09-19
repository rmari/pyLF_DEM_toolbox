#!/usr/bin/env python3

import argparse
import numpy as np
from scipy import special
import re
import pyLFDEM as lfdem
from shutil import copyfile
import sys


def getArgParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    parser.add_argument('params_file', type=str)
    parser.add_argument('-n', '--binary_conf', action='store_true')
    parser.add_argument('-f', '--overwrite', action='store_true')
    parser.add_argument('-r', '--rate_primary', type=str, required=False)
    parser.add_argument('-s', '--shear_stress', type=str, required=False)
    parser.add_argument('-rosp', '--rate_OSP_max_ratio',
                        type=float, required=True)
    parser.add_argument('-aosp', '--amplitude_OSP', type=float, required=True)

    return parser


def echoInput(in_args, simu):
    echofile = open("input_data_"+simu.simu_name+".dat", "w")
    configfilename = "input_config_"+simu.simu_name+".dat"
    echofile.write("pyLF_DEM version "+simu.gitVersion()+"\n" +
                   "called as:\n" +
                   ' '.join(sys.argv)+"\n\n")

    echofile.write(in_args['config_file'] + " is copied in " +
                   configfilename + "\n\n")

    echofile.write("Script file "+__file__+" is:\n\n")
    for line in open(__file__, "r"):
        echofile.write(line)
    echofile.write("\n")

    echofile.write("Parameter file "+in_args['params_file']+" is:\n\n")

    for line in open(in_args['params_file'], "r"):
        echofile.write(line)
    echofile.close()
    print("echoed")
    if in_args['binary_conf']:
        configfilename += ".bin"
    else:
        configfilename += ".dat"
    copyfile(in_args['config_file'], configfilename)


def getSimuName(in_args, control_var):
    conf_name = in_args['config_file'].replace(".dat", "").replace(".bin", "")
    conf_name = conf_name[conf_name.rfind("/")+1:]
    param_name = in_args['params_file'].replace(".txt", "").replace(".dat", "")
    param_name = param_name[param_name.rfind("/")+1:]

    simu_name = conf_name + "_" + param_name
    if control_var == "rate":
        simu_name += "_rate"+str(in_args['rate_primary'])
    else:
        simu_name += "_stress"+str(in_args['shear_stress'])
    simu_name += "_rate"+str(in_args['rate_primary'])\
                 + "_ratioOSP"+str(in_args['rate_OSP_max_ratio'])\
                 + "_amplitudeOSP"+str(in_args['amplitude_OSP'])
    return simu_name


def getSuffix(value_str):
    split = re.search('[a-zA-Z]', value_str).start()
    return float(value_str[:split]), value_str[split:]


def setupSimulation(in_args, simu):
    """ Setup for weaving simulations """

    print("Weaving simulation.")
    print("Based on pyLFDEM version", simu.gitVersion())

    simu.force_to_run = in_args['overwrite']

    unit = str()
    primary_rate = str()
    shear_stress = str()
    try:
        primary_rate, unit = getSuffix(in_args['rate_primary'])
        control_var = "rate"
    except TypeError:  # stress control
        try:
            shear_stress, unit = getSuffix(in_args['shear_stress'])
            control_var = "stress"
        except TypeError:
            print("No rate nor stress specified.")
            exit(1)

    simu.setControlVariable(control_var)

    simu.setDefaultParameters(unit)
    simu.readParameterFile(in_args['params_file'])
    simu.p.cross_shear = True
    simu.tagStrainParameters()

    if control_var == "rate":
        simu.setupNonDimensionalization(primary_rate, unit)
    else:
        simu.setupNonDimensionalization(shear_stress, unit)

    simu.assertParameterCompatibility()
    simu.resolveTimeOrStrainParameters()
    if in_args['binary_conf']:
        np_np_fixed = simu.get_np_Binary(in_args['config_file'])
        is2d = simu.isTwoDimensionBinary(in_args['config_file'])
    else:
        np_np_fixed = simu.get_np(in_args['config_file'])
        is2d = simu.isTwoDimension(in_args['config_file'])

    system = simu.getSys()
    system.set_np(np_np_fixed[0])
    if np_np_fixed[1] > 0:
        simu.p.np_fixed = np_np_fixed[1]

    system.setupSystemPreConfiguration("rate", is2d)
    if in_args['binary_conf']:
        simu.importConfigurationBinary(in_args['config_file'])
    else:
        simu.importConfiguration(in_args['config_file'])

    system.setupSystemPostConfiguration()
    simu.simu_name = getSimuName(in_args, control_var)
    simu.openOutputFiles(simu.simu_name)
    echoInput(in_args, simu)
    in_args['rate_primary'] = system.get_shear_rate()
    print("Simulation setup [ok]")


def calcShearRateAndDirection(rate_primary, rate_max_OSP, amplitude_OSP, time):
    omega = rate_max_OSP/amplitude_OSP
    rate_OSP = rate_max_OSP*np.cos(omega*time)
    total_rate = np.sqrt(rate_OSP**2 + rate_primary**2)

    theta = np.arctan2(rate_OSP, rate_primary)
    return total_rate, theta


def outputData(tk, simu, binconf_counter):
    system = simu.getSys()
    output_events =\
        tk.getElapsedClocks(system.get_time(),
                            np.abs(system.get_curvilinear_strain()))
    simu.generateOutput(output_events, binconf_counter)


def updateShearDirection(system, in_args):
    rate_OSP = in_args['rate_OSP_max_ratio']*in_args['rate_primary']
    rate, theta =\
        calcShearRateAndDirection(in_args['rate_primary'],
                                  rate_OSP,
                                  in_args['amplitude_OSP'],
                                  system.get_time())
    system.set_shear_rate(rate)

    system.setShearDirection(theta)


def halfPeriodStrain(rate_ratio_OSP, amplitude_OSP):
    strain = special.ellipe(rate_ratio_OSP**2/(1+rate_ratio_OSP**2))
    strain *= 2*np.sqrt(1+rate_ratio_OSP**2)*amplitude_OSP/rate_ratio_OSP
    return strain


def weaving_simu(in_args):
    simu = lfdem.Simulation()
    system = simu.getSys()

    setupSimulation(in_args, simu)

    sine_arch_strain = halfPeriodStrain(in_args['rate_OSP_max_ratio'],
                                        in_args['amplitude_OSP'])
    simu.p.time_interval_output_data = 0.01*sine_arch_strain
    simu.p.time_interval_output_config = 0.5*sine_arch_strain
    tk = simu.initTimeKeeper()
    binconf_counter = 0
    while simu.keepRunning():
        updateShearDirection(system, in_args)
        simu.timeEvolutionUntilNextOutput(tk)
        outputData(tk, simu, binconf_counter)
        simu.printProgress()
    print("Time evolution done")

if __name__ == '__main__':
    p = getArgParser()
    args = vars(p.parse_args())
    try:
        weaving_simu(args)
    except RuntimeError as e:
        print(e)
