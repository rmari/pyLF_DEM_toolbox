#!/usr/bin/env python3

import argparse
import numpy as np
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
    parser.add_argument('-a', '--amplitude', type=float, required=True)

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
    simu_name += "_echo"\
                 + "_amplitude"+str(in_args['amplitude'])
    return simu_name


def getSuffix(value_str):
    split = re.search('[a-zA-Z]', value_str).start()
    return float(value_str[:split]), value_str[split:]


def importConf(simu, in_args):
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


def setupSimulation(in_args, simu):
    """ Setup for echo simulations """

    print("Echo simulation.")
    print("Based on pyLFDEM version", simu.gitVersion())

    simu.force_to_run = in_args['overwrite']

    simu.setControlVariable("rate")

    simu.setDefaultParameters("h")
    simu.readParameterFile(in_args['params_file'])
    simu.p.cross_shear = True
    simu.tagStrainParameters()

    simu.setupNonDimensionalization(1, "h")

    simu.assertParameterCompatibility()
    simu.resolveTimeOrStrainParameters()

    importConf(simu, in_args)

    system = simu.getSys()

    system.setupSystemPostConfiguration()
    simu.simu_name = getSimuName(in_args, "rate")
    simu.openOutputFiles(simu.simu_name)
    echoInput(in_args, simu)

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


def echo_simu(in_args):
    simu = lfdem.Simulation()
    system = simu.getSys()

    setupSimulation(in_args, simu)

    in_args['amplitude'] = abs(in_args['amplitude'])

    print("[Note] Overriding time_interval_output_* in user parameter file.")
    arch_ratio_data = 200
    arch_ratio_par = 10
    simu.p.time_interval_output_data = 2*in_args['amplitude']/arch_ratio_data
    simu.p.time_interval_output_config = 2*in_args['amplitude']/arch_ratio_par

    tk = lfdem.TimeKeeper()
    tk.addClock("data",
                lfdem.LinearClock(simu.p.time_interval_output_data, True))
    tk.addClock("data",
                lfdem.LinearClock(simu.p.time_interval_output_config, True))
    tk.addClock("reverse",
                lfdem.LinearClock(in_args["amplitude"],
                                  2*in_args["amplitude"],
                                  True))
    print(in_args["amplitude"])
    binconf_counter = 0
    shear_dir = 0
    system.setShearDirection(shear_dir)

    while simu.keepRunning():
        simu.timeEvolutionUntilNextOutput(tk)
        output_events =\
            tk.getElapsedClocks(system.get_time(),
                                np.abs(system.get_curvilinear_strain()))
        simu.generateOutput(output_events, binconf_counter)
        if "reverse" in output_events:
            shear_dir = np.pi - shear_dir
            system.setShearDirection(shear_dir)
        simu.printProgress()

    print("Time evolution done")

if __name__ == '__main__':
    p = getArgParser()
    args = vars(p.parse_args())
    try:
        echo_simu(args)
    except RuntimeError as e:
        print(e)
