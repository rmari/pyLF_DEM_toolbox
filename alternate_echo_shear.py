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
    parser.add_argument('-g', '--shear_strain', type=float, required=True)
    parser.add_argument('-p', '--oscillation_nb', type=int, required=True)
    parser.add_argument('-r', '--rate_primary',
                        type=str, required=False, default="1h")
    parser.add_argument('-t', '--time_ratio', type=float)

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
    simu_name += "_echo_and_shear"\
                 + "_amplitude"+str(in_args['amplitude'])\
                 + "_"+str(in_args['oscillation_nb'])+"periods"\
                 + "_"+str(in_args['shear_strain'])+"straight"
    if in_args['time_ratio'] is not None:
        simu_name += "_rate"+str(in_args['rate_primary'])
        simu_name += "_ratio"+str(in_args['time_ratio'])
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

    print("Alternate Echo/Shear simulation.")
    print("Based on pyLFDEM version", simu.gitVersion())

    simu.force_to_run = in_args['overwrite']

    if in_args['time_ratio'] is not None\
            and in_args['rate_primary'] == "1h":
        print(" time ratio useless in rate independent rheology.",
              file=sys.stderr)
        sys.exit(1)
    if in_args['time_ratio'] is None\
            and in_args['rate_primary'] != "1h":
        print("needs time ratio in rate dependent rheology.",
              file=sys.stderr)
        sys.exit(1)

    primary_rate, unit = getSuffix(in_args['rate_primary'])
    control_var = "rate"
    simu.setControlVariable(control_var)

    simu.setDefaultParameters(unit)
    simu.readParameterFile(in_args['params_file'])
    simu.p.cross_shear = True
    simu.tagStrainParameters()

    simu.setupNonDimensionalization(primary_rate, unit)

    simu.assertParameterCompatibility()
    simu.resolveTimeOrStrainParameters()

    importConf(simu, in_args)

    system = simu.getSys()

    system.setupSystemPostConfiguration()
    simu.simu_name = getSimuName(in_args, control_var)
    simu.openOutputFiles(simu.simu_name)
    echoInput(in_args, simu)

    print("Simulation setup [ok]")


def outputData(tk, simu, binconf_counter):
    system = simu.getSys()
    output_events =\
        tk.getElapsedClocks(system.get_time(),
                            np.abs(system.get_curvilinear_strain()))
    simu.generateOutput(output_events, binconf_counter)


def getShearDirection(events, shear_dir):
    if shear_dir != 0:  # we are in echo mode
        if "reverse" in events:
            shear_dir = -shear_dir
        if "shear_on" in events:
            shear_dir = 0
    else:
        if "echo_on" in events:
            shear_dir = np.pi/2
    return shear_dir


def alternate_simu(in_args):
    simu = lfdem.Simulation()
    system = simu.getSys()

    setupSimulation(in_args, simu)

    in_args['amplitude'] = abs(in_args['amplitude'])

    print("[Note] Overriding time_interval_output_* in user parameter file.")
    arch_ratio_data = 10
    arch_ratio_par = 1
    simu.p.time_interval_output_data = 2*in_args['amplitude']/arch_ratio_data
    simu.p.time_interval_output_config = 2*in_args['amplitude']/arch_ratio_par

    tk = lfdem.TimeKeeper()
    tk.addClock("data",
                lfdem.LinearClock(simu.p.time_interval_output_data, True))
    tk.addClock("config",
                lfdem.LinearClock(simu.p.time_interval_output_config, True))

    n = in_args["oscillation_nb"]
    gamma_OSP = in_args["amplitude"]
    gamma_para = in_args["shear_strain"]
    echoing_strain = 4*n*gamma_OSP
    alternating_period = echoing_strain + gamma_para
    tk.addClock("echo_on",
                lfdem.LinearClock(alternating_period,
                                  True))
    tk.addClock("shear_on",
                lfdem.LinearClock(echoing_strain,
                                  alternating_period,
                                  True))
    tk.addClock("reverse",
                lfdem.LinearClock(gamma_OSP,
                                  2*gamma_OSP,
                                  True))
    if in_args["time_ratio"] is not None:
        echo_ratio = in_args["time_ratio"]
        avg_rate = system.get_shear_rate()
        para_rate = avg_rate/(1-echo_ratio)
        echoing_rate = 4*n*gamma_OSP*avg_rate/gamma_para/echo_ratio
    else:
        para_rate = system.get_shear_rate()
        echoing_rate = system.get_shear_rate()
    binconf_counter = 0
    shear_dir = np.pi/2
    system.setShearDirection(shear_dir)

    while simu.keepRunning():
        simu.timeEvolutionUntilNextOutput(tk)
        events =\
            tk.getElapsedClocks(system.get_time(),
                                np.abs(system.get_curvilinear_strain()))
        simu.generateOutput(events, binconf_counter)
        shear_dir = getShearDirection(events, shear_dir)
        system.setShearDirection(shear_dir)

        if "shear_on" in events:
            tk.removeClock("reverse")
            system.set_shear_rate(para_rate)
        if "echo_on" in events:
            current_strain = system.get_curvilinear_strain()
            reverse_clock =\
                lfdem.LinearClock(current_strain+gamma_OSP,
                                  2*gamma_OSP,
                                  True)
            tk.addClock("reverse", reverse_clock)
            system.set_shear_rate(echoing_rate)
        simu.printProgress()

    print("Time evolution done")


if __name__ == '__main__':
    p = getArgParser()
    args = vars(p.parse_args())
    try:
        alternate_simu(args)
    except RuntimeError as e:
        print(e)
