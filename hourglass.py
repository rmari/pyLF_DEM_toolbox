#!/usr/bin/env python3

import argparse
import numpy as np
import re
import pyLFDEM as lf
import pyLFDEM_helper as lfh
import sys


def getArgParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    parser.add_argument('params_file', type=str)
    parser.add_argument('-n', '--binary_conf', action='store_true')
    parser.add_argument('-f', '--overwrite', action='store_true')
    parser.add_argument('-r', '--rate', type=str, required=True)
    parser.add_argument('-a', '--amplitude',
                        type=float, required=True)
    parser.add_argument('-t', '--angle',
                        type=float, required=True)

    return parser


def getCallString():
    call_str = "Hourglass script called as:\n"+' '.join(sys.argv)+"\n\n"

    call_str += "Script file "+__file__+" is:\n\n"
    for line in open(__file__, "r"):
        call_str += line
    call_str += "\n"
    return call_str

def getSimuId(in_args):
    simu_name = "_hourglass"
    simu_name += "_amplitude"+str(in_args['amplitude'])\
                 + "_angle"+str(in_args['angle'])
    return simu_name

def getSuffix(value_str):
    split = re.search('[a-zA-Z]', value_str).start()
    return float(value_str[:split]), value_str[split:]


def getArchStrainData(in_args, n):
    aosp = in_args['amplitude_OSP']
    rosp = in_args['rate_OSP_max_ratio']

    x = np.linspace(0, 2*np.pi, 4*n+5)*aosp/rosp
    z = aosp*np.sin(rosp*x/aosp)
    dx = np.diff(x)
    dz = np.diff(z)
    thetas = np.arctan2(dz, dx)
    strain_steps = np.sqrt(dx**2+dz**2)

    rate_primary = in_args['rate_primary']
    rate_OSP = rosp*rate_primary*np.cos(rosp*x/aosp)
    shear_rates = np.sqrt(rate_OSP**2 + rate_primary**2)
    return shear_rates, thetas, strain_steps


def updateDirection(in_args, tk, events, simu):
    strain = simu.getSys().get_cumulated_strain()

    diagonal_amplitude = in_args["amplitude"]
    transverse_amplitude = 2*np.sin(in_args["angle"])*in_args["amplitude"]
    cycle_amplitude = 4*diagonal_amplitude + 2*transverse_amplitude
    if "end_diagonal1" in events:
        tk.removeClock()
        tk.addClock("data",
                    lf.LinearClock(transverse_amplitude/10, True))
        tk.addClock("config",
                    lf.LinearClock(transverse_amplitude, True))
        tk.addClock("end_transverse1",
                    lf.LinearClock(transverse_amplitude, True))
        simu.getSys().p.theta_shear = -90     
        simu.setupFlow()

    if "end_transverse1" in events:
        tk.removeClock()
        tk.addClock("data",
                    lf.LinearClock(diagonal_amplitude/10, True))
        tk.addClock("config",
                    lf.LinearClock(diagonal_amplitude, True))
        tk.addClock("end_diagonal2",
                    lf.LinearClock(2*diagonal_amplitude, True))
        simu.getSys().p.theta_shear = -in_args["angle"]     
        simu.setupFlow()

    if "end_diagonal2" in events:
        tk.removeClock()
        tk.addClock("data",
                    lf.LinearClock(transverse_amplitude/10, True))
        tk.addClock("config",
                    lf.LinearClock(transverse_amplitude, True))
        tk.addClock("end_transverse1",
                    lf.LinearClock(transverse_amplitude, True))
        simu.getSys().p.theta_shear = in_args["angle"]   
        simu.setupFlow()

    if "end_transverse2" in events:
        tk.removeClock()
        tk.addClock("data",
                    lf.LinearClock(diagonal_amplitude/10, True))
        tk.addClock("config",
                    lf.LinearClock(diagonal_amplitude, True))
        tk.addClock("end_diagonal1",
                    lf.LinearClock(2*diagonal_amplitude, True))
        simu.getSys().p.theta_shear = in_args["angle"]     
        simu.setupFlow()

    return tk

def initTimeKeeper(in_args):
    print("[Note] Overriding time_interval_output_* in user parameter file.")

    tk = lf.TimeKeeper()
    diagonal_amplitude = in_args["amplitude"]
    transverse_amplitude = 2*np.sin(in_args["angle"])*in_args["amplitude"]
    cycle_amplitude = 4*diagonal_amplitude + 2*transverse_amplitude
    tk.addClock("data",
                lf.LinearClock(diagonal_amplitude/10, True))
    tk.addClock("config",
                lf.LinearClock(diagonal_amplitude, True))
    tk.addClock("end_diagonal1",
                lf.LinearClock(diagonal_amplitude, True))
    return tk

def hourglass_shear(in_args):
    simu = lfh.setupBasicSimu(ctrl_str="r "+in_args['rate'], 
                              conf_file=in_args['config_file'], 
                              binary_conf=in_args['binary_conf'], 
                              params_file=in_args['params_file'], 
                              call_str=getCallString(),                              
                              identifier=getSimuId(in_args),
                              overwrite=in_args['overwrite']) 
    
    binconf_counter = 0
    simu.generateOutput(lf.SetString(["data", "config"]), binconf_counter)

    system = simu.getSys()
    tk = initTimeKeeper(in_args)
    while simu.keepRunning():
        simu.timeEvolutionUntilNextOutput(tk)
        events = tk.getElapsedClocks(system.get_time(), system.get_cumulated_strain())
        simu.generateOutput(events, binconf_counter)
        simu.printProgress()
        tk = updateDirection(in_args, tk, events, simu)

    print("Time evolution done")


if __name__ == '__main__':
    p = getArgParser()
    args = vars(p.parse_args())
    try:
        hourglass_shear(args)
    except RuntimeError as e:
        print(e)
