#!/usr/bin/env python
#
# Make a succession of imposed shear stress and relaxation.
#

import numpy as np
import sys
import Simulation_LFDEM as lfdem


def setupSimulation(sim, conf_file, param_file, binary_conf,
                    dimensionlessnumber, input_scale):
    mysys = sim.getSys()

    sim.setupNonDimensionalization(dimensionlessnumber, input_scale)

    sim.assertParameterCompatibility()

    sim.resolveTimeOrStrainParameters()

    if binary_conf:
        mysys.set_np(sim.get_np_Binary(conf_file))
        is2d = sim.isTwoDimensionBinary(conf_file)
    else:
        mysys.set_np(sim.get_np(conf_file))
        is2d = sim.isTwoDimension(conf_file)

    mysys.setupSystemPreConfiguration("stress", is2d)

    if binary_conf:
        sim.importConfigurationBinary(conf_file)
    else:
        sim.importConfiguration(conf_file)

    mysys.setupSystemPostConfiguration()

    sim.prepareSimulationName(binary_conf, conf_file, param_file,
                              "hammer", dimensionlessnumber, input_scale)
    sim.openOutputFiles()
    in_args = ' '.join(sys.argv)
    in_files = lfdem.StringVector(2)
    in_files[0] = conf_file
    in_files[1] = param_file

    sim.echoInputFiles(in_args, in_files)


def hammer_simu():
    conf_file, param_file, stress, period = sys.argv[1:]
    stress = float(stress)
    period = float(period)
    is_binary = True

    mysim = lfdem.Simulation()

    output_times = np.linspace(0, period, 11)
    output_times = output_times[1:]
    scale = "kn"
    mysim.setControlVariable("stress")
    mysim.setDefaultParameters(scale)
    mysim.readParameterFile(param_file)
    mysim.p.out_data_particle = True
    mysim.p.out_data_interaction = True

    mysim.p.time_output_data = period
    mysim.p.time_output_config = period

    setupSimulation(mysim, conf_file, param_file, is_binary, stress, scale)
    mysim.setupEvents()
    mysys = mysim.getSys()
    running_stress = mysys.target_stress

    mysim.generateOutput(0, 0)
    while output_times[0] < 100*period:
        for t in output_times:
            mysys.timeEvolution("time", t)
            mysim.evaluateData()
            mysim.outputData()
            mysim.outputConfigurationBinary()
            print("time: ", mysys.get_time_in_simulation_units(),
                  " , strain: ", mysys.get_shear_strain())
        mysim.outputConfigurationData()
        output_times += period
        mysys.target_stress = running_stress - mysys.target_stress

if len(sys.argv) != 5:
    print("usage: ", sys.argv[0], "binconf params stress period")
    exit(1)

try:
    hammer_simu()
except RuntimeError as e:
    print(e)
