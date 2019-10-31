//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.cli;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.optimize.anneal.AnnealingSchedule;
import ffx.algorithms.optimize.anneal.FlatEndAnnealSchedule;
import ffx.algorithms.optimize.anneal.SimulatedAnnealing;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Option;

import java.io.File;
import java.util.logging.Logger;

/**
 * Represents command line options for scripts that utilize simulated annealing.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class AnnealOptions {

    private static final Logger logger = Logger.getLogger(AnnealOptions.class.getName());

    /**
     * -w or --windows Number of annealing windows (10).
     */
    @Option(names = {"-W", "--windows"}, paramLabel = "10",
            description = "Number of annealing windows.")
    int windows = 10;

    /**
     * --tmS or --temperingSchedule sets the schedule to be used.
     */
    @Option(names = {"--tmS", "--temperingSchedule"}, paramLabel = "EXP",
            description = "Tempering schedule: choose between EXP (exponential) or LINEAR")
    private String temperString = "EXP";

    /**
     * --tmB or --temperingBefore sets the number of annealing windows to hold
     * flat at the high temperature (in addition to normal windows).
     */
    @Option(names = {"--tmB", "--temperingBefore"}, paramLabel = "0",
            description = "Number of (annealing, not MD/MC) steps to remain at the high temperature")
    private int temperBefore = 0;

    /**
     * --tmA or --temperingAfter sets the number of annealing windows to hold
     * flat at the low temperature (in addition to normal windows).
     */
    @Option(names = {"--tmA", "--temperingAfter"}, paramLabel = "0",
            description = "Number of (annealing, not MD/MC) steps to remain at the low temperature")
    private int temperAfter = 0;

    /**
     * -l or --low Low temperature limit in degrees Kelvin (10.0).
     */
    @Option(names = {"--tl", "--temperatureLow"}, paramLabel = "10.0",
            description = "Low temperature limit (Kelvin).")
    double low = 10.0;

    /**
     * -u or --upper Upper temperature limit in degrees Kelvin (1000.0).
     */
    @Option(names = {"--tu", "--temperatureUpper"}, paramLabel = "1000.0",
            description = "High temperature limit (Kelvin).")
    double upper = 1000.0;

    @Option(names = {"--rv", "--reinitVelocities"}, paramLabel = "false",
            description = "Re-initialize velocities before each round of annealing.")
    boolean reinitV = false;

    /**
     * Constructs an AnnealingSchedule.
     *
     * @return An AnnealingSchedule.
     */
    public AnnealingSchedule getSchedule() {
        SimulatedAnnealing.Schedules sch = SimulatedAnnealing.Schedules.parse(temperString);
        AnnealingSchedule as = sch.generate(windows, low, upper);
        if (temperBefore > 0 || temperAfter > 0) {
            as = new FlatEndAnnealSchedule(as, low, upper, temperBefore, temperAfter);
        }
        return as;
    }

    /**
     * Creates a SimulatedAnnealing object.
     *
     * @param dynOpts   Dynamics options to use.
     * @param mola      MolecularAssembly
     * @param potential Potential
     * @param props     Properties
     * @param alist     AlgorithmListener
     * @return          SimulatedAnnealing
     */
    public SimulatedAnnealing createAnnealer(DynamicsOptions dynOpts, MolecularAssembly mola,
                                             Potential potential, CompositeConfiguration props,
                                             AlgorithmListener alist) {
        return createAnnealer(dynOpts, mola, potential, props, alist, null);
    }

    /**
     * Creates a SimulatedAnnealing object.
     *
     * @param dynOpts   Dynamics options to use.
     * @param mola      MolecularAssembly
     * @param potential Potential
     * @param props     Properties
     * @param alist     AlgorithmListener
     * @param dynFile   Dynamics restart file.
     * @return          SimulatedAnnealing
     */
    public SimulatedAnnealing createAnnealer(DynamicsOptions dynOpts, MolecularAssembly mola,
                                             Potential potential, CompositeConfiguration props,
                                             AlgorithmListener alist, File dynFile) {
        AnnealingSchedule schedule = getSchedule();
        double totNormLen = schedule.totalWindowLength();
        int totSteps = dynOpts.steps;
        int perWindowSteps = (int) (totSteps / totNormLen);

        if (totSteps > perWindowSteps * totNormLen) {
            ++perWindowSteps;
        }

        logger.info(String.format(" Each of %d simulated annealing windows will have %d " +
                "steps each for a total annealing duration of %d timesteps",
                schedule.getNumWindows(), perWindowSteps, perWindowSteps * schedule.getNumWindows()));

        return new SimulatedAnnealing(mola, potential, props, alist, dynOpts.thermostat,
                dynOpts.integrator, schedule, perWindowSteps, dynOpts.dt, reinitV, dynFile);
    }
}
