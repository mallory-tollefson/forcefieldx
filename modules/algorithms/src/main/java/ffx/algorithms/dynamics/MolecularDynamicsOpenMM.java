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
package ffx.algorithms.dynamics;

import java.io.File;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setVelocitiesToTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPotentialEnergy;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.parsers.DYNFilter;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.NS2SEC;
import static ffx.utilities.Constants.kB;

/**
 * Runs Molecular Dynamics using OpenMM implementation
 *
 * @author Hernan V. Bernabe
 */
public class MolecularDynamicsOpenMM extends MolecularDynamics {

    private static final Logger logger = Logger.getLogger(MolecularDynamicsOpenMM.class.getName());

    /**
     * OpenMM ForceFieldEnergy.
     */
    private ForceFieldEnergyOpenMM forceFieldEnergyOpenMM;
    /**
     * OpenMM Context.
     */
    private PointerByReference context;
    /**
     * OpenMM Integrator.
     */
    private PointerByReference integrator;
    /**
     * Number of OpenMM Particles (i.e. the number of FFX atoms).
     */
    private int numParticles;
    /**
     * Number of OpenMM Degrees of Freedom.
     */
    private int dof;
    /**
     * Integrator Type.
     */
    private final IntegratorEnum integratorType;
    /**
     * Thermostat Type.
     */
    private final ThermostatEnum thermostatType;
    /**
     * Integrator String.
     */
    private String integratorString;
    /**
     * Number of OpenMM MD steps per iteration.
     */
    private int intervalSteps;
    /**
     * Flag to indicate OpenMM MD interactions are running.
     */
    private boolean running;
    /**
     * Run time.
     */
    private long time;
    /**
     * Total energy at the end of a molecular dynamics move.
     */
    private double endTotalEnergy;
    /**
     * Total number of atoms in the system, typically used in for loops.
     */
    private int nAtoms;
    /**
     * Boolean used to signify molecular dynamics in the NPT
     * (isothermal-isobaric) ensemble.
     */
    private boolean constantPressure = false;
    /**
     * Double that holds the target pressure for the barostat under NPT
     * dynamics.
     */
    private double pressure;
    /**
     * Frequency of collisions for the barostat under NPT dynamics.
     */
    private int barostatFrequency;
    /**
     * Composite properties for the system used in the simulation.
     */
    private CompositeConfiguration properties;
    /**
     * Random number generator used to psuedo seed the OpenMM velocity generator
     * method.
     */
    private Random random;
    /**
     * Boolean to signify that we are updating the system (post-MD move).
     */
    private boolean update = false;

    /**
     * Constructs an MolecularDynamicsOpenMM object, to perform molecular
     * dynamics using native OpenMM routines, avoiding the cost of communicating
     * coordinates, gradients, and energies back and forth across the PCI bus.
     *
     * @param assembly               MolecularAssembly to operate on
     * @param forceFieldEnergyOpenMM ForceFieldEnergyOpenMM Potential. Cannot be any other type of Potential.
     * @param properties             Associated properties
     * @param listener               a {@link ffx.algorithms.AlgorithmListener} object.
     * @param thermostat             May have to be slightly modified for native OpenMM routines
     * @param integratorMD           May have to be slightly modified for native OpenMM routines
     */
    public MolecularDynamicsOpenMM(MolecularAssembly assembly, ForceFieldEnergyOpenMM forceFieldEnergyOpenMM,
                                   CompositeConfiguration properties, AlgorithmListener listener,
                                   ThermostatEnum thermostat, IntegratorEnum integratorMD) {
        super(assembly, forceFieldEnergyOpenMM, properties, listener, thermostat, integratorMD);

        // Initialization specific to MolecularDynamicsOpenMM
        running = false;
        this.forceFieldEnergyOpenMM = forceFieldEnergyOpenMM;
        numParticles = forceFieldEnergyOpenMM.getNumParticles();
        forceFieldEnergyOpenMM.addCOMMRemover(false);
        thermostatType = thermostat;
        integratorType = integratorMD;
        integratorToString(integratorType);

        this.properties = properties;

        random = new Random();
        if (properties.containsKey("velRandomSeed")) {
            random.setSeed(properties.getInt("velRandomSeed", 0));
        } else {
            random.setSeed(0);
        }
    }

    /**
     * Get the ForceFieldEnergyOpenMM instance used to run MD.
     *
     * @return a {@link ffx.potential.ForceFieldEnergyOpenMM} object.
     */
    public ForceFieldEnergyOpenMM getForceFieldEnergyOpenMM() {
        return forceFieldEnergyOpenMM;
    }

    /**
     * <p>
     * Setter for the field <code>pressure</code>.</p>
     *
     * @param pressure a double.
     */
    public void setPressure(double pressure) {
        this.pressure = pressure;
    }

    /**
     * <p>
     * Setter for the field <code>barostatFrequency</code>.</p>
     *
     * @param barostatFrequency a int.
     */
    public void setBarostatFrequency(int barostatFrequency) {
        this.barostatFrequency = barostatFrequency;
    }

    /**
     * <p>
     * setNPTDynamics.</p>
     */
    public void setNPTDynamics() {
        constantPressure = true;
    }

    /**
     * <p>
     * Setter for the field <code>intervalSteps</code>.</p>
     *
     * @param intervalSteps a int.
     */
    @Override
    public void setIntervalSteps(int intervalSteps) {
        this.intervalSteps = intervalSteps;
        logger.info(String.format(" Interval Steps set at %d", intervalSteps));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void init(int numSteps, double timeStep, double printInterval, double saveInterval,
                     String fileType, double restartFrequency, double temperature, boolean initVelocities, File dyn) {
        this.targetTemperature = temperature;
        this.dt = timeStep;
        this.printFrequency = (int) printInterval;
        this.restartFile = dyn;
        this.initVelocities = initVelocities;

        switch (thermostatType) {
            case BUSSI:
            case BERENDSEN:
                if (!integratorString.equalsIgnoreCase("LANGEVIN")) {
                    logger.log(basicLogging, String.format(" Replacing FFX thermostat %s with OpenMM Andersen thermostat", thermostatType));
                    forceFieldEnergyOpenMM.addAndersenThermostat(targetTemperature);
                    if (constantPressure) {
                        setMonteCarloBarostat(pressure, targetTemperature, barostatFrequency);
                    }
                } else {
                    logger.log(basicLogging, " Langevin/Stochastic dynamics already has temperature control, will not be adding thermostat!");
                }
                break;
            case ADIABATIC:
                if (integratorString.equalsIgnoreCase("LANGEVIN") && constantPressure) {
                    setMonteCarloBarostat(pressure, targetTemperature, barostatFrequency);
                }
            default:
                break;
            // No thermostat.
        }

        updateContext();

        // Convert the print interval to a print frequency.
        printFrequency = 100;
        printInterval *= 1000; // Time step is in fsec, so convert printInterval to fsec.
        if (printInterval >= dt) {
            printFrequency = (int) (printInterval / dt);
        }

        // Convert save interval to a save frequency.
        saveSnapshotFrequency = 1000;
        // saveInterval *= 1000; // Time step is in fsec, so convert saveInterval to fsec.
        if (saveInterval >= dt) {
            saveSnapshotFrequency = (int) (saveInterval / dt);
        }

        // Convert restart interval to frequency.
        saveRestartFileFrequency = 1000;
        // restartFrequency *= 1000; // Time step is in fsec, so converting restartFrequency
        if (restartFrequency >= dt) {
            saveRestartFileFrequency = (int) (restartFrequency / dt);
        }

        done = false;

        assemblyInfo();

        String firstFileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath());

        if (dyn == null) {
            restartFile = new File(firstFileName + ".dyn");
            loadRestart = false;
        } else {
            restartFile = dyn;
            loadRestart = true;
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        dof = forceFieldEnergyOpenMM.calculateDegreesOfFreedom();

        boolean setPositions = true;
        boolean setVelocities = true;

        if (!initialized) {
            if (loadRestart) {
                Crystal crystal = molecularAssembly.getCrystal();
                if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                } else {
                    logger.log(basicLogging, " Continuing from " + dyn.getAbsolutePath());
                    forceFieldEnergyOpenMM.setCrystal(crystal);
                    // Load positions into the main FFX data structure, move into primary unit cell, then load to OpenMM.
                    Atom[] atoms = molecularAssembly.getAtomArray();
                    double[] xyz = new double[3];
                    double[] vel = new double[3];
                    double[] acc = new double[3];
                    double[] accPrev = new double[3];
                    for (int i = 0; i < atoms.length; i++) {
                        Atom atom = atoms[i];
                        int i3 = i * 3;
                        for (int j = 0; j < 3; j++) {
                            xyz[j] = x[i3 + j];
                            vel[j] = v[i3 + j];
                            acc[j] = a[i3 + j];
                            accPrev[j] = aPrevious[i3 + j];
                        }
                        atom.setXYZ(xyz);
                        atom.setVelocity(vel);
                        atom.setAcceleration(acc);
                        atom.setPreviousAcceleration(accPrev);
                    }
                    molecularAssembly.moveAllIntoUnitCell();
                }
            } else {
                if (initVelocities) {
                    if (properties.containsKey("openMMInitVel")) {
                        int randomSeed = random.nextInt();
                        long ommSetVelTime = -System.nanoTime();
                        OpenMM_Context_setVelocitiesToTemperature(context, targetTemperature, randomSeed);
                        logger.log(intermediateLogging, format(" OpenMM set velocities to target temperature %f with random seed %d",
                                targetTemperature, randomSeed));
                        ommSetVelTime += System.nanoTime();
                        logger.log(intermediateLogging, format(" Set velocites with OpenMM in %6.3f", ommSetVelTime * NS2SEC));
                        setVelocities = false;
                    } else {
                        getThermostat().setQuiet(true);
                        getThermostat().maxwell(targetTemperature);
                        Atom[] atoms = molecularAssembly.getAtomArray();
                        double[] vel = new double[3];
                        long ffxSetVelTime = -System.nanoTime();
                        for (int i = 0; i < atoms.length; i++) {
                            Atom atom = atoms[i];
                            int i3 = i * 3;
                            for (int j = 0; j < 3; j++) {
                                vel[j] = v[i3 + j];
                            }
                            atom.setVelocity(vel);
                        }
                        ffxSetVelTime += System.nanoTime();
                        logger.fine(String.format(" Set velocities with FFX in %6.3f", ffxSetVelTime * NS2SEC));
                    }
                }
            }

            long setPosVel = -System.nanoTime();
            setOpenMMState(setPositions, setVelocities);
            setPosVel += System.nanoTime();

            logger.fine(String.format(" Set positions and velocities in %6.3f seconds", setPosVel * NS2SEC));

            // Call to retrieve the starting kinetic energy for the system.
            long retrieveEnergyTime = -System.nanoTime();

            forceFieldEnergyOpenMM.setLambda(forceFieldEnergyOpenMM.getLambda());

            getOpenMMEnergies();

            retrieveEnergyTime += System.nanoTime();

            logger.fine(String.format(" Retrieved energies in %6.3f seconds", retrieveEnergyTime * NS2SEC));

            startingKineticEnergy = currentKineticEnergy;
        }

        saveSnapshotAsPDB = true;
        if (fileType.equalsIgnoreCase("XYZ")) {
            saveSnapshotAsPDB = false;
        } else if (!fileType.equalsIgnoreCase("PDB")) {
            logger.warning(" Snapshot file type unrecognized; saving snaphshots as PDB.\n");
        }

        Atom[] atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        running = false;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Start sets up context, write out file name, restart file name, sets the
     * integrator and determines whether the simulation is starting out from a
     * previous molecular dynamics run (.dyn) or if the initial velocities are
     * determined by a Maxwell Boltzmann distribution. This method then calls
     * methods openMMUpdate and takeOpenMMSteps to run the molecular dynamics
     * simulation.
     */
    @Override
    public void dynamic(int numSteps, double timeStep, double printInterval, double saveInterval, double temperature, boolean initVelocities, File dyn) {

        long initTime = -System.nanoTime();
        init(numSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency, temperature, initVelocities, dyn);
        initTime += System.nanoTime();

        logger.fine(String.format("\n Initialized system in %6.3f sec.", initTime * NS2SEC));

        storeState();

        if (intervalSteps == 0 || intervalSteps > numSteps) {
            intervalSteps = numSteps;
        }
        running = true;

        // Update the time step in Picoseconds.
        OpenMM_Integrator_setStepSize(integrator, dt * 1.0e-3);

        int i = 0;
        time = -System.nanoTime();

        // Initial update from OpenMM.
        updateFromOpenMM(i, running);

        while (i < numSteps) {

            // Take MD steps in OpenMM.
            long takeStepsTime = -System.nanoTime();
            takeOpenMMSteps(intervalSteps);
            takeStepsTime += System.nanoTime();
            logger.fine(String.format("\n Took steps in %6.3f", takeStepsTime * NS2SEC));

            // Update the total step count.
            i += intervalSteps;

            update = true;
            long secondUpdateTime = -System.nanoTime();
            updateFromOpenMM(i, running);
            secondUpdateTime += System.nanoTime();

            logger.fine(String.format("\n Update finished in %6.3f", secondUpdateTime * NS2SEC));
            update = false;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getStartingTotalEnergy() {
        return startingTotalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getEndTotalEnergy() {
        return endTotalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTimeStep() {
        return dt;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getIntervalSteps() {
        return intervalSteps;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumAtoms() {
        return nAtoms;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setFileType(String fileType) {
        this.fileType = fileType;
    }

    /**
     * {@inheritDoc}
     * <p>
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }

    /**
     * {@inheritDoc}
     * <p>
     * UNSUPPORTED: MolecularDynamicsOpenMM is not presently capable of handling
     * extended system variables. Will throw an UnsupportedOperationException.
     */
    @Override
    public void detachExtendedSystem() {
        throw new UnsupportedOperationException(" MolecularDynamicsOpenMM does not support extended system variables!");
    }

    /**
     * <p>
     * updateContext.</p>
     */
    private void updateContext() {
        if (context == null) {
            forceFieldEnergyOpenMM.createContext(integratorString, dt, targetTemperature);
        } else {
            String currentIntegrator = forceFieldEnergyOpenMM.getIntegratorString();
            double currentTimeStp = forceFieldEnergyOpenMM.getTimeStep();
            double currentTemperature = forceFieldEnergyOpenMM.getTemperature();
            if (currentTemperature != targetTemperature || currentTimeStp != dt
                    || !currentIntegrator.equalsIgnoreCase(integratorString)) {
                forceFieldEnergyOpenMM.createContext(integratorString, dt, targetTemperature);
            }
        }
        context = forceFieldEnergyOpenMM.getContext();
        integrator = forceFieldEnergyOpenMM.getIntegrator();
    }

    /**
     * takeOpenMMSteps moves the simulation forward in time a user defined
     * number of steps and integrates the equations of motion for each step.
     * This method ensures that the algorithm reports back only when the time
     * interval (steps) specified by the user is completed.
     *
     * @param intervalSteps Number of MD steps to take.
     */
    private void takeOpenMMSteps(int intervalSteps) {
        OpenMM_Integrator_step(integrator, intervalSteps);
    }

    /**
     * Get OpenMM Energies.
     */
    private void getOpenMMEnergies() {
        context = forceFieldEnergyOpenMM.getContext();

        PointerByReference state = OpenMM_Context_getState(context, OpenMM_State_Energy, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);

        OpenMM_State_destroy(state);
    }

    /**
     * Get OpenMM Energies and Positions.
     */
    private void getOpenMMEnergiesAndPositions() {
        context = forceFieldEnergyOpenMM.getContext();

        int infoMask = OpenMM_State_Positions + OpenMM_State_Energy;
        PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);

        currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        currentKineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        currentTotalEnergy = currentPotentialEnergy + currentKineticEnergy;
        currentTemperature = 2.0 * currentKineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);

        PointerByReference positions = OpenMM_State_getPositions(state);
        forceFieldEnergyOpenMM.getOpenMMPositions(positions, numParticles * 3, x);
        forceFieldEnergyOpenMM.getPeriodicBoxVectors(state);

        OpenMM_State_destroy(state);
    }

    /**
     * Set to OpenMM positions and velocities.
     */
    private void setOpenMMState(boolean setPositions, boolean setVelocities) {
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;

        if (setPositions) {
            if (x == null || x.length < nAtoms * 3) {
                logger.severe(" Position vector has not been allocated.");
            }

            for (int i = 0; i < nAtoms; i++) {
                int index = i * 3;
                Atom atom = atoms[i];
                x[index] = atom.getX();
                x[index + 1] = atom.getY();
                x[index + 2] = atom.getZ();
            }
            forceFieldEnergyOpenMM.setOpenMMPositions(x, numberOfVariables);
            forceFieldEnergyOpenMM.setOpenMMPeriodicBoxVectors();
        }

        if (setVelocities) {
            double[] velocity = new double[3];
            for (int i = 0; i < nAtoms; i++) {
                int index = i * 3;
                Atom atom = atoms[i];
                atom.getVelocity(velocity);
                v[index] = velocity[0];
                v[index + 1] = velocity[1];
                v[index + 2] = velocity[2];
            }
            forceFieldEnergyOpenMM.setOpenMMVelocities(v, numberOfVariables);
        }
    }

    /**
     * updateFromOpenMM obtains the state of the simulation from OpenMM,
     * completes some logging, and saves restart files.
     *
     * @param i       Number of OpenMM MD rounds.
     * @param running True if OpenMM MD rounds have begun running.
     */
    private void updateFromOpenMM(int i, boolean running) {

        double priorPE = currentPotentialEnergy;

        if (update) {
            getOpenMMEnergiesAndPositions();
        }

        double defaultDeltaPEThresh = 1.0E6;
        detectAtypicalEnergy(priorPE, defaultDeltaPEThresh);

        if (running) {
            if (i == 0) {
                logger.log(basicLogging, format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
                logger.log(basicLogging, format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
                logger.log(basicLogging, format("  %8s %12.4f %12.4f %12.4f %8.2f",
                        "", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));
                startingKineticEnergy = currentKineticEnergy;
                startingTotalEnergy = currentTotalEnergy;
            } else if (i % printFrequency == 0) {
                double simTime = i * dt * 1.0e-3;
                time += System.nanoTime();
                logger.log(basicLogging, format(" %7.3e %12.4f %12.4f %12.4f %8.2f %8.2f",
                        simTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * NS2SEC));

                endTotalEnergy = currentTotalEnergy;

                time = -System.nanoTime();
            }

            if (saveSnapshotFrequency > 0 && i % (saveSnapshotFrequency * 1000) == 0 && i != 0) {
                for (AssemblyInfo ai : assemblies) {
                    if (ai.archiveFile != null && !saveSnapshotAsPDB) {
                        if (ai.xyzFilter.writeFile(ai.archiveFile, true)) {
                            logger.log(basicLogging, String.format(" Appended snap shot to %s", ai.archiveFile.getName()));
                        } else {
                            logger.warning(String.format(" Appending snap shot to %s failed", ai.archiveFile.getName()));
                        }
                    } else if (saveSnapshotAsPDB) {
                        if (ai.pdbFilter.writeFile(ai.pdbFile, false)) {
                            logger.log(basicLogging, String.format(" Wrote PDB file to %s", ai.pdbFile.getName()));
                        } else {
                            logger.warning(String.format(" Writing PDB file to %s failed.", ai.pdbFile.getName()));
                        }
                    }
                }
            }

            // Write out restart files every saveRestartFileFrequency steps.
            if (saveRestartFileFrequency > 0 && i % (saveRestartFileFrequency * 1000) == 0 && i != 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.log(basicLogging, format(" Wrote dynamics restart file to %s", restartFile.getName()));
                    if (constantPressure) {
                        Crystal crystal = molecularAssembly.getCrystal();
                        double currentDensity = crystal.getDensity(molecularAssembly.getTotalMass());
                        logger.info(format(" Density %6.3f (g/cc) with unit cell %s.",
                                currentDensity, crystal.toShortString()));
                    }
                } else {
                    logger.info(format(" Writing dynamics restart file to %s failed.", restartFile.getName()));
                }
            }
        }
    }

    /**
     * <p>
     * integratorToString.</p>
     *
     * @param integrator a {@link ffx.algorithms.dynamics.integrators.IntegratorEnum}
     *                   object.
     */
    private void integratorToString(IntegratorEnum integrator) {
        if (integrator == null) {
            integratorString = "VERLET";
            logger.info(" No specified integrator, will use Verlet");
        } else {
            switch (integratorType) {
                case STOCHASTIC:
                    integratorString = "LANGEVIN";
                    break;
                case VERLET:
                case VELOCITYVERLET:
                    integratorString = "VERLET";
                    break;
                case RESPA:
                    integratorString = "RESPA";
                    break;
                default:
                    integratorString = "VERLET";
                    logger.warning(String.format(" Integrator %s incompatible with "
                            + "OpenMM MD integration; defaulting to %s", integratorType, integratorString));
                    break;
            }
        }
    }

    /**
     * <p>
     * setMonteCarloBarostat.</p>
     *
     * @param pressure    a double.
     * @param temperature a double.
     * @param frequency   a int.
     */
    private void setMonteCarloBarostat(double pressure, double temperature, int frequency) {
        forceFieldEnergyOpenMM.addMonteCarloBarostat(pressure, temperature, frequency);
    }
}
