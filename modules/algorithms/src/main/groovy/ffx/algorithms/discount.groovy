/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */

// MOLECULAR & STOCHASTIC DYNAMICS

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.DiscountPh;
import ffx.algorithms.DiscountPh.Mode;
import ffx.algorithms.Protonate.DynamicsLauncher;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedVariable;
import ffx.potential.extended.TitrationESV;

// Number of molecular dynamics steps
int nSteps = 1000000;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 0.1;

// Temperature in degrees Kelvin.
double temperature = 298.15;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = null;

// Integrators [ BEEMAN, RESPA, STOCHASTIC, VELOCITYVERLET]
Integrators integrator = null;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Interval to write out restart file (psec)
double restartFrequency = 1000;

// File type of snapshots.
String fileType = "PDB";

// Monte-Carlo step frequencies for titration and rotamer moves.
int titrationFrequency = 10;
int titrationDuration = 1000;
int rotamerFrequencyRatio = 2;

// Simulation pH
double pH = 7.4;

// Single-residue titration option.
Character chainID = ' ';
int resID = -1;
List<String> resList = new ArrayList<>();
double window = 2.0;
boolean titrateTermini = false;
Mode mode = DiscountPh.Mode.USE_CURRENT;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc discount [options] <filename>');
cli.h(longOpt:'help', 'Print this message.');
//cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
//cli.i(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic / VELOCITYVERLET]');
cli.l(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
cli.s(longOpt:'restart', args:1, argName:'0.1', 'Interval to write out restart file (psec).');
cli.f(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');
//cli.ra(longOpt:'resAll', 'Titrate all residues.');
cli.rl(longOpt:'resList', args:1, 'Titrate a list of residues (eg A4.A8.B2.B34)');
//cli.rn(longOpt:'resName', args:1, 'Titrate a list of residue names (eg "LYS,TYR,HIS")');
//cli.rw(longOpt:'resWindow', args:1, 'Titrate all residues with intrinsic pKa within [arg] units of simulation pH.');
cli.pH(longOpt:'pH', args:1, argName:'7.4', 'Constant simulation pH.');
cli.mc(longOpt:'mcStepFreq', args:1, argName:'10', 'Number of steps between Monte-Carlo proton attempts.')
cli.mcr(longOpt:'rotamerStepFreq', args:1, argName:'0', 'Number of steps between Monte-Carlo rotamer attempts.')
cli.mcmd(longOpt:'mcRunTime', args:1, argName:'1000', 'Number of steps for which to run continuous proton dynamics during MC move.');
cli.a(longOpt:'mode', args:1, argName:'useCurrent', 'Controls starting ESV lambda values: [random, halfLambda, useCurrent]');
//cli.tt(longOpt:'titrateTermini', args:1, argName:'false', 'Titrate amino acid chain ends.');
def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

if ((options.rw && (options.ra || options.rl)) || (options.ra && options.rl)) {
    return cli.usage();
    logger.info(" Must specify one of the following: -ra, -rl, or -rw.");
}

if (!options.ra && !options.rl && !options.rw && !options.rn) {
    return cli.usage();
    logger.info(" Must specify one of the following: -ra, -rl, -rn, or -rw.");
}

if (options.a) {
    if ((options.a).equalsIgnoreCase("halfLambda")) {
        mode = DiscountPh.Mode.HALF_LAMBDA;
    } else if ((options.a).equalsIgnoreCase("useCurrent")) {
        mode = DiscountPh.Mode.USE_CURRENT;
    } else {
        mode = DiscountPh.Mode.valueOf(options.a);
    }
}

if (!options.pH) {
    return cli.usage();
    logger.info(" Must specify a solution pH.");
}

if (options.rl) {
    def tok = (options.rl).tokenize('.');
    for (String t : tok) {
        resList.add(t);
    }
}

if (options.tt) {
    titrateTermini = true;
    System.setProperty("cphmd-termini","true");
}

if (options.mcD) {
    dynamics = Boolean.parseBoolean(options.mcD);
}

if (options.rw) {
    window = Double.parseDouble(options.rw);
}

if (options.mc) {
    mcStepFrequency = Integer.parseInt(options.mc);
}

if (options.mcr) {
    rotamerStepFrequency = Integer.parseInt(options.mcr);
}

if (options.pH) {
    pH = Double.parseDouble(options.pH);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Write dyn interval in picoseconds
if (options.s) {
    restartFrequency = Double.parseDouble(options.s);
}

//
if (options.f) {
    fileType = options.f.toUpperCase();
}
// Load the time steps in femtoseconds.
if (options.d) {
    timeStep = Double.parseDouble(options.d);
}

// Report interval in picoseconds.
if (options.l) {
    printInterval = Double.parseDouble(options.l);
}

// Write snapshot interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

if (options.p) {
    System.setProperty("polarization", options.p);
}



// Thermostat.
if (options.b) {
    try {
        thermostat = Thermostats.valueOf(options.b.toUpperCase());
    } catch (Exception e) {
        thermostat = null;
    }
}

if (options.mcmd) {
    mcRunTime = Integer.parseInt(options.mcmd);
}

// Integrator.
if (options.i) {
    try {
        integrator = Integrators.valueOf(options.i.toUpperCase());
    } catch (Exception e) {
        integrator = null;
    }
}

System.setProperty("forcefield","AMOEBA_PROTEIN_2013");
System.setProperty("mpoleterm","false");
System.setProperty("pme-qi","true");
System.setProperty("cphmd-mode","USE_CURRENT");
// TODO: add flag which sets up the runs necessary to fit a model compound reference


/*
    String zeroReferenceEnergies = System.getProperty("cphmd-zeroReferences");  // [no-arg]
    String overrideFlag = System.getProperty("cphmd-override");                 // MCOverride.ACCEPT|REJECT|ONCE
    String beforeAfter = System.getProperty("cphmd-snapshots");                 // SnapshotsType.SEPARATE|INTERLEAVED|NONE
    String histidineMode = System.getProperty("cphmd-histidineMode");
*/

List<String> arguments = options.arguments();
String modelfilename = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    modelfilename = arguments.get(0);
    open(modelfilename);
} else if (active == null) {
    return cli.usage();
} else {
    modelfilename = active.getFile();
}

logger.info("\n Running molecular dynamics on " + modelfilename);

// Restart File
File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}

// create the MD object
MolecularDynamics molDyn = new MolecularDynamics(active, active.getPotentialEnergy(), active.getProperties(), sh, thermostat, integrator);
molDyn.setFileType(fileType);
molDyn.setRestartFrequency(restartFrequency);

// create the Monte-Carlo listener and connect it to the MD
DiscountPh cphmd = new DiscountPh(active, molDyn, pH, temperature, dyn);

// set residues to be titrated
if (options.ra) {
    cphmd.findTitrations();
} else if (options.rl) {
    cphmd.findTitrations(resList);
} else if (options.rw) {
    cphmd.findTitrations(pH, window);
} else if (options.rn) {
    cphmd.findTitrations(options.rn);
}

// finalize the ESV machinery and ready the MD (re-)launcher
//Object[] opt = new Object[9];
//opt[0] = new Integer(mcRunTime);    opt[1] = new Double(1.0);           opt[2] = new Double(1.0);
//opt[3] = new Double(0.0);           opt[4] = new Double(temperature);   opt[5] = new Boolean(false);
//opt[6] = fileType;                  opt[7] = new Double(0.0);           opt[8] = dyn;

// launch
cphmd.readyup();
//mcProt.launch(molDyn, nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities, fileType, restartFrequency);
cphmd.dynamic(nSteps, titrationFrequency, titrationDuration, timeStep, printInterval, saveInterval, temperature, initVelocities, fileType, restartFrequency);