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

import edu.uiowa.jopenmm.OpenMMLibrary
import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.bonded.Residue
import ffx.potential.bonded.ResidueEnumerations
import org.apache.commons.lang3.ObjectUtils

import java.awt.Point
import java.lang.reflect.Array
import java.util.logging.Level

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setForceGroup
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_append
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_create
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_destroy
import static java.lang.String.format

import com.google.common.collect.MinMaxPriorityQueue
import com.sun.jna.ptr.PointerByReference

import org.apache.commons.io.FilenameUtils
import ffx.numerics.Potential
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter

import edu.uiowa.jopenmm.GKNPOpenMMLibrary;

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.math3.util.FastMath.pow

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Energy &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy.", name = "ffxc VolumeOpenMM")
class VolumeOpenMM extends PotentialScript {
    /**
     * -v or --verbose enables printing out all energy components for multi-snapshot files (
     * the first snapshot is always printed verbosely).
     */
    @Option(names = ['-v', '--verbose'], paramLabel = "false",
            description = "Print out all energy components for multi-snapshot files")
    private boolean verbose = false;

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private List<String> filenames = null


    public double energy = 0.0
    public ForceFieldEnergy forceFieldEnergy = null


    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    private AssemblyState assemblyState = null

    private class StateContainer implements Comparable<StateContainer> {
        private final AssemblyState state
        private final double e

        StateContainer(AssemblyState state, double e) {
            this.state = state
            this.e = e
        }

        AssemblyState getState() {
            return state
        }

        double getEnergy() {
            return e
        }

        @Override
        int compareTo(StateContainer o) {
            return Double.compare(e, o.getEnergy())
        }
    }

    /**
     * Execute the script.
     */
    VolumeOpenMM run() {

        if (!init()) {
            return this
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        String filename = activeAssembly.getFile().getAbsolutePath()
        logger.info(" Running VolumeOpenMM on " + filename)

        File saveDir = baseDir
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        //EXAMPLE GKNP FORCE
        PointerByReference system = OpenMMLibrary.OpenMM_System_create()
        PointerByReference nonBondedForce = OpenMMLibrary.OpenMM_NonbondedForce_create()
        PointerByReference GKNPForce = GKNPOpenMMLibrary.OpenMM_GKNPForce_create()
        GKNPOpenMMLibrary.OpenMM_GKNPForce_setNonbondedMethod(GKNPForce, 0);
        GKNPOpenMMLibrary.OpenMM_GKNPForce_setCutoffDistance(GKNPForce, 1.0)
        OpenMMLibrary.OpenMM_System_addForce(system, nonBondedForce)
        OpenMMLibrary.OpenMM_System_addForce(system, GKNPForce)

        Atom[] atoms = activeAssembly.getAtomArray()
        int nAtoms = atoms.length

        //Constants for unit/energy conversions
        double solventPressure = 0.11337
        double volumeOffsetVdwToSEV = 27.939
        double surfaceAreaOffsetVdwToSASA = 46.111
        double surfaceTension = 0.16
        double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0)
        double ang2nm = 0.1
        double kcalmol2kjmol = 4.184
        double sigmaw = 3.15365*ang2nm /* LJ sigma of TIP4P water oxygen */
        double epsilonw = 0.155*kcalmol2kjmol       /* LJ epsilon of TIP4P water oxygen */
        double rho = 0.033428/pow(ang2nm,3)   /* water number density */
        double epsilon_LJ = 0.155*kcalmol2kjmol
        double sigma_LJ

        // Input
        boolean[] isHydrogen = new boolean[nAtoms]
        double[] radii = new double[nAtoms]
        double[] volume = new double[nAtoms]
        double[] gamma = new double[nAtoms]
        double[][] positions = new double[nAtoms][3]

        Arrays.fill(gamma, 1.0)
        double fourThirdsPI = 4.0 / 3.0 * Math.PI
        int index = 0
        for (Atom atom : atoms) {
            OpenMMLibrary.OpenMM_System_addParticle(system, atom.getMass())
            isHydrogen[index] = atom.isHydrogen()
            if (includeHydrogen) {
                isHydrogen[index] = false
            }
            radii[index] = atom.getVDWType().radius / 2.0
            if (sigma) {
                radii[index] *= rminToSigma;
            }
            radii[index] += probe
            volume[index] = fourThirdsPI * pow(radii[index], 3)
            positions[index][0] = atom.getX()
            positions[index][1] = atom.getY()
            positions[index][2] = atom.getZ()
            index++
        }




        // A forceGroup of 0 is for bonded forces; cavitation forces are non-bond-like.
        int forceGroup = 0;
        OpenMM_Force_setForceGroup(GKNPForce, forceGroup)

        logger.info("Adding GKNP Force to ForceFieldEnergyOpenMM")
        ForceFieldEnergyOpenMM forceFieldEnergyOpenMM = (ForceFieldEnergyOpenMM) forceFieldEnergy
        forceFieldEnergyOpenMM.addForce(GKNPForce)
        forceFieldEnergyOpenMM.reinitContext()

        energy = forceFieldEnergy.energy(x, true)

        return this
    }
}
