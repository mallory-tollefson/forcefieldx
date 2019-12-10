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
import edu.uiowa.jopenmm.OpenMMUtils
import edu.uiowa.jopenmm.OpenMM_Vec3
import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.bonded.Residue
import ffx.potential.bonded.ResidueEnumerations
import org.apache.commons.lang3.ObjectUtils
import picocli.CommandLine

import java.awt.Point
import java.lang.reflect.Array
import java.util.logging.Level

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setForceGroup
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_append
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_create
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_destroy
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_append
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_create
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
import static org.apache.commons.math3.util.FastMath.sqrt

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
     * -y or --includeHydrogen leaves in hydrogen when calculating the overlap tree.
     */
    @CommandLine.Option(names = ['-y', '--includeHydrogen'], paramLabel = "false",
            description = "Include Hydrogen in calculation of overlaps and volumes")
    private boolean includeHydrogen = false

    /**
     * -s or --sigma Use sigma radii instead of Rmin.
     */
    @CommandLine.Option(names = ['-s', '--sigma'], paramLabel = "false",
            description = "Use sigma radii instead of Rmin")
    private boolean sigma = false
    
    /**
     * -v or --verbose enables printing out all energy components for multi-snapshot files (
     * the first snapshot is always printed verbosely).
     */
    @CommandLine.Option(names = ['-v', '--verbose'], paramLabel = "false",
            description = "Print out all components of volume of molecule and offset")
    private boolean verbose = false

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    public double energy = 0.0
    public ForceFieldEnergy forceFieldEnergy = null
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
            saveDir = new File(FilenameUtils.getFullPath(filename))
        }
        //EXAMPLE GKNP FORCE
        //PointerByReference system = OpenMMLibrary.OpenMM_System_create()
        //PointerByReference nonBondedForce = OpenMMLibrary.OpenMM_NonbondedForce_create()
        PointerByReference GKNPForce = GKNPOpenMMLibrary.OpenMM_GKNPForce_create()
        GKNPOpenMMLibrary.OpenMM_GKNPForce_setNonbondedMethod(GKNPForce, 0);
        GKNPOpenMMLibrary.OpenMM_GKNPForce_setCutoffDistance(GKNPForce, 1.0)
        //OpenMMLibrary.OpenMM_System_addForce(system, nonBondedForce)
        //OpenMMLibrary.OpenMM_System_addForce(system, GKNPForce)

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
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
        PointerByReference positions = OpenMM_Vec3Array_create(0)
        OpenMM_Vec3.ByValue pos = new OpenMM_Vec3.ByValue()
        Arrays.fill(gamma, 0.117*kcalmol2kjmol/(ang2nm*ang2nm))

        int index = 0
        for (Atom atom : atoms) {
            //OpenMMLibrary.OpenMM_System_addParticle(system, atom.getMass())

            isHydrogen[index] = atom.isHydrogen()
            if (includeHydrogen) {
                isHydrogen[index] = false
            }
            radii[index] = (atom.getVDWType().radius / 2.0)*ang2nm
            if (sigma) {
                radii[index] *= rminToSigma
            }
            pos.x = atom.getX()*ang2nm
            pos.y = atom.getY()*ang2nm
            pos.z = atom.getZ()*ang2nm
            OpenMM_Vec3Array_append(positions, pos)

            sigma_LJ = 2.0 * atom.getVDWType().radius
            double sij = sqrt(sigmaw*sigma_LJ)
            double eij = sqrt(epsilonw*epsilon_LJ)
            double alpha = - 16.0 * Math.PI * rho * eij * pow(sij,6) / 3.0
            //OpenMMLibrary.OpenMM_NonbondedForce_addParticle(nonBondedForce, 0, 0, 0)
            GKNPOpenMMLibrary.OpenMM_GKNPForce_addParticle(GKNPForce, radii[index], gamma[index], alpha,
                    atom.getCharge(), isHydrogen[index])
            index++
        }

        //PointerByReference verletIntegrator = OpenMMLibrary.OpenMM_VerletIntegrator_create(1.0)
        //PointerByReference platform = OpenMMLibrary.OpenMM_Platform_getPlatformByName("CUDA")

        //PointerByReference context = OpenMMLibrary.OpenMM_Context_create_2(system, verletIntegrator, platform)
        //OpenMMLibrary.OpenMM_Context_setPositions(context, positions)

        //PointerByReference state = OpenMMLibrary.OpenMM_Context_getState(context, OpenMM_State_Energy | OpenMM_State_Forces | OpenMM_State_Positions, 0)

        //double energy = 0
        //energy = OpenMMLibrary.OpenMM_State_getPotentialEnergy(state)

        // A forceGroup of 1 is for nonbonded forces; GKNP forces are nonbond-like.
        int forceGroup = 1;
        OpenMM_Force_setForceGroup(GKNPForce, forceGroup)

        logger.info("Adding GKNP Force to ForceFieldEnergyOpenMM")
        ForceFieldEnergyOpenMM forceFieldEnergyOpenMM = (ForceFieldEnergyOpenMM) forceFieldEnergy
        forceFieldEnergyOpenMM.addForce(GKNPForce)
        forceFieldEnergyOpenMM.reinitContext()

        energy = forceFieldEnergy.energy(x, true)

        return this
    }
}
