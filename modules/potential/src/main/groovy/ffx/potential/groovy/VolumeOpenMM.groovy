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

import static java.lang.String.format
import static java.lang.String.format
import static java.lang.String.format

import com.sun.jna.ptr.PointerByReference

import static org.apache.commons.math3.util.FastMath.pow

import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KJPerKcal
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_NmPerAngstrom
import static edu.uiowa.jopenmm.GKNPOpenMMLibrary.OpenMM_GKNPForce_addParticle
import static edu.uiowa.jopenmm.GKNPOpenMMLibrary.OpenMM_GKNPForce_create
import static edu.uiowa.jopenmm.GKNPOpenMMLibrary.OpenMM_GKNPForce_setCutoffDistance
import static edu.uiowa.jopenmm.GKNPOpenMMLibrary.OpenMM_GKNPForce_setNonbondedMethod
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setForceGroup

import ffx.potential.ForceFieldEnergyOpenMM
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.nonbonded.GeneralizedKirkwood
import ffx.potential.parameters.ForceField
import static ffx.potential.nonbonded.GeneralizedKirkwood.DEFAULT_GAUSSVOL_RADII_OFFSET

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

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

    /**
     * Execute the script.
     */
    VolumeOpenMM run() {

        if (!init()) {
            return this
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.openAll(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        String filename = activeAssembly.getFile().getAbsolutePath()
        logger.info(" Running VolumeOpenMM on " + filename)

        // Compute the energy before adding the GKNP force.
        Atom[] atoms = activeAssembly.getAtomArray()
        ForceFieldEnergyOpenMM forceFieldEnergyOpenMM = (ForceFieldEnergyOpenMM) activeAssembly.getPotentialEnergy()
        double[] x = new double[atoms.length * 3]
        forceFieldEnergyOpenMM.getCoordinates(x)
        double initialEnergy = forceFieldEnergyOpenMM.energy(x, true)

        PointerByReference gknpForce = OpenMM_GKNPForce_create()
        OpenMM_GKNPForce_setNonbondedMethod(gknpForce, 0)
        OpenMM_GKNPForce_setCutoffDistance(gknpForce, 1.0)

        ForceField forceField = activeAssembly.getForceField()
        double probe = forceField.getDouble("PROBE_RADIUS", DEFAULT_GAUSSVOL_RADII_OFFSET)
        double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0)
        double surfTen = forceField.getDouble("SURFACE_TENSION", GeneralizedKirkwood.DEFAULT_CAVDISP_SURFACE_TENSION)
        double gamma = surfTen * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom)
        double alpha = 0.0
        double charge = 0.0
        for (Atom atom : atoms) {
            int isHydrogen = (!atom.isHydrogen()) ? 0 : 1
            if (includeHydrogen) {
                isHydrogen = 0
            }
            double radii = (atom.getVDWType().radius / 2.0)
            if (sigma) {
                radii *= rminToSigma
            }
            radii = (radii + probe) * OpenMM_NmPerAngstrom
            OpenMM_GKNPForce_addParticle(gknpForce, radii, gamma, alpha, charge, isHydrogen)
        }
        int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
        OpenMM_Force_setForceGroup(gknpForce, forceGroup)
        forceFieldEnergyOpenMM.addForce(gknpForce)
        forceFieldEnergyOpenMM.reinitContext()

        // Compute the energy after adding the GKNP force.
        double energy = forceFieldEnergyOpenMM.energy(x, true)
        double de = energy - initialEnergy

        logger.info(format("\n Surface Area:        %8.4f (Ang^2)", de/surfTen))
        logger.info(format(" Surface Tension:     %8.4f (kcal/mol/Ang^2)", surfTen))
        logger.info(format(" Surface Area Energy: %8.4f (kcal/mol)", de))

        return this
    }
}
