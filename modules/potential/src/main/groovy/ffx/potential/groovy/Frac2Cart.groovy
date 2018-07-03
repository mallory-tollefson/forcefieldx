package ffx.potential.groovy

import org.apache.commons.io.FilenameUtils

import ffx.crystal.Crystal
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The Frac2Cart script converts from Fractional to Cartesian coordinates.
 * <br>
 * Usage:
 * <br>
 * ffxc Frac2Cart &lt;filename&gt;
 */
@Command(description = " Convert fractional to Cartesian coordinates.", name = "ffxc Frac2Cart")
class Frac2Cart extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    /**
     * Execute the script.
     */
    @Override
    Frac2Cart run() {

        if (!init()) {
            return
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            assemblies = [activeAssembly]
        }

        String filename = activeAssembly.getFile().getAbsolutePath()
        logger.info("\n Converting from fractional to Cartesian coordinates for " + filename)

        // Loop over each system.
        for (int i = 0; i < assemblies.length; i++) {
            def system = assemblies[i]
            Crystal crystal = system.getCrystal().getUnitCell()

            List<Atom> atoms = system.getAtomList()
            double[] frac = new double[3]
            double[] cart = new double[3]

            for (Atom atom in atoms) {
                atom.getXYZ(frac)
                crystal.toCartesianCoordinates(frac, cart)
                atom.moveTo(cart)
            }
        }

        String ext = FilenameUtils.getExtension(filename)
        filename = FilenameUtils.removeExtension(filename)
        if (ext.toUpperCase().contains("XYZ")) {
            potentialFunctions.saveAsXYZ(assemblies[0], new File(filename + ".xyz"))
        } else {
            potentialFunctions.saveAsPDB(assemblies, new File(filename + ".pdb"))
        }

        return this
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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