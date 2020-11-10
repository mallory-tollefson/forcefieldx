//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.potential.groovy.test

import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.GradientOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.extended.ExtendedSystem
import ffx.potential.extended.ExtendedVariable
import ffx.potential.extended.TitrationESV
import ffx.potential.extended.TitrationUtils
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import java.util.stream.IntStream

import static ffx.potential.extended.TitrationUtils.inactivateResidue
import static ffx.utilities.StringUtils.parseAtomRanges
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.sqrt
import static org.apache.commons.math3.util.FastMath.abs

/**
 * The Gradient script evaluates the consistency of the energy and gradient.
 * <br>
 * Usage:
 * <br>
 * ffxc test.Gradient [options] &lt;filename&gt;
 */
@Command(description = " Test the potential energy gradient.", name = "ffxc test.Gradient")
class PhGradient extends PotentialScript {

  @Mixin
  GradientOptions gradientOptions

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  /**
   * --pH or --constantPH Constant pH value for the test.
   */
  @CommandLine.Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
      description = 'Constant pH value for the test.')
  double pH = 7.4

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

  private ForceFieldEnergy energy
  //For unit test of endstate energies.
  public HashMap<String,double[]> endstateEnergyMap = new HashMap<String,double[]>()
  public int nFailures = 0
  public int nESVFailures = 0

  /**
   * Gradient constructor.
   */
  PhGradient() {
    this(new Binding())
  }

  /**
   * Gradient constructor.
   * @param binding The Groovy Binding to use.
   */
  PhGradient(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  PhGradient run() {

    if (!init()) {
      return this
    }

    TitrationUtils.initDiscountPreloadProperties()

    String modelFilename
    if (filenames != null && filenames.size() > 0) {
      MolecularAssembly molecularAssembly =
          TitrationUtils.openFullyProtonated(new File(filenames.get(0)))
      MolecularAssembly[] assemblies = [molecularAssembly]
      activeAssembly = assemblies[0]
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      modelFilename = activeAssembly.getFile().getAbsolutePath()
    }

    logger.info("\n Testing the atomic coordinate gradient of " + modelFilename + "\n")

    // Select all possible titrating residues.
    List<Residue> titrating = TitrationUtils.chooseTitratables(activeAssembly)
    energy = activeAssembly.getPotentialEnergy()

    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly)
    List<ExtendedVariable> titratingESVs = new ArrayList<>()
    esvSystem.setConstantPh(pH)
    for (Residue res : titrating) {
      ffx.potential.bonded.MultiResidue multi = TitrationUtils.titratingMultiresidueFactory(activeAssembly, res)
      TitrationESV esv = new TitrationESV(esvSystem, multi)
      titratingESVs.add(esv)
      for (Residue background : multi.getInactive()) {
        inactivateResidue(background)
      }
      esvSystem.addVariable(esv)
    }
    energy.attachExtendedSystem(esvSystem)
    logger.info(format(" Extended system with %d residues.", titratingESVs.size()));

    Atom[] atoms = activeAssembly.getAtomArray()
    int nAtoms = atoms.length

    // Apply atom selections
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    // Finite-difference step size in Angstroms.
    double step = gradientOptions.dx
    logger.info(" Finite-difference step size:\t" + step)

    // Print out the energy for each step.
    boolean print = gradientOptions.verbose
    logger.info(" Verbose printing:\t\t" + print)

    // Collect atoms to test.
    List<Integer> atomsToTest
    if (gradientOptions.gradientAtoms.equalsIgnoreCase("NONE")) {
      logger.info(" The gradient of no atoms will be evaluated.")
      return this
    } else if (gradientOptions.gradientAtoms.equalsIgnoreCase("ALL")) {
      logger.info(" Checking gradient for all active atoms.\n")
      atomsToTest = new ArrayList<>()
      IntStream.range(0, nAtoms).forEach(val -> atomsToTest.add(val))
    } else {
      atomsToTest = parseAtomRanges(" Gradient atoms", gradientOptions.gradientAtoms, nAtoms)
      logger.info(
          " Checking gradient for active atoms in the range: " + gradientOptions.gradientAtoms +
              "\n")
    }

    // Map selected atom numbers to active atom numbers
    Map<Integer, Integer> allToActive = new HashMap<>()
    int nActive = 0
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i]
      if (atom.isActive()) {
        allToActive.put(i, nActive)
        nActive++
      }
    }

    // Collect analytic gradient.
    double n = energy.getNumberOfVariables()
    double[] x = new double[n]
    double[] g = new double[n]
    energy.getCoordinates(x)
    energy.energyAndGradient(x, g)
    int index = 0
    double[][] allAnalytic = new double[nAtoms][3]
    for (Atom a : atoms) {
      a.getXYZGradient(allAnalytic[index++])
    }

    // Upper bound for a typical atomic gradient.
    double expGrad = 1000.0
    double gradientTolerance = gradientOptions.tolerance
    double width = 2.0 * step
    double avLen = 0.0
    double avGrad = 0.0
    double expGrad2 = expGrad * expGrad

    int nTested = 0
    for (int k : atomsToTest) {
      Atom a0 = atoms[k]
      if (!a0.isActive()) {
        continue
      }

      nTested++
      double[] analytic = allAnalytic[k]

      // Coordinate array only includes active atoms.
      int i = allToActive.get(k)
      int i3 = i * 3
      int i0 = i3 + 0
      int i1 = i3 + 1
      int i2 = i3 + 2
      double[] numeric = new double[3]

      // Find numeric dX
      double orig = x[i0]
      x[i0] = x[i0] + step
      double e = energy.energy(x)
      x[i0] = orig - step
      e -= energy.energy(x)
      x[i0] = orig
      numeric[0] = e / width

      // Find numeric dY
      orig = x[i1]
      x[i1] = x[i1] + step
      e = energy.energy(x)
      x[i1] = orig - step
      e -= energy.energy(x)
      x[i1] = orig
      numeric[1] = e / width

      // Find numeric dZ
      orig = x[i2]
      x[i2] = x[i2] + step
      e = energy.energy(x)
      x[i2] = orig - step
      e -= energy.energy(x)
      x[i2] = orig
      numeric[2] = e / width

      double dx = analytic[0] - numeric[0]
      double dy = analytic[1] - numeric[1]
      double dz = analytic[2] - numeric[2]
      double len = dx * dx + dy * dy + dz * dz
      avLen += len
      len = sqrt(len)

      double grad2 =
          analytic[0] * analytic[0] + analytic[1] * analytic[1] + analytic[2] * analytic[2]
      avGrad += grad2

      if (len > gradientTolerance) {
        logger.info(format(" %s\n Failed: %10.6f\n", a0.toString(), len)
            + format(" Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
            + format(" Numeric:  (%12.4f, %12.4f, %12.4f)\n", numeric[0], numeric[1], numeric[2]))
        ++nFailures
      } else {
        logger.info(format(" %s\n Passed: %10.6f\n", a0.toString(), len)
            + format(" Analytic: (%12.4f, %12.4f, %12.4f)\n", analytic[0], analytic[1], analytic[2])
            + format(" Numeric:  (%12.4f, %12.4f, %12.4f)", numeric[0], numeric[1], numeric[2]))
      }

      if (grad2 > expGrad2) {
        logger.info(format(" Atom %d has an unusually large gradient: %10.6f", i + 1, Math.sqrt(
            grad2)))
      }
      logger.info("\n")
    }

    avLen = avLen / nTested
    avLen = sqrt(avLen)
    if (avLen > gradientTolerance) {
      logger.info(format(" Test failure: RMSD from analytic solution is %10.6f > %10.6f", avLen,
          gradientTolerance))
    } else {
      logger.info(format(" Test success: RMSD from analytic solution is %10.6f < %10.6f", avLen,
          gradientTolerance))
    }
    logger.info(format(" Number of atoms failing analytic test: %d", nFailures))

    avGrad = avGrad / nTested
    avGrad = sqrt(avGrad)
    if (avGrad > expGrad) {
      logger.info(format(" Unusually large RMS gradient: %10.6f > %10.6f", avGrad, expGrad))
    } else {
      logger.info(format(" RMS gradient: %10.6f", avGrad))
    }

    // Set all ESV variables to 0.5
    int numESVs = titratingESVs.size()
    for (int i=0; i<numESVs; i++) {
      esvSystem.setLambda(i, 0.5)
    }

    energy.getCoordinates(x)
    double e = energy.energyAndGradient(x, g)
    double[] esvDerivs = esvSystem.derivatives

    // Check the dU/dL_i analytic results vs. finite-differences for extended system variables.
    // Loop over extended system variables
    for (int i=0; i<numESVs; i++) {
      esvSystem.setLambda(i, 0.5 + step)
      double ePlus = energy.energy(x)
      esvSystem.setLambda(i, 0.5 - step)
      double eMinus = energy.energy(x)
      esvSystem.setLambda(i, 0.5)

      double fdDeriv = (ePlus - eMinus) / width
      double error = abs(fdDeriv - esvDerivs[i])

      if (error > gradientTolerance) {
        logger.info(format(" ESV %d\n Failed: %10.6f\n", i, error)
            + format(" Analytic: %12.4f vs. Numeric: %12.4f\n", esvDerivs[i], fdDeriv))
        ++nESVFailures
      } else {
        logger.info(format(" ESV %d\n Passed: %10.6f\n", i, error)
            + format(" Analytic: %12.4f vs. Numeric: %12.4f\n", esvDerivs[i], fdDeriv))
      }
    }

    if(print){
      if(numESVs < 4){
        printPermutations(esvSystem, numESVs, energy, x)
      }
    }

    return this
  }

  private void printPermutations(ExtendedSystem esvSystem, int numESVs, ForceFieldEnergy energy, double[] x){
    for (int i=0; i<numESVs; i++) {
      esvSystem.setLambda(i, 0.0)
    }
    energy.getCoordinates(x)
    printPermutationsR(esvSystem, numESVs-1, energy, x)
  }

  private void printPermutationsR(ExtendedSystem esvSystem, int esvID, ForceFieldEnergy energy, double[] x){
    for (int i=0; i<=1; i++){
      esvSystem.setLambda(esvID,(double) i)
      if(esvID!=0){
        printPermutationsR(esvSystem, esvID-1, energy, x)
      }
      else{
        double[] energyAndInteractionList = new double [26]
        String lambdaList = esvSystem.getLambdaList()
        logger.info(format("Lambda List: %s", lambdaList))
        //Add ForceFieldEnergy to hashmap for testing. Protonation endstates used as key in map.
        energy.energy(x,true)

        // Bond Energy
        energyAndInteractionList[0] = energy.getBondEnergy()
        energyAndInteractionList[1] = (double) energy.getNumberofBonds()
        // Angle Energy
        energyAndInteractionList[2] = energy.getAngleEnergy()
        energyAndInteractionList[3] = (double) energy.getNumberofAngles()
        // Stretch-Bend Energy
        energyAndInteractionList[4] = energy.getStrenchBendEnergy()
        energyAndInteractionList[5] = (double) energy.getNumberofStretchBends()
        // Urey-Bradley Energy,
        energyAndInteractionList[6] = energy.getUreyBradleyEnergy()
        energyAndInteractionList[7] = (double) energy.getNumberofUreyBradleys()
        // Out-of-Plane Bend
        energyAndInteractionList[8] = energy.getOutOfPlaneBendEnergy()
        energyAndInteractionList[9] = (double) energy.getNumberofOutOfPlaneBends()
        // Torsional Angle
        energyAndInteractionList[10] = energy.getTorsionEnergy()
        energyAndInteractionList[11] = (double) energy.getNumberofTorsions()
        // Improper Torsional Angle
        energyAndInteractionList[12] = energy.getImproperTorsionEnergy()
        energyAndInteractionList[13] = (double) energy.getNumberofImproperTorsions()
        // Pi-Orbital Torsion
        energyAndInteractionList[14] = energy.getPiOrbitalTorsionEnergy()
        energyAndInteractionList[15] = (double) energy.getNumberofPiOrbitalTorsions()
        // Torsion-Torsion
        energyAndInteractionList[16] = energy.getTorsionTorsionEnergy()
        energyAndInteractionList[17] = (double) energy.getNumberofTorsionTorsions()
        // van Der Waals
        energyAndInteractionList[18] = energy.getVanDerWaalsEnergy()
        energyAndInteractionList[19] = (double) energy.getVanDerWaalsInteractions()
        // Permanent Multipoles
        energyAndInteractionList[20] = energy.getPermanentMultipoleEnergy()
        energyAndInteractionList[21] = (double) energy.getPermanentInteractions()
        // Polarization Energy
        energyAndInteractionList[22] = energy.getPolarizationEnergy()
        energyAndInteractionList[23] = (double) energy.getPermanentInteractions()
        // Extended System Bias
        energyAndInteractionList[24] = energy.getEsvBiasEnergy()
        // Total Energy
        energyAndInteractionList[25] = energy.getTotalEnergy()

        endstateEnergyMap.put(lambdaList,energyAndInteractionList)
        logger.info(format("\n"))
      }
    }
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (energy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(energy)
    }
    return potentials
  }

}
