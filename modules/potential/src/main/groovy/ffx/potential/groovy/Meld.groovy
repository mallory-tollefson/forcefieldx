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
import ffx.potential.bonded.Residue
import ffx.potential.bonded.ResidueEnumerations

import java.util.logging.Level
import static java.lang.String.format

import com.google.common.collect.MinMaxPriorityQueue
import com.sun.jna.ptr.PointerByReference

import org.apache.commons.io.FilenameUtils
import org.apache.commons.lang3.StringUtils

import ffx.numerics.Potential
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter

import edu.uiowa.jopenmm.MeldOpenMMLibrary

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Energy &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy.", name = "ffxc Energy")
class Energy extends PotentialScript {

    /**
     * -g or --gradient to print out gradients.
     */
    @Option(names = ['-g', '--gradient'], paramLabel = "false",
            description = 'Print out atomic gradients as well as energy.')
    private boolean gradient = false

    /**
     * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
     */
    @Option(names = ['--es1', '--noElecStart1'], paramLabel = "1",
            description = 'Starting no-electrostatics atom for 1st topology')
    private int es1 = 1
    /**
     * * --fl or --findLowest Return the n lowest energy structures from an ARC or PDB file.
     */

    @Option(names = ['--fl', '--findLowest'], paramLabel = "0",
            description = 'Return the n lowest energy structures from an ARC or PDB file.')
    private int fl = 0

    /**
     * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
     */
    @Option(names = ['--ef1', '--noElecFinal1'], paramLabel = "-1",
            description = 'Final no-electrostatics atom for 1st topology')
    private int ef1 = -1

    /**
     * --ss or --secondaryStructure defines the predicted presence of secondary structure.
     */
    @Option(names = ['--ss', '--secondaryStructure'], paramLabel = "",
            description = 'Secondary structure prediction. H = Helix, E = Beta Sheet, . = Coil')
    private String secondaryStructure = ""

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
    Energy run() {

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
        logger.info(" Running Energy on " + filename)

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
        Atom[] atoms = activeAssembly.getAtomArray()

        // Apply the no electrostatics atom selection
        int noElecStart = es1
        noElecStart = (noElecStart < 1) ? 1 : noElecStart

        int noElecStop = ef1
        noElecStop = (noElecStop > atoms.length) ? atoms.length : noElecStop

        if (noElecStart <= noElecStop) {
            logger.info(format(" Disabling electrostatics for atoms %d (%s) to %d (%s).",
                    noElecStart, atoms[noElecStart - 1], noElecStop, atoms[noElecStop - 1]))
        }
        for (int i = noElecStart; i <= noElecStop; i++) {
            Atom ai = atoms[i - 1]
            ai.setElectrostatics(false)
            ai.print(Level.FINE)
        }

        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)

        if (gradient) {
            double[] g = new double[nVars]
            int nAts = nVars / 3
            energy = forceFieldEnergy.energyAndGradient(x, g, true)
            logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"))
            for (int i = 0; i < nAts; i++) {
                int i3 = 3 * i
                logger.info(format(" %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]))
            }
        } else {
            energy = forceFieldEnergy.energy(x, true)
        }

        SystemFilter systemFilter = potentialFunctions.getFilter()

        if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {

            int numSnaps = fl
            double lowestEnergy = energy
            assemblyState = new AssemblyState(activeAssembly)
            int index = 1

            // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
            MinMaxPriorityQueue<StateContainer> lowestEnergyQueue = null
            if (fl > 0) {
                lowestEnergyQueue = MinMaxPriorityQueue
                        .maximumSize(numSnaps)
                        .create()
                lowestEnergyQueue.add(new StateContainer(assemblyState, lowestEnergy))
            }

            while (systemFilter.readNext()) {
                index++
                forceFieldEnergy.getCoordinates(x)
                if (verbose) {
                    logger.info(format(" Snapshot %4d", index));
                    energy = forceFieldEnergy.energy(x, true);
                } else {
                    energy = forceFieldEnergy.energy(x, false)
                    logger.info(format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy))
                }

                if (fl > 0) {
                    lowestEnergyQueue.add(new StateContainer(new AssemblyState(activeAssembly), energy))
                }
            }

            if (fl > 0) {
                if (numSnaps > index) {
                    logger.warning(format(
                            " Requested %d snapshots, but file %s has only %d snapshots. All %d energies will be reported",
                            numSnaps, filename, index, index))
                    numSnaps = index
                }

                File saveDir = baseDir
                String modelFilename = activeAssembly.getFile().getAbsolutePath()
                if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
                    saveDir = new File(FilenameUtils.getFullPath(modelFilename))
                }
                String dirName = saveDir.toString() + File.separator
                String fileName = FilenameUtils.getName(modelFilename)

                for (int i = 0; i < numSnaps - 1; i++) {
                    StateContainer savedState = lowestEnergyQueue.removeLast()
                    AssemblyState finalAssembly = savedState.getState()
                    finalAssembly.revertState()
                    double finalEnergy = savedState.getEnergy()
                    logger.info(format(" The potential energy found is %16.8f (kcal/mol)", finalEnergy))
                    File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
                    MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
                    potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
                }

                StateContainer savedState = lowestEnergyQueue.removeLast()
                AssemblyState lowestAssembly = savedState.getState()
                lowestEnergy = savedState.getEnergy()

                assemblyState.revertState()
                logger.info(format(" The lowest potential energy found is %16.8f (kcal/mol)", lowestEnergy))

                // Prints our final energy (which will be the lowest energy
                File saveFile = potentialFunctions.versionFile(new File(dirName + fileName))
                MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
                potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
            }
        }

        // Validate that secondary structure restraints and automatically modify as necessary.
        secondaryStructure = validateSecondaryStructurePrediction()
        checkForAppropriateResidueIdentities()

        //EXAMPLE MELD FORCE
        PointerByReference meldForce = MeldOpenMMLibrary.OpenMM_MeldForce_create()

        // Set up MELD restraints for secondary structure. Force constants and quadratic cut values were set
        // by the Dr. Ken Dill Research Group
        meldForce = setUpSecondaryMeldRestraints(meldForce, 2.48, 2.48, 2)

        return this

    }

    /**
     * This method validates that the user-supplied secondary structure predictions are the correct length and contain
     * the correct characters.
     * @return String containing the validated and edited secondary structure restraints.
     */
    String validateSecondaryStructurePrediction() {
        // The only characters that should be present in secondary structure restraint string are 'H' for helix, 'E'
        // for beta sheet and '.' for coil.
        if (!secondaryStructure.matches("^[HE.]+") && !secondaryStructure.isEmpty()) {
            logger.severe(" Secondary structure restraints may only contain characters 'H', 'E' and '.'")
        }

        MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
        int numResidues = molecularAssembly.getResidueList().size()
        int numSecondaryStructure = secondaryStructure.length()

        // Only one secondary structure restraint should exist per residue.
        if (numSecondaryStructure == 0) {
            logger.warning(" No secondary structure restraints have been provided. Simulation will proceed " +
                    "with all residues having random coil secondary structure restraints.")
            String randomCoil = StringUtils.leftPad("", numResidues, ".")
            return randomCoil
        } else if (numSecondaryStructure < numResidues) {
            logger.warning(" Too few secondary structure restraints exist for number of residues present. " +
                    "Random coil will be added to end residues without provided secondary structure restraints.")
            String extraCoil = StringUtils.rightPad(secondaryStructure, numResidues, '.')
            return extraCoil
        } else if (numSecondaryStructure == numResidues) {
            logger.info(" Secondary structure restraints will be added for all residues.")
            return secondaryStructure
        } else if (numSecondaryStructure > numResidues) {
            logger.warning(" Too many secondary structure restraints exist for number of residues present."
                    + " Provided secondary structure restraints will be truncated.")
            String truncated = secondaryStructure.substring(0, numResidues)
            return truncated
        } else {
            logger.severe(" Secondary structure restraints or residues do not exist.")
        }
    }

    /**
     * This method checks that secondary structure assignments are appropriate for the residue identity. ACE and NME
     * residues do not have alpha carbons, so they are not compatible with the alpha helix or beta sheet MELD
     * restraints.
     */
    void checkForAppropriateResidueIdentities(){
        MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
        ArrayList<Residue> residues = molecularAssembly.getResidueList()
        for(int i=0; i<secondaryStructure.length(); i++){
            Residue residue = residues.get(i)
            ResidueEnumerations.AminoAcid3 aminoAcid3 = residue.getAminoAcid3()

            String aminoAcidString = aminoAcid3.toString()
            String NMEString = Residue.AA3.NME.toString()
            String ACEString = Residue.AA3.ACE.toString()

            if(aminoAcidString.equals(NMEString) || aminoAcidString.equals((ACEString))){
                if(secondaryStructure[i].equals('H')){
                    logger.info(" Secondary structure was modified to accommodate non-standard amino acid residue.")
                    secondaryStructure = secondaryStructure.substring(0,i) + '.' + secondaryStructure.substring(i+1)
                } else if (secondaryStructure[i].equals('E')){
                    logger.info(" Secondary structure was modified to accommodate non-standard amino acid residue.")
                    secondaryStructure = secondaryStructure.substring(0,i) + '.' + secondaryStructure.substring(i+1)
                }
            }
        }
    }

    /**
     * This method sets up MELD torsion and distance restraints for secondary structure elements.
     * @param torsionForceConstant A float with value supplied by The Dill Group. In kJ/mol/(10 degree)^2
     * @param distanceForceConstant A float with value supplied by The Dill Group. In kJ/mol/Angstrom^2.
     * @param quadraticCut A float with value supplied by The Dill Group. This tells where to begin having a quadratic
     * bottom rather than a flat bottom on the force. In Angstroms.
     * @return The meldForce PointerByReference.
     */
    PointerByReference setUpSecondaryMeldRestraints(PointerByReference meldForce, float torsionForceConstant, float distanceForceConstant, float quadraticCut) {
        torsionForceConstant /= 100
        distanceForceConstant *= 100 // Convert to kJ/mol/nm^2
        quadraticCut *= 10 //Convert to nm.
        int minNumResForSecondary = 5

        String helixChar = 'H'
        String sheetChar = 'E'
        String coilChar = '.'

        ArrayList<ArrayList<Integer>> helices = extractSecondaryElement(secondaryStructure, helixChar, minNumResForSecondary)
        ArrayList<ArrayList<Integer>> sheets = extractSecondaryElement(secondaryStructure, sheetChar, minNumResForSecondary)

        float phi = -62.5
        float deltaPhi = 17.5
        float psi = -42.5
        float deltaPsi = 17.5
        meldForce = setRestraintsByElementType(meldForce, helices, phi, deltaPhi, psi, deltaPsi, quadraticCut, torsionForceConstant, distanceForceConstant)

        phi = -117.5
        deltaPhi = 27.5
        psi = 145.0
        deltaPsi = 25.0
        meldForce = setRestraintsByElementType(meldForce, sheets, phi, deltaPhi, psi, deltaPsi, quadraticCut, torsionForceConstant, distanceForceConstant)

        return meldForce
    }

    /**
     * This method sets the torsion and distance restraints for a secondary element type. For example, if the starting
     * and ending points for all helices are provided, then the torsion and distance restraints appropriate for helices
     * will be set.
     * @param meldForce A PointerByReference containing the meld forces.
     * @param elements An ArrayList<ArrayList<Integer>> that contains the starting and ending points for secondary
     * structure elements (helices, sheets, coils, etc).
     * @param phi A float with value supplied by The Dill Group.
     * @param deltaPhi A float with value supplied by The Dill Group.
     * @param psi A float with value supplied by The Dill Group.
     * @param deltaPsi A float with value supplied by The Dill Group.
     * @param quadraticCut A float with value supplied by The Dill Group.
     * @param torsionForceConstant A float with value supplied by The Dill Group.
     * @param distanceForceConstant A float with value supplied by The Dill Group.
     * @return The meldForce PointerByReference.
     */
    PointerByReference setRestraintsByElementType(PointerByReference meldForce, ArrayList<ArrayList<Integer>> elements, float phi, float deltaPhi, float psi, float deltaPsi, float quadraticCut, float torsionForceConstant, float distanceForceConstant){
        if(!elements.isEmpty()){
            for (int i = 0; i <= elements.size(); i++) { // For every secondary element.
                ArrayList<Integer> element = elements.pop()
                int elementStart = (int) element.pop()
                int elementEnd = (int) element.pop()
                ArrayList<Residue> residues = activeAssembly.getResidueList()
                for (int j = elementStart + 1; j < elementEnd - 1; j++) { // For every residue in an element.
                    Residue residueJminus1 = residues.get(j - 1)
                    Residue residueJ = residues.get(j)
                    Residue residueJplus1 = residues.get(j + 1)

                    Atom carbon = (Atom) residueJ.getAtomNode("C")
                    Atom alphaC = (Atom) residueJ.getAtomNode("CA")
                    Atom nitrogen = (Atom) residueJ.getAtomNode("N")
                    Atom lastCarbon = (Atom) residueJminus1.getAtomNode("C")
                    Atom nextNitrogen = (Atom) residueJplus1.getAtomNode("N")

                    //Phi torsion restraint.
                    MeldOpenMMLibrary.OpenMM_MeldForce_addTorsionRestraint(meldForce, lastCarbon.getArrayIndex(), nitrogen.getArrayIndex(), alphaC.getArrayIndex(), carbon.getArrayIndex(), phi, deltaPhi, torsionForceConstant)

                    //Psi torsion restraint.
                    MeldOpenMMLibrary.OpenMM_MeldForce_addTorsionRestraint(meldForce, nitrogen.getArrayIndex(), alphaC.getArrayIndex(), carbon.getArrayIndex(), nextNitrogen.getArrayIndex(), psi, deltaPsi, torsionForceConstant)
                }

                // DISTANCE RESTRAINTS FOR HELICES
                Residue elementStartRes = residues.get(elementStart)
                Residue elementStartResPlus1 = residues.get(elementStart+1)
                Residue elementStartResPlus3 = residues.get(elementStart+3)
                Residue elementStartResPlus4 =  residues.get(elementStart+4)

                Atom alphaC = (Atom) elementStartRes.getAtomNode("CA")
                Atom alphaCPlus1 = (Atom) elementStartResPlus1.getAtomNode("CA")
                Atom alphaCPlus3 = (Atom) elementStartResPlus3.getAtomNode("CA")
                Atom alphaCPlus4 = (Atom) elementStartResPlus4.getAtomNode("CA")

                //The four floats (r1-r4) are in nanometers and were provided by the Dill Research Group
                try {
                    int alphaCIndex = alphaC.getArrayIndex()
                    int alphaCPlus3Index = alphaCPlus3.getArrayIndex()
                    float r1 = 0 // Units are nanometers for all r constants.
                    float r2 = 0.785
                    float r3 = 1.063
                    float r4 = 1.063 + quadraticCut
                    MeldOpenMMLibrary.OpenMM_MeldForce_addDistanceRestraint(meldForce, alphaCIndex, alphaCPlus3Index, r1, r2, r3, r4, distanceForceConstant)
                }catch(Exception e){
                    logger.severe(" Meld distance restraint cannot be added.\n"+e.printStackTrace())
                }

                try {
                    int alphaCPlus1Index = alphaCPlus1.getArrayIndex()
                    int alphaCPlus4Index = alphaCPlus4.getArrayIndex()
                    float r1 = 0 // Units are nanometers for all r constants.
                    float r2 = 0.785
                    float r3 = 1.063
                    float r4 = 1.063 + quadraticCut
                    MeldOpenMMLibrary.OpenMM_MeldForce_addDistanceRestraint(meldForce, alphaCPlus1Index, alphaCPlus4Index, r1, r2, r3, r4, distanceForceConstant)
                }catch(Exception e){
                    logger.severe(" Meld distance restraint cannot be added.\n"+e.printStackTrace())
                }

                try{
                    int alphaCIndex = alphaC.getArrayIndex()
                    int alphaCPlus4Index = alphaCPlus4.getArrayIndex()
                    float r1 = 0 // Units are nanometers for all r constants.
                    float r2 = 1.086
                    float r3 = 1.394
                    float r4 = 1.394 + quadraticCut
                    MeldOpenMMLibrary.OpenMM_MeldForce_addDistanceRestraint(meldForce, alphaCIndex, alphaCPlus4Index, r1, r2, r3, r4, distanceForceConstant)
                }catch(Exception e){
                    logger.severe(" Meld distance restraint cannot be added.\n"+e.printStackTrace())
                }
            }
        }
        return meldForce
    }

    /**
     * This method determines the starting and ending indices for secondary elements of the requested type based
     * on the user-supplied secondary structure restraint predictions. The Dill Group requires that secondary
     * elements have at least three consecutive residues to be considered a secondary element.
     * @param ss A string of the secondary structure prediction.
     * @param elementType Character indicating type of secondary element being searched (helix, coil, sheet).
     * @param minNumResidues Integer minimum of consecutive secondary structure predictions
     * to create a secondary element.
     * @return ArrayList<ArrayList<Integer>     > Contains starting and ending residues for each secondary element.
     */
    ArrayList<ArrayList<Integer>> extractSecondaryElement(String ss, String elementType, int minNumResidues) {
        ArrayList<ArrayList<Integer>> allElements = new ArrayList<ArrayList<Integer>>()
        //Will hold starting and ending indices for all found secondary elements of the requested type.
        int lastMatch = 0 //Track of the most recent index to have a character matching the requested elementType.
        int i = 0 //Iterates through each index in the secondary structure string.
        while (i < ss.length()) {
            if (ss[i].equals(elementType)) {
                int elementStartIndex = i
                //Set the starting index for the secondary element as soon as the value at the ith index matches the
                //requested element type.
                for (int j = i + 1; j <= ss.length(); j++) {
                    //Use the jth index to iterate through the secondary structure prediction until the end of the
                    // secondary element is found.
                    if (j < ss.length()) {
                        if (!ss[j].equals(elementType)) {
                            if (j == lastMatch + 1) {
                                //If the most recent lastMatch is only one index away, then check and store the
                                // starting and ending indices of the secondary element.
                                i = j
                                //Set i=j so that i begins searching for the next element at the end of the most recent
                                // secondary element.
                                int elementLength = j - elementStartIndex
                                if (elementLength > minNumResidues) {
                                    //If secondary element is above minimum length, store starting and ending indices
                                    // of secondary element.
                                    ArrayList<Integer> currentElement = new ArrayList<Integer>()
                                    currentElement.add((Integer) elementStartIndex)
                                    currentElement.add((Integer) lastMatch)
                                    allElements.add(currentElement)
                                }
                                j = ss.length() + 1
                                //Since end of current secondary element has been found, exit inner j loop.
                            }
                        } else {
                            lastMatch = j
                            i++
                            //If the jth index in the secondary structure string matches the requested element,
                            // increment j until the end of the secondary element is found.
                        }
                    }
                    if (j == ss.length()) {
                        //Handle the case when a secondary element is at the very end of the secondary structure string.
                        i = ss.length() + 1
                        if (j == lastMatch + 1) {
                            int elementLength = j - elementStartIndex
                            if (elementLength > minNumResidues) {
                                ArrayList<Integer> currentElement = new ArrayList<Integer>()
                                currentElement.add((Integer) elementStartIndex)
                                currentElement.add((Integer) lastMatch)
                                allElements.add(currentElement)
                            }
                            j = ss.length() + 1
                            //Since end of current secondary element has been found, exit inner j loop.
                        }
                    }
                }
            } else {
                i++
                //Increment i until end of secondary structure prediction or until requested secondary element is found.
            }
        }
        return allElements
    }

    @Override
    List<Potential> getPotentials() {
        return forceFieldEnergy == null ? Collections.emptyList() : Collections.singletonList(forceFieldEnergy)
    }
}

