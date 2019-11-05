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
@Command(description = " Compute the force field potential energy.", name = "ffxc Meld")
class Meld extends PotentialScript {

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
    Meld run() {

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

        //EXAMPLE MELD FORCE
        PointerByReference meldForce = MeldOpenMMLibrary.OpenMM_MeldForce_create()

        // ArrayList to hold the collections of restraints (currently hydrophobic and secondary structure restraints)
        ArrayList<SelectivelyActiveCollection> collections = new ArrayList<SelectivelyActiveCollection>()

        // Validate that secondary structure restraints and automatically modify as necessary.
        secondaryStructure = validateSecondaryStructurePrediction()
        checkForAppropriateResidueIdentities()

        // Set up MELD restraints for secondary structure. Force constants and quadratic cut values were set
        // by the Dr. Ken Dill Research Group
        SelectivelyActiveCollection collectionSecondary = new SelectivelyActiveCollection()
        LinearRamp ramp = new LinearRamp(0.0,100.0,0.0,1.0)

        ConstantScaler constantScaler = new ConstantScaler()
        meldForce = setUpSecondaryMeldRestraints(meldForce, constantScaler, ramp, 2.48, 2.48, 2, collectionSecondary)
        collections.add(collectionSecondary)

        // Set up hydrophobic contact restraints.
        SelectivelyActiveCollection collectionHydrophobic = new SelectivelyActiveCollection()
        // Factor of 4.0 is the factor used by the Dill Group.
        NonLinearScaler nonLinearScaler = new NonLinearScaler(4.0)
        float contactsPerHydrophobe = 1.3
        meldForce = setUpHydrophobicRestraints(meldForce, nonLinearScaler, ramp, collectionHydrophobic, contactsPerHydrophobe)
        collections.add(collectionHydrophobic)

        //Likely create array of restraintEnergies for the meldForce here.
        ArrayList<Restraint> alwaysActiveRestraints = new ArrayList<Restraint>()
        MeldRestraintTransformer transformer = new MeldRestraintTransformer(meldForce, collections, alwaysActiveRestraints)
        transformer.addInteractions()

        // A forceGroup of 0 is for bonded forces; meld forces are bond-like.
        int forceGroup = 0;
        OpenMM_Force_setForceGroup(meldForce, forceGroup)

        logger.info("Adding Meld Force to ForceFieldEnergyOpenMM")
        ForceFieldEnergyOpenMM forceFieldEnergyOpenMM = (ForceFieldEnergyOpenMM) forceFieldEnergy
        forceFieldEnergyOpenMM.addForce(meldForce)
        forceFieldEnergyOpenMM.reinitContext()

        energy = forceFieldEnergy.energy(x, true)

        return this
    }

    /**
     * This method sets distance restraints between atoms of hydrophobic residues.
     * @param meldForce The meld force PointerByReference.
     * @param collectionHydrophobic A SelectivelyActiveCollection of hydrophobic RestraintGroups.
     * @param contactsPerHydrophobe A float representing the number of hydrophobic contacts that occur per hydrophobic
     * residue. The Dill Group uses 1.3.
     * @return The meld force PointerByReference
     */
    PointerByReference setUpHydrophobicRestraints(PointerByReference meldForce, NonLinearScaler nonLinearScaler, LinearRamp ramp, SelectivelyActiveCollection collectionHydrophobic, float contactsPerHydrophobe) {
        MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
        ArrayList<Residue> residueList = molecularAssembly.getResidueList()
        ArrayList<Residue> hydrophobicResidues1 = new ArrayList<Residue>()
        for (int i = 0; i < residueList.size(); i++) {
            Residue residue = residueList.get(i)
            ResidueEnumerations.AminoAcid3 aminoAcid3 = residue.getAminoAcid3()
            String aminoAcidString = aminoAcid3.toString()

            // Hydrophobic residues identities.
            String ALAString = Residue.AA3.ALA.toString()
            String VALString = Residue.AA3.VAL.toString()
            String LEUString = Residue.AA3.LEU.toString()
            String ILEString = Residue.AA3.ILE.toString()
            String PHEString = Residue.AA3.PHE.toString()
            String TRPString = Residue.AA3.TRP.toString()
            String METString = Residue.AA3.MET.toString()
            String PROString = Residue.AA3.PRO.toString()

            if (aminoAcidString.equals(ALAString) || aminoAcidString.equals(VALString) || aminoAcidString.equals(LEUString)
                    || aminoAcidString.equals(ILEString) || aminoAcidString.equals(PHEString) || aminoAcidString.equals(TRPString)
                    || aminoAcidString.equals(METString) || aminoAcidString.equals(PROString)) {
                hydrophobicResidues1.add(residue)
            }
        }

        // If no hydrophobic residues are present, then return from the method. Else, set up hydrophobic restraints.
        if (hydrophobicResidues1.isEmpty()) {
            logger.warning(" No hydrophobic residues in sequence. No hydrophobic restraints will be added to the system.")
            return meldForce
        } else {
            // The pairs ArrayList keeps track of residue indices that should get hydrophobic restraints.
            ArrayList<ArrayList<Integer>> pairs = new ArrayList<ArrayList<Integer>>()
            for (int i = 0; i < hydrophobicResidues1.size(); i++) {
                for (int j = 0; j < hydrophobicResidues1.size(); j++) {
                    if (i == j) {
                        continue
                    }
                    int residueNumber1 = hydrophobicResidues1.get(i).getResidueNumber()
                    int residueNumber2 = hydrophobicResidues1.get(j).getResidueNumber()
                    if (Math.abs(residueNumber1 - residueNumber2) < 7) {
                        continue
                    }

                    // Check that a pair is not being added twice.
                    boolean alreadyStored = false
                    for (ArrayList<Integer> pair : pairs) {
                        Integer index1 = pair.get(0)
                        Integer index2 = pair.get(1)
                        if ((index1.intValue() == i && index2.intValue() == j) || (index1.intValue() == j && index2.intValue() == i)) {
                            alreadyStored = true
                        }
                    }
                    if(alreadyStored == false) {
                        ArrayList<Integer> newPair = new ArrayList<Integer>()
                        newPair.add((Integer) i)
                        newPair.add((Integer) j)
                        pairs.add(newPair)
                    }
                }
            }

            // restraintGroups will hold each RestraintGroup and will be added to the collection at the end
            ArrayList<RestraintGroup> restraintGroups = new ArrayList<RestraintGroup>()
            // Set up hydrophobic restraints for each pair.
            for(ArrayList<Integer> pair: pairs){
                int index1 = (int) pair.get(0)
                int index2 = (int) pair.get(1)
                Residue residue1 = hydrophobicResidues1.get(index1)
                Residue residue2 = hydrophobicResidues1.get(index2)
                ArrayList<Atom> atoms1 = residue1.getAtomList()
                ArrayList<Atom> atoms2 = residue2.getAtomList()

                // Names of atoms that can have hydrophobic restraints.
                ArrayList<String> atomNamesForRestraints = new ArrayList<String>()
                atomNamesForRestraints.add("CA")
                atomNamesForRestraints.add("CB")
                atomNamesForRestraints.add("CD")
                atomNamesForRestraints.add("CD1")
                atomNamesForRestraints.add("CD2")
                atomNamesForRestraints.add("CE")
                atomNamesForRestraints.add("CE1")
                atomNamesForRestraints.add("CE2")
                atomNamesForRestraints.add("CE3")
                atomNamesForRestraints.add("CG")
                atomNamesForRestraints.add("CG1")
                atomNamesForRestraints.add("CG2")
                atomNamesForRestraints.add("CG3")
                atomNamesForRestraints.add("CH2")
                atomNamesForRestraints.add("CZ")
                atomNamesForRestraints.add("CZ2")
                atomNamesForRestraints.add("CZ3")
                atomNamesForRestraints.add("NE1")
                atomNamesForRestraints.add("SD")

                // hydrophobicRestraints holds all Restraints that are created for a particular pair of residues.
                ArrayList<Restraint> hydrophobicRestraints = new ArrayList<>()
                for(Atom atom1: atoms1) {
                    String name1 = atom1.getName()
                    for (Atom atom2 : atoms2) {
                        String name2 = atom2.getName()
                        if (atomNamesForRestraints.contains(name1) && atomNamesForRestraints.contains(name2)) {
                            int atom1Index = atom1.getArrayIndex()
                            int atom2Index = atom2.getArrayIndex()
                            //Create distance restraints between specified atoms of both residues.
                            DistanceRestraint distanceRestraint = new DistanceRestraint(meldForce, nonLinearScaler, ramp, atom1Index, atom2Index, 0.0, 0.0, 0.5, 0.7, 250.0)
                            hydrophobicRestraints.add(distanceRestraint)
                        }
                    }
                }
                if (!hydrophobicRestraints.isEmpty()) {
                    RestraintGroup restraintGroup = new RestraintGroup(hydrophobicRestraints, 1)
                    restraintGroups.add(restraintGroup)
                }
            }
            // Add restraintGroups to the collection and set the numActive for the collection.
            collectionHydrophobic.setRestraintGroups(restraintGroups)
            int active = (int) contactsPerHydrophobe * hydrophobicResidues1.size()
            collectionHydrophobic.setNumActive(active)
        }
        return meldForce
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
            String randomCoil = org.apache.commons.lang3.StringUtils.leftPad("", numResidues, ".")
            return randomCoil
        } else if (numSecondaryStructure < numResidues) {
            logger.warning(" Too few secondary structure restraints exist for number of residues present. " +
                    "Random coil will be added to end residues without provided secondary structure restraints.")
            String extraCoil = org.apache.commons.lang3.StringUtils.rightPad(secondaryStructure, numResidues, '.')
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
    void checkForAppropriateResidueIdentities() {
        MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
        ArrayList<Residue> residues = molecularAssembly.getResidueList()
        for (int i = 0; i < secondaryStructure.length(); i++) {
            Residue residue = residues.get(i)
            ResidueEnumerations.AminoAcid3 aminoAcid3 = residue.getAminoAcid3()

            String aminoAcidString = aminoAcid3.toString()
            String NMEString = Residue.AA3.NME.toString()
            String ACEString = Residue.AA3.ACE.toString()

            if (aminoAcidString.equals(NMEString) || aminoAcidString.equals((ACEString))) {
                if (secondaryStructure[i].equals('H')) {
                    logger.info(" Secondary structure was modified to accommodate non-standard amino acid residue.")
                    secondaryStructure = secondaryStructure.substring(0, i) + '.' + secondaryStructure.substring(i + 1)
                } else if (secondaryStructure[i].equals('E')) {
                    logger.info(" Secondary structure was modified to accommodate non-standard amino acid residue.")
                    secondaryStructure = secondaryStructure.substring(0, i) + '.' + secondaryStructure.substring(i + 1)
                }
            }
        }
    }

    /**
     * This method sets up MELD torsion and distance restraints for secondary structure elements and adds them to groups
     * and collections.
     * @param torsionForceConstant A float with value supplied by The Dill Group. In kJ/mol/(10 degree)^2
     * @param distanceForceConstant A float with value supplied by The Dill Group. In kJ/mol/Angstrom^2.
     * @param quadraticCut A float with value supplied by The Dill Group. This tells where to begin having a quadratic
     * bottom rather than a flat bottom on the force. In Angstroms.
     * @return The meldForce PointerByReference.
     */
    PointerByReference setUpSecondaryMeldRestraints(PointerByReference meldForce, ConstantScaler constantScaler, LinearRamp ramp, float torsionForceConstant, float distanceForceConstant, float quadraticCut, SelectivelyActiveCollection collection) {
        torsionForceConstant /= 100
        distanceForceConstant *= 100 // Convert to kJ/mol/nm^2
        quadraticCut *= 10 //Convert to nm.
        int minNumResForSecondary = 4
        int runLength = 5

        String helixChar = 'H'
        String sheetChar = 'E'

        ArrayList<ArrayList<Integer>> helices = extractSecondaryElement(secondaryStructure, helixChar, minNumResForSecondary, runLength)
        ArrayList<ArrayList<Integer>> sheets = extractSecondaryElement(secondaryStructure, sheetChar, minNumResForSecondary, runLength)

        float phi = -62.5
        float deltaPhi = 17.5
        float psi = -42.5
        float deltaPsi = 17.5
        ArrayList<RestraintGroup> helixRestraintGroupList = new ArrayList<RestraintGroup>()
        meldForce = setRestraintsByElementType(meldForce, constantScaler, ramp, helixRestraintGroupList, helixChar, helices, phi, deltaPhi, psi, deltaPsi, quadraticCut, torsionForceConstant, distanceForceConstant)

        phi = -117.5
        deltaPhi = 27.5
        psi = 145.0
        deltaPsi = 25.0
        ArrayList<RestraintGroup> sheetRestraintGroupList = new ArrayList<RestraintGroup>()
        meldForce = setRestraintsByElementType(meldForce, constantScaler, ramp, sheetRestraintGroupList, sheetChar, sheets, phi, deltaPhi, psi, deltaPsi, quadraticCut, torsionForceConstant, distanceForceConstant)

        // Set up the collection.
        ArrayList<RestraintGroup> allRestraintsGroupList = new ArrayList<RestraintGroup>()
        allRestraintsGroupList.addAll(helixRestraintGroupList)
        allRestraintsGroupList.addAll(sheetRestraintGroupList)
        int keep = (int) allRestraintsGroupList.size() * 0.7
        collection.setNumActive(keep)
        collection.setRestraintGroups(allRestraintsGroupList)

        return meldForce
    }

    /**
     * This method determines the starting and ending indices for secondary elements of the requested type based
     * on the user-supplied secondary structure restraint predictions. The Dill Group requires that secondary
     * elements are 5 residues long (runLength) and at least 4 residues must match the requested secondary element
     * type (minNumResidues) in order to be a secondary element.
     * @param ss A string of the secondary structure prediction.
     * @param elementType Character indicating type of secondary element being searched (helix, coil, sheet).
     * @param minNumResidues Integer minimum of matching secondary structure predictions within a run length
     * to create a secondary element. The Dill Group uses 3.
     * @param runLength The minimum number of residues necessary to create a secondary element. The Dill Group uses 5.
     * @return
     */
    ArrayList<ArrayList<Integer>> extractSecondaryElement(String ss, String elementType, int minNumResidues, int runLength) {
        //Will hold starting and ending indices for all found secondary elements of the requested type.
        ArrayList<ArrayList<Integer>> allElements = new ArrayList<ArrayList<Integer>>()
        //Iterate through the secondary structure prediction string in increments of 5 residues.
        for (int i = 0; i <= ss.length() - runLength; i++) {
            String substring = ss.substring(i, i + runLength)
            //MatchSum holds the number of times that the 5 residue segment matches the requested element type.
            int matchSum = 0;
            for (int j = 0; j < runLength; j++) {
                if (substring[j].equals(elementType)) {
                    matchSum++
                }
            }
            //If enough residues in the 5 residue segment match the requested element type, then store the starting
            // and ending indices for the secondary element.
            if (matchSum >= minNumResidues) {
                ArrayList<Integer> currentElement = new ArrayList<Integer>()
                currentElement.add((Integer) i)
                currentElement.add((Integer) i + runLength)
                allElements.add(currentElement)
            }
        }
        return allElements
    }

    /**
     * This method sets the torsion and distance restraints for a secondary element type. For example, if the starting
     * and ending points for all helices are provided, then the torsion and distance restraints appropriate for helices
     * will be set.
     * @param meldForce A PointerByReference containing the meld forces.
     * @param restraintGroupList An empty ArrayList of RestraintGroup objects that gets filled during the method.
     * @param elementType A character indicating the type of secondary structure being searched for (H for helices,
     * E for sheets)
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
    PointerByReference setRestraintsByElementType(PointerByReference meldForce, ConstantScaler constantScaler, LinearRamp ramp, ArrayList<RestraintGroup> restraintGroupList, String elementType, ArrayList<ArrayList<Integer>> elements, float phi, float deltaPhi, float psi, float deltaPsi, float quadraticCut, float torsionForceConstant, float distanceForceConstant) {
        if (!elements.isEmpty()) {
            for (int i = 0; i < elements.size(); i++) { // For every secondary element.
                ArrayList<Restraint> restraintList = new ArrayList<Restraint>()
                ArrayList<Integer> element = elements.get(i)
                int elementStart = (int) element.get(0)
                int elementEnd = (int) element.get(1)
                ArrayList<Residue> residues = activeAssembly.getResidueList()

                // TORSION RESTRAINTS
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
                    try {
                        int lastCarbonIndex = lastCarbon.getArrayIndex()
                        int nitrogenIndex = nitrogen.getArrayIndex()
                        int alphaCIndex = alphaC.getArrayIndex()
                        int carbonIndex = carbon.getArrayIndex()
                        TorsionRestraint torsionRestraint = new TorsionRestraint(meldForce, constantScaler, ramp, lastCarbonIndex, nitrogenIndex, alphaCIndex, carbonIndex, phi, deltaPhi, torsionForceConstant)
                        restraintList.add(torsionRestraint)
                    } catch (Exception e) {
                        logger.severe(" Meld phi torsion restraint cannot be created.\n" + e.printStackTrace())
                    }

                    //Psi torsion restraint.
                    try {
                        int nitrogenIndex = nitrogen.getArrayIndex()
                        int alphaCIndex = alphaC.getArrayIndex()
                        int carbonIndex = carbon.getArrayIndex()
                        int nextNitrogenIndex = nextNitrogen.getArrayIndex()
                        TorsionRestraint torsionRestraint = new TorsionRestraint(meldForce, constantScaler, ramp, nitrogenIndex, alphaCIndex, carbonIndex, nextNitrogenIndex, psi, deltaPsi, torsionForceConstant)
                        restraintList.add(torsionRestraint)
                    } catch (Exception e) {
                        logger.severe(" Meld psi torsion restraint cannot be created.\n" + e.printStackTrace())
                    }
                }

                // DISTANCE RESTRAINTS
                Residue elementStartRes = residues.get(elementStart)
                Residue elementStartResPlus1 = residues.get(elementStart + 1)
                Residue elementStartResPlus3 = residues.get(elementStart + 3)
                Residue elementStartResPlus4 = residues.get(elementStart + 4)

                Atom alphaC = (Atom) elementStartRes.getAtomNode("CA")
                Atom alphaCPlus1 = (Atom) elementStartResPlus1.getAtomNode("CA")
                Atom alphaCPlus3 = (Atom) elementStartResPlus3.getAtomNode("CA")
                Atom alphaCPlus4 = (Atom) elementStartResPlus4.getAtomNode("CA")

                //The four floats (r1-r4) are in nanometers and were provided by the Dill Research Group
                try {
                    int alphaCIndex = alphaC.getArrayIndex()
                    int alphaCPlus3Index = alphaCPlus3.getArrayIndex()
                    float r1
                    float r2
                    float r3
                    float r4
                    if (elementType.equals('H')) {
                        r1 = 0 // Units are nanometers for all r constants.
                        r2 = 0.485
                        r3 = 0.561
                        r4 = 0.561 + quadraticCut
                    } else if (elementType.equals('E')) {
                        r1 = 0 // Units are nanometers for all r constants.
                        r2 = 0.785
                        r3 = 1.063
                        r4 = 1.063 + quadraticCut
                    }
                    DistanceRestraint distanceRestraint = new DistanceRestraint(meldForce, constantScaler, ramp, alphaCIndex, alphaCPlus3Index, r1, r2, r3, r4, distanceForceConstant)
                    restraintList.add(distanceRestraint)
                } catch (Exception e) {
                    logger.severe(" Meld distance restraint cannot be created.\n" + e.printStackTrace())
                }

                try {
                    int alphaCPlus1Index = alphaCPlus1.getArrayIndex()
                    int alphaCPlus4Index = alphaCPlus4.getArrayIndex()
                    float r1
                    float r2
                    float r3
                    float r4
                    if (elementType.equals('H')) {
                        r1 = 0 // Units are nanometers for all r constants.
                        r2 = 0.485
                        r3 = 0.561
                        r4 = 0.561 + quadraticCut
                    } else if (elementType.equals('E')) {
                        r1 = 0 // Units are nanometers for all r constants.
                        r2 = 0.785
                        r3 = 1.063
                        r4 = 1.063 + quadraticCut
                    }
                    DistanceRestraint distanceRestraint = new DistanceRestraint(meldForce, constantScaler, ramp, alphaCPlus1Index, alphaCPlus4Index, r1, r2, r3, r4, distanceForceConstant)
                    restraintList.add(distanceRestraint)
                } catch (Exception e) {
                    logger.severe(" Meld distance restraint cannot be created.\n" + e.printStackTrace())
                }

                try {
                    int alphaCIndex = alphaC.getArrayIndex()
                    int alphaCPlus4Index = alphaCPlus4.getArrayIndex()
                    float r1
                    float r2
                    float r3
                    float r4
                    if (elementType.equals('H')) {
                        r1 = 0 // Units are nanometers for all r constants.
                        r2 = 0.581
                        r3 = 0.684
                        r4 = 0.684 + quadraticCut
                    } else if (elementType.equals('E')) {
                        r1 = 0 // Units are nanometers for all r constants.
                        r2 = 1.086
                        r3 = 1.394
                        r4 = 1.394 + quadraticCut
                    }
                    DistanceRestraint distanceRestraint = new DistanceRestraint(meldForce, constantScaler, ramp, alphaCIndex, alphaCPlus4Index, r1, r2, r3, r4, distanceForceConstant)
                    restraintList.add(distanceRestraint)
                } catch (Exception e) {
                    logger.severe(" Meld distance restraint cannot be created.\n" + e.printStackTrace())
                }
                RestraintGroup restraintGroup = new RestraintGroup(restraintList, restraintList.size())
                restraintGroupList.add(restraintGroup)
            }
        }
        return meldForce
    }

    /**
     * Restraint objects are added to the meldForce.
     */
    private class Restraint {
        PointerByReference meldForce
    }

    /**
     * SelectableRestraints can be turned on or off throughout the simulation.
     */
    private class SelectableRestraint extends Restraint {
    }

    /**
     * AlwaysOnRestraints cannot be turned on/off during the simulation.
     */
    private class AlwaysOnRestraint extends Restraint {
    }

    /**
     * DistanceRestraints restrain the distance beetween two atoms.
     */
    private class DistanceRestraint extends SelectableRestraint {
        RestraintScaler scaler
        LinearRamp ramp
        int alphaCIndex
        int alphaCPlus3Index
        float r1 // nanometers
        float r2 // nanometers
        float r3 // nanometers
        float r4 // nanometers
        float distanceForceConstant // kJ/mol/nm^2

        private DistanceRestraint(PointerByReference meldForce, RestraintScaler scaler, LinearRamp ramp, int alphaCIndex, int alphaCPlus3Index, float r1, float r2, float r3, float r4, float distanceForceConstant) {
            this.scaler = scaler
            this.ramp = ramp
            this.meldForce = meldForce
            this.alphaCIndex = alphaCIndex
            this.alphaCPlus3Index = alphaCPlus3Index
            this.r1 = r1
            this.r2 = r2
            this.r3 = r3
            this.r4 = r4
            this.distanceForceConstant = distanceForceConstant
        }
    }

    /**
     * TorsionRestraints restrain the angle created by selected atoms.
     */
    private class TorsionRestraint extends SelectableRestraint {
        RestraintScaler scaler
        LinearRamp ramp
        int atom1Index
        int atom2Index
        int atom3Index
        int atom4Index
        float angle // Degrees
        float deltaAngle // Degrees
        float torsionForceConstant // kJ/mol/Degree^2

        private TorsionRestraint(PointerByReference meldForce, RestraintScaler scaler, LinearRamp ramp, int atom1Index, int atom2Index, int atom3Index, int atom4Index, float angle, float deltaAngle, float torsionForceConstant) {
            this.scaler = scaler
            this.ramp = ramp
            this.meldForce = meldForce
            this.atom1Index = atom1Index
            this.atom2Index = atom2Index
            this.atom3Index = atom3Index
            this.atom4Index = atom4Index
            this.angle = angle
            this.deltaAngle = deltaAngle
            this.torsionForceConstant = torsionForceConstant
        }
    }

    /**
     * A RestraintGroup is a collection of Restraints that belong together (for example, all DistanceRestraints and
     * TorsionRestraints that are part of a single alpha helix).
     */
    private class RestraintGroup {
        ArrayList<Restraint> restraints
        int numActive
        int numRestraints

        private RestraintGroup(ArrayList<Restraint> restraints, int numActive) {
            if (restraints.isEmpty()) {
                logger.severe(" Restraint list cannot be empty.")
            }
            if (numActive < 0) {
                logger.severe(" NumActive cannot be less than 0.")
            }
            this.restraints = restraints
            this.numActive = numActive
            this.numRestraints = restraints.size()
            if (numActive > numRestraints) {
                logger.severe(" Number of active restraints must be <= number of total restraints.")
            }
        }

        private int getNumRestraints() {
            return numRestraints
        }

        private ArrayList<Restraint> getRestraints() {
            return restraints
        }
    }

    /**
     * A collection of multiple RestraintGroups.
     */
    private class SelectivelyActiveCollection {
        ArrayList<RestraintGroup> restraintGroups
        int numActive

        private SelectivelyActiveCollection() {
        }

        private SelectivelyActiveCollection(ArrayList<RestraintGroup> restraintGroups, int numActive) {
            this.restraintGroups = restraintGroups
            this.numActive = numActive

            if (restraintGroups.isEmpty()) {
                logger.severe("SelectivelyActiveCollection cannot have empty restraint list.")
            }
            int numGroups = restraintGroups.size()
            if (numActive > numGroups) {
                logger.severe("Number of active restraint groups must be less than number of groups.")
            }
        }

        private setRestraintGroups(ArrayList<RestraintGroup> restraintGroups) {
            this.restraintGroups = restraintGroups
        }

        private setNumActive(int numActive) {
            this.numActive = numActive
        }

        private addRestraint(RestraintGroup restraintGroup) {
            restraintGroups.add(restraintGroup)
        }

        private addRestraint(Restraint restraint) {
            ArrayList<Restraint> restraintList = new ArrayList<Restraint>()
            restraintList.add(restraint)
            RestraintGroup restraintGroup = new RestraintGroup(restraint, 1)
            restraintGroups.add(restraintGroup)
        }

        private ArrayList<RestraintGroup> getRestraintGroups() {
            return restraintGroups
        }

        private int getNumRestraintGroups() {
            return restraintGroups.size()
        }
    }

    /**
     * Contains the collections, groups and restraints for a system. Controls the addition/updating of restraints
     * during a simulation.
     */
    private class MeldRestraintTransformer {
        PointerByReference meldForce
        ArrayList<SelectivelyActiveCollection> selectivelyActiveCollections
        ArrayList<AlwaysOnRestraint> alwaysActiveRestraints

        private MeldRestraintTransformer(PointerByReference meldForce, ArrayList<SelectivelyActiveCollection> selectivelyActiveCollections, ArrayList<AlwaysOnRestraint> alwaysActiveRestraints) {
            this.meldForce = meldForce
            this.selectivelyActiveCollections = selectivelyActiveCollections
            this.alwaysActiveRestraints = alwaysActiveRestraints
        }

        private addInteractions() {
            for (int collectionInd = 0; collectionInd < selectivelyActiveCollections.size(); collectionInd++) {
                SelectivelyActiveCollection collection = selectivelyActiveCollections.get(collectionInd)
                PointerByReference meldGroupIndices = OpenMM_IntArray_create(0)
                ArrayList<RestraintGroup> restraintGroups = collection.getRestraintGroups()
                for (int groupInd = 0; groupInd < collection.getNumRestraintGroups(); groupInd++) {
                    RestraintGroup group = restraintGroups.get(groupInd)
                    PointerByReference meldRestraintIndices = OpenMM_IntArray_create(0)
                    ArrayList<Restraint> restraints = group.getRestraints()
                    for (int restraintInd = 0; restraintInd < group.getNumRestraints(); restraintInd++) {
                        Restraint restraint = restraints.get(restraintInd)
                        int meldRestraintIndex = (int) addMeldRestraint(meldForce, restraint, 0.0, 0.0)
                        OpenMM_IntArray_append(meldRestraintIndices, meldRestraintIndex)
                    }
                    int meldGroupIndex = MeldOpenMMLibrary.OpenMM_MeldForce_addGroup(meldForce, meldRestraintIndices, group.numActive)
                    OpenMM_IntArray_append(meldGroupIndices, meldGroupIndex)
                    OpenMM_IntArray_destroy(meldRestraintIndices)
                }
                MeldOpenMMLibrary.OpenMM_MeldForce_addCollection(meldForce, meldGroupIndices, collection.numActive)
                OpenMM_IntArray_destroy(meldGroupIndices)
            }
            //TODO: Add the meld force to the system and return the system.
        }

        private update(float alpha, float timestep) {
            for (int collectionInd = 0; collectionInd < selectivelyActiveCollections.size(); collectionInd++) {
                SelectivelyActiveCollection collection = selectivelyActiveCollections.get(collectionInd)
                ArrayList<RestraintGroup> restraintGroups = collection.getRestraintGroups()
                for (int groupInd = 0; collection.getNumRestraintGroups(); groupInd++) {
                    RestraintGroup group = restraintGroups.get(groupInd)
                    ArrayList<Restraint> restraints = group.getRestraints()
                    for (int restraintInd = 0; restraintInd < group.getNumRestraints(); restraintInd++) {
                        Restraint restraint = restraints.get(restraintInd)
                        updateMeldRestraint(meldForce, restraint, alpha, timestep)
                    }
                }
            }
            //TODO: updateParametersInContext
            //MeldOpenMMLibrary.OpenMM_MeldForce_updateParametersInContext(meldForce, context)
        }
    }

    /**
     * This method adds a MELD restraint to the system.
     * @param meldForce A PointerByReference containing the meld forces.
     * @param restraint A Restraint.
     * @param alpha A float between [0,1] indicating the lambda value (MC-OST) or temperature (Replica Exchange)
     * of the system.
     * @param timestep An integer indicating the molecular dynamics timestep.
     * @return The index of the restraint.
     */
    private int addMeldRestraint(PointerByReference meldForce, Restraint restraint, float alpha, float timestep) {
        int restIndex
        if (restraint instanceof DistanceRestraint) {
            float scale = restraint.scaler.call(alpha) * restraint.ramp.call(timestep)
            float scaledForceConstant = restraint.distanceForceConstant * scale
            restIndex = MeldOpenMMLibrary.OpenMM_MeldForce_addDistanceRestraint(meldForce, restraint.alphaCIndex, restraint.alphaCPlus3Index, restraint.r1, restraint.r2, restraint.r3, restraint.r4, scaledForceConstant)
        } else if (restraint instanceof TorsionRestraint) {
            float scale = restraint.scaler.call(alpha) * restraint.ramp.call(timestep)
            float scaledForceConstant = restraint.torsionForceConstant * scale
            restIndex = MeldOpenMMLibrary.OpenMM_MeldForce_addTorsionRestraint(meldForce, restraint.atom1Index, restraint.atom2Index, restraint.atom3Index, restraint.atom4Index, restraint.angle, restraint.deltaAngle, scaledForceConstant)
        } else {
            logger.severe("Restraint type cannot be handled.")
        }
        return restIndex
    }

    /**
     * This method updates a MELD restraint in the system.
     * @param meldForce A PointerByReference containing the meld forces.
     * @param restraint A Restraint.
     * @param alpha A float between [0,1] indicating the lambda value (MC-OST) or temperature (Replica Exchange)
     * of the system.
     * @param timestep An integer indicating the molecular dynamics timestep.
     * @param index The index of the restraint.
     * @return The new index.
     */
    private int updateMeldRestraint(PointerByReference meldForce, Restraint restraint, float alpha, float timestep, int index) {
        if (restraint instanceof DistanceRestraint) {
            float scale = restraint.scaler.call(alpha) * restraint.ramp.call(timestep)
            float scaledForceConstant = restraint.distanceForceConstant * scale
            MeldOpenMMLibrary.OpenMM_MeldForce_modifyDistanceRestraint(meldForce, restraint.alphaCIndex, restraint.alphaCPlus3Index, restraint.r1, restraint.r2, restraint.r3, restraint.r4, scaledForceConstant)
            index++
        } else if (restraint instanceof TorsionRestraint) {
            float scale = restraint.scaler.call(alpha) * restraint.ramp.call(timestep)
            float scaledForceConstant = restraint.torsionForceConstant * scale
            MeldOpenMMLibrary.OpenMM_MeldForce_modifyTorsionRestraint(meldForce, restraint.atom1Index, restraint.atom2Index, restraint.atom3Index, restraint.atom4Index, restraint.angle, restraint.deltaAngle, scaledForceConstant)
            index++
        } else {
            logger.severe("Restraint type cannot be handled.")
        }
        return index
    }

    private abstract class RestraintScaler{
        float alphaMin = 0.0
        float alphaMax = 1.0

        void checkAlphaRange(alpha){
            if(alpha<0 || alpha>1.0) {
                logger.severe("Alpha must be in range [0,1].")
            }
        }

        abstract float call(float alpha)
    }

    private class ConstantScaler extends RestraintScaler{
        private ConstantScaler(){
        }

        float call(float alpha){
            checkAlphaRange(alpha)
            return 1.0
        }
    }

    private class NonLinearScaler extends RestraintScaler{
        float factor
        float strengthAtAlphaMax = 0.0
        float strengthAtAlphaMin = 1.0

        private NonLinearScaler(float factor){
            this.factor = factor
            if(factor<1){
                logger.severe("Factor must be greater than 1.")
            }
        }

        private float handleBoundaries(float alpha){
            if(alpha<=alphaMin){
                return 1.0
            } else if (alpha>=alphaMax){
                return 0.0
            } else{
                return null
            }
        }

        float call(float alpha){
            checkAlphaRange(alpha)
            float scale = handleBoundaries(alpha)
            if(scale.is(null)){
                float delta = (alpha - alphaMin)/(alphaMax - alphaMin)
                float norm = 1.0/(Math.exp(factor)-1.0)
                scale = norm * (Math.exp(factor*(1.0-delta))-1.0)
            }
            scale = (1.0-scale)*(strengthAtAlphaMax-strengthAtAlphaMin) + strengthAtAlphaMin
            return scale
        }
    }

    private class LinearRamp{
        float startTime
        float endTime
        float startWeight
        float endWeight

        private LinearRamp(float startTime, float endTime, float startWeight, float endWeight){
            this.startTime = startTime
            this.endTime = endTime
            this.startWeight = startWeight
            this.endWeight = endWeight
        }

        float call(float timeStep){
            if(timeStep < 0){
                logger.severe("Timestep cannot be less than 0.")
            }
            if(timeStep<startTime){
                return startWeight
            }else if(timeStep<endTime){
                return startWeight + (endWeight-startWeight) * (timeStep - startTime) / (endTime - startTime)
            }else{
                return endWeight
            }
        }
    }

    @Override
    List<Potential> getPotentials() {
        return forceFieldEnergy == null ? Collections.emptyList() : Collections.singletonList(forceFieldEnergy)
    }
}

