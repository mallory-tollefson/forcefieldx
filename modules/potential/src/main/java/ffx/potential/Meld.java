// ******************************************************************************
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
// ******************************************************************************
package ffx.potential;

import static edu.uiowa.jopenmm.OpenMMLibrary.*;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_destroy;

import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMMMeldLibrary;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations;
import java.util.ArrayList;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * The Meld class adds meld restraints (distance and torsion restraints) for protein folding.
 *
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class Meld {
  private static final Logger logger = Logger.getLogger(Meld.class.getName());
  private String secondaryStructure = "";
  PointerByReference meldForce;
  MeldRestraintTransformer transformer;

  /** Execute the script. */
  Meld(CompositeConfiguration properties, MolecularAssembly molecularAssembly) {
    meldForce = OpenMMMeldLibrary.OpenMM_MeldForce_create();

    // ArrayList to hold the collections of restraints (currently hydrophobic and secondary
    // structure restraints)
    ArrayList<SelectivelyActiveCollection> collections =
        new ArrayList<SelectivelyActiveCollection>();

    // Validate that secondary structure restraints and automatically modify as necessary.
    secondaryStructure = properties.getString("secondaryStructure", "");
    if (secondaryStructure.isEmpty()) {
      logger.severe(
          "Secondary structure string is empty. Must supply secondary structure prediction.");
    }
    secondaryStructure = validateSecondaryStructurePrediction(molecularAssembly);
    checkForAppropriateResidueIdentities(molecularAssembly);

    // Set up MELD restraints for secondary structure. Force constants and quadratic cut values were
    // set
    // by the Dr. Ken Dill Research Group
    SelectivelyActiveCollection collectionSecondary = new SelectivelyActiveCollection();
    double startWeight = properties.getDouble("MeldRampStartWeight", 0.0);
    LinearRamp ramp = new LinearRamp(0.0, 100.0, startWeight, 1.0);

    // The Dill Group uses a constant scaler for secondary structure restraints. We have decided to
    // use no scaler instead.
    // ConstantScaler constantScaler = new ConstantScaler();
    NoScaler noScaler = new NoScaler();
    setUpSecondaryMeldRestraints(
        molecularAssembly, noScaler, ramp, 2.48, 2.48, 2, collectionSecondary);
    collections.add(collectionSecondary);

    // Set up hydrophobic contact restraints.
    SelectivelyActiveCollection collectionHydrophobic = new SelectivelyActiveCollection();
    // Factor of 4.0 is the factor used by the Dill Group.
    // The Dill group uses a non-linear scaler for hydrophobic restraints. We have decided to use no
    // scaler instead.
    // NonLinearScaler nonLinearScaler = new NonLinearScaler(4.0);
    double contactsPerHydrophobe = 1.3;
    setUpHydrophobicRestraints(
        molecularAssembly, noScaler, ramp, collectionHydrophobic, contactsPerHydrophobe);
    collections.add(collectionHydrophobic);

    // Likely create array of restraintEnergies for the meldForce here.
    ArrayList<AlwaysOnRestraint> alwaysActiveRestraints = new ArrayList<>();
    transformer = new MeldRestraintTransformer(collections, alwaysActiveRestraints);
    transformer.addInteractions();
  }

  public PointerByReference getMeldForce() {
    return meldForce;
  }

  /**
   * This method sets distance restraints between atoms of hydrophobic residues.
   *
   * @param molecularAssembly The molecular assembly.
   * @param scaler Scales the value that controls the strength of the restraint force constants.
   * @param ramp Slowly turns restraints on/off during the simulation.
   * @param collectionHydrophobic A SelectivelyActiveCollection of hydrophobic RestraintGroups.
   * @param contactsPerHydrophobe A float representing the number of hydrophobic contacts that occur
   *     per hydrophobic residue. The Dill Group uses 1.3.
   * @return The meld force PointerByReference
   */
  void setUpHydrophobicRestraints(
      MolecularAssembly molecularAssembly,
      RestraintScaler scaler,
      LinearRamp ramp,
      SelectivelyActiveCollection collectionHydrophobic,
      double contactsPerHydrophobe) {
    ArrayList<Residue> residueList = (ArrayList<Residue>) molecularAssembly.getResidueList();
    ArrayList<Residue> hydrophobicResidues1 = new ArrayList<>();
    for (int i = 0; i < residueList.size(); i++) {
      Residue residue = residueList.get(i);
      ResidueEnumerations.AminoAcid3 aminoAcid3 = residue.getAminoAcid3();
      String aminoAcidString = aminoAcid3.toString();

      // Hydrophobic residues identities.
      String ALAString = Residue.AA3.ALA.toString();
      String VALString = Residue.AA3.VAL.toString();
      String LEUString = Residue.AA3.LEU.toString();
      String ILEString = Residue.AA3.ILE.toString();
      String PHEString = Residue.AA3.PHE.toString();
      String TRPString = Residue.AA3.TRP.toString();
      String METString = Residue.AA3.MET.toString();
      String PROString = Residue.AA3.PRO.toString();

      if (aminoAcidString.equals(ALAString)
          || aminoAcidString.equals(VALString)
          || aminoAcidString.equals(LEUString)
          || aminoAcidString.equals(ILEString)
          || aminoAcidString.equals(PHEString)
          || aminoAcidString.equals(TRPString)
          || aminoAcidString.equals(METString)
          || aminoAcidString.equals(PROString)) {
        hydrophobicResidues1.add(residue);
      }
    }

    // If no hydrophobic residues are present, then return from the method. Else, set up hydrophobic
    // restraints.
    if (hydrophobicResidues1.isEmpty()) {
      logger.warning(
          " No hydrophobic residues in sequence. No hydrophobic restraints will be added to the system.");
    } else {
      // The pairs ArrayList keeps track of residue indices that should get hydrophobic restraints.
      ArrayList<ArrayList<Integer>> pairs = new ArrayList<>();
      for (int i = 0; i < hydrophobicResidues1.size(); i++) {
        for (int j = 0; j < hydrophobicResidues1.size(); j++) {
          if (i == j) {
            continue;
          }
          int residueNumber1 = hydrophobicResidues1.get(i).getResidueNumber();
          int residueNumber2 = hydrophobicResidues1.get(j).getResidueNumber();
          if (Math.abs(residueNumber1 - residueNumber2) < 7) {
            continue;
          }

          // Check that a pair is not being added twice.
          boolean alreadyStored = false;
          for (ArrayList<Integer> pair : pairs) {
            Integer index1 = pair.get(0);
            Integer index2 = pair.get(1);
            if ((index1.intValue() == i && index2.intValue() == j)
                || (index1.intValue() == j && index2.intValue() == i)) {
              alreadyStored = true;
            }
          }
          if (alreadyStored == false) {
            ArrayList<Integer> newPair = new ArrayList<>();
            newPair.add(i);
            newPair.add(j);
            pairs.add(newPair);
          }
        }
      }

      // restraintGroups will hold each RestraintGroup and will be added to the collection at the
      // end
      ArrayList<RestraintGroup> restraintGroups = new ArrayList<>();
      // Set up hydrophobic restraints for each pair.
      for (ArrayList<Integer> pair : pairs) {
        int index1 = pair.get(0);
        int index2 = pair.get(1);
        Residue residue1 = hydrophobicResidues1.get(index1);
        Residue residue2 = hydrophobicResidues1.get(index2);
        ArrayList<Atom> atoms1 = (ArrayList<Atom>) residue1.getAtomList();
        ArrayList<Atom> atoms2 = (ArrayList<Atom>) residue2.getAtomList();

        // Names of atoms that can have hydrophobic restraints.
        ArrayList<String> atomNamesForRestraints = new ArrayList<>();
        atomNamesForRestraints.add("CA");
        atomNamesForRestraints.add("CB");
        atomNamesForRestraints.add("CD");
        atomNamesForRestraints.add("CD1");
        atomNamesForRestraints.add("CD2");
        atomNamesForRestraints.add("CE");
        atomNamesForRestraints.add("CE1");
        atomNamesForRestraints.add("CE2");
        atomNamesForRestraints.add("CE3");
        atomNamesForRestraints.add("CG");
        atomNamesForRestraints.add("CG1");
        atomNamesForRestraints.add("CG2");
        atomNamesForRestraints.add("CG3");
        atomNamesForRestraints.add("CH2");
        atomNamesForRestraints.add("CZ");
        atomNamesForRestraints.add("CZ2");
        atomNamesForRestraints.add("CZ3");
        atomNamesForRestraints.add("NE1");
        atomNamesForRestraints.add("SD");

        // hydrophobicRestraints holds all Restraints that are created for a particular pair of
        // residues.
        ArrayList<Restraint> hydrophobicRestraints = new ArrayList<>();
        for (Atom atom1 : atoms1) {
          String name1 = atom1.getName();
          for (Atom atom2 : atoms2) {
            String name2 = atom2.getName();
            if (atomNamesForRestraints.contains(name1) && atomNamesForRestraints.contains(name2)) {
              int atom1Index = atom1.getArrayIndex();
              int atom2Index = atom2.getArrayIndex();
              // Create distance restraints between specified atoms of both residues.
              DistanceRestraint distanceRestraint =
                  new DistanceRestraint(
                      scaler, ramp, atom1Index, atom2Index, 0.0, 0.0, 0.5, 0.7, 250.0);
              hydrophobicRestraints.add(distanceRestraint);
            }
          }
        }
        if (!hydrophobicRestraints.isEmpty()) {
          RestraintGroup restraintGroup = new RestraintGroup(hydrophobicRestraints, 1);
          restraintGroups.add(restraintGroup);
        }
      }
      // Add restraintGroups to the collection and set the numActive for the collection.
      collectionHydrophobic.setRestraintGroups(restraintGroups);
      int active = (int) contactsPerHydrophobe * hydrophobicResidues1.size();
      collectionHydrophobic.setNumActive(active);
    }
  }

  /**
   * This method validates that the user-supplied secondary structure predictions are the correct
   * length and contain the correct characters.
   *
   * @param molecularAssembly The molecular assembly.
   * @return String containing the validated and edited secondary structure restraints.
   */
  String validateSecondaryStructurePrediction(MolecularAssembly molecularAssembly) {
    // The only characters that should be present in secondary structure restraint string are 'H'
    // for helix, 'E'
    // for beta sheet and '.' for coil.
    if (!secondaryStructure.matches("^[HE.]+") && !secondaryStructure.isEmpty()) {
      logger.severe(" Secondary structure restraints may only contain characters 'H', 'E' and '.'");
    }

    int numResidues = molecularAssembly.getResidueList().size();
    int numSecondaryStructure = secondaryStructure.length();

    // Only one secondary structure restraint should exist per residue.
    if (numSecondaryStructure == 0) {
      logger.warning(
          " No secondary structure restraints have been provided. Simulation will proceed "
              + "with all residues having random coil secondary structure restraints.");
      String randomCoil = org.apache.commons.lang3.StringUtils.leftPad("", numResidues, ".");
      return randomCoil;
    } else if (numSecondaryStructure < numResidues) {
      logger.warning(
          " Too few secondary structure restraints exist for number of residues present. "
              + "Random coil will be added to end residues without provided secondary structure restraints.");
      String extraCoil =
          org.apache.commons.lang3.StringUtils.rightPad(secondaryStructure, numResidues, '.');
      return extraCoil;
    } else if (numSecondaryStructure == numResidues) {
      logger.info(" Secondary structure restraints will be added for all residues.");
      return secondaryStructure;
    } else if (numSecondaryStructure > numResidues) {
      logger.warning(
          " Too many secondary structure restraints exist for number of residues present."
              + " Provided secondary structure restraints will be truncated.");
      String truncated = secondaryStructure.substring(0, numResidues);
      return truncated;
    } else {
      logger.severe(" Secondary structure restraints or residues do not exist.");
      return null;
    }
  }

  /**
   * This method checks that secondary structure assignments are appropriate for the residue
   * identity. ACE and NME residues do not have alpha carbons, so they are not compatible with the
   * alpha helix or beta sheet MELD restraints.
   *
   * @param molecularAssembly The molecular assembly.
   */
  void checkForAppropriateResidueIdentities(MolecularAssembly molecularAssembly) {
    ArrayList<Residue> residues = (ArrayList<Residue>) molecularAssembly.getResidueList();
    for (int i = 0; i < secondaryStructure.length(); i++) {
      Residue residue = residues.get(i);
      ResidueEnumerations.AminoAcid3 aminoAcid3 = residue.getAminoAcid3();

      String aminoAcidString = aminoAcid3.toString();
      String NMEString = Residue.AA3.NME.toString();
      String ACEString = Residue.AA3.ACE.toString();

      if (aminoAcidString.equals(NMEString) || aminoAcidString.equals((ACEString))) {
        char character = secondaryStructure.charAt(i);
        if (character == 'H') {
          logger.info(
              " Secondary structure was modified to accommodate non-standard amino acid residue.");
          secondaryStructure =
              secondaryStructure.substring(0, i) + '.' + secondaryStructure.substring(i + 1);
        } else if (character == 'E') {
          logger.info(
              " Secondary structure was modified to accommodate non-standard amino acid residue.");
          secondaryStructure =
              secondaryStructure.substring(0, i) + '.' + secondaryStructure.substring(i + 1);
        }
      }
    }
  }

  /**
   * This method sets up MELD torsion and distance restraints for secondary structure elements and
   * adds them to groups and collections.
   *
   * @param molecularAssembly The molecular assembly.
   * @param scaler Scales the value that controls the strength of the restraint force constants.
   * @param ramp Slowly turns restraints on/off during the simulation.
   * @param torsionForceConstant A float with value supplied by The Dill Group. In kJ/mol/(10
   *     degree)^2
   * @param distanceForceConstant A float with value supplied by The Dill Group. In
   *     kJ/mol/Angstrom^2.
   * @param quadraticCut A float with value supplied by The Dill Group. This tells where to begin
   *     having a quadratic bottom rather than a flat bottom on the force. In Angstroms.
   * @param collection
   */
  void setUpSecondaryMeldRestraints(
      MolecularAssembly molecularAssembly,
      RestraintScaler scaler,
      LinearRamp ramp,
      double torsionForceConstant,
      double distanceForceConstant,
      double quadraticCut,
      SelectivelyActiveCollection collection) {
    torsionForceConstant /= 100;
    distanceForceConstant *= 100; // Convert to kJ/mol/nm^2
    quadraticCut *= 10; // Convert to nm.
    int minNumResForSecondary = 4;
    int runLength = 5;

    char helixChar = 'H';
    char sheetChar = 'E';

    ArrayList<ArrayList<Integer>> helices =
        extractSecondaryElement(secondaryStructure, helixChar, minNumResForSecondary, runLength);
    ArrayList<ArrayList<Integer>> sheets =
        extractSecondaryElement(secondaryStructure, sheetChar, minNumResForSecondary, runLength);

    double phi = -62.5;
    double deltaPhi = 17.5;
    double psi = -42.5;
    double deltaPsi = 17.5;
    ArrayList<RestraintGroup> helixRestraintGroupList = new ArrayList<>();
    setRestraintsByElementType(
        molecularAssembly,
        scaler,
        ramp,
        helixRestraintGroupList,
        helixChar,
        helices,
        phi,
        deltaPhi,
        psi,
        deltaPsi,
        quadraticCut,
        torsionForceConstant,
        distanceForceConstant);

    phi = -117.5;
    deltaPhi = 27.5;
    psi = 145.0;
    deltaPsi = 25.0;
    ArrayList<RestraintGroup> sheetRestraintGroupList = new ArrayList<>();
    setRestraintsByElementType(
        molecularAssembly,
        scaler,
        ramp,
        sheetRestraintGroupList,
        sheetChar,
        sheets,
        phi,
        deltaPhi,
        psi,
        deltaPsi,
        quadraticCut,
        torsionForceConstant,
        distanceForceConstant);

    // Set up the collection.
    ArrayList<RestraintGroup> allRestraintsGroupList = new ArrayList<>();
    allRestraintsGroupList.addAll(helixRestraintGroupList);
    allRestraintsGroupList.addAll(sheetRestraintGroupList);
    double keepDouble = allRestraintsGroupList.size() * 1.0;
    int keep = (int) keepDouble;
    collection.setNumActive(keep);
    collection.setRestraintGroups(allRestraintsGroupList);
  }

  /**
   * This method determines the starting and ending indices for secondary elements of the requested
   * type based on the user-supplied secondary structure restraint predictions. The Dill Group
   * requires that secondary elements are 5 residues long (runLength) and at least 4 residues must
   * match the requested secondary element type (minNumResidues) in order to be a secondary element.
   *
   * @param ss A string of the secondary structure prediction.
   * @param elementType Character indicating type of secondary element being searched (helix, coil,
   *     sheet).
   * @param minNumResidues Integer minimum of matching secondary structure predictions within a run
   *     length to create a secondary element. The Dill Group uses 3.
   * @param runLength The minimum number of residues necessary to create a secondary element. The
   *     Dill Group uses 5.
   * @return ArrayList<ArrayList < Integer>> that contains the starting and ending indices for
   *     helices or sheets.
   */
  ArrayList<ArrayList<Integer>> extractSecondaryElement(
      String ss, char elementType, int minNumResidues, int runLength) {
    // Will hold starting and ending indices for all found secondary elements of the requested type.
    ArrayList<ArrayList<Integer>> allElements = new ArrayList<>();
    // Iterate through the secondary structure prediction string in increments of 5 residues.
    for (int i = 0; i <= ss.length() - runLength; i++) {
      String substring = ss.substring(i, i + runLength);
      // MatchSum holds the number of times that the 5 residue segment matches the requested element
      // type.
      int matchSum = 0;
      for (int j = 0; j < runLength; j++) {
        if (substring.charAt(j) == elementType) {
          matchSum++;
        }
      }
      // If enough residues in the 5 residue segment match the requested element type, then store
      // the starting
      // and ending indices for the secondary element.
      if (matchSum >= minNumResidues) {
        ArrayList<Integer> currentElement = new ArrayList<>();
        currentElement.add(i);
        currentElement.add(i + runLength);
        allElements.add(currentElement);
      }
    }
    return allElements;
  }

  /**
   * This method sets the torsion and distance restraints for a secondary element type. For example,
   * if the starting and ending points for all helices are provided, then the torsion and distance
   * restraints appropriate for helices will be set.
   *
   * @param molecularAssembly The molecular assembly.
   * @param scaler Scales the value that controls the strength of the restraint force constants.
   * @param ramp Slowly turns restraints on/off during the simulation.
   * @param restraintGroupList An empty ArrayList of RestraintGroup objects that gets filled during
   *     the method.
   * @param elementType A character indicating the type of secondary structure being searched for (H
   *     for helices, E for sheets)
   * @param elements An ArrayList<ArrayList<Integer>> that contains the starting and ending points
   *     for secondary structure elements (helices, sheets, coils, etc).
   * @param phi A float with value supplied by The Dill Group.
   * @param deltaPhi A float with value supplied by The Dill Group.
   * @param psi A float with value supplied by The Dill Group.
   * @param deltaPsi A float with value supplied by The Dill Group.
   * @param quadraticCut A float with value supplied by The Dill Group.
   * @param torsionForceConstant A float with value supplied by The Dill Group.
   * @param distanceForceConstant A float with value supplied by The Dill Group.
   */
  void setRestraintsByElementType(
      MolecularAssembly molecularAssembly,
      RestraintScaler scaler,
      LinearRamp ramp,
      ArrayList<RestraintGroup> restraintGroupList,
      char elementType,
      ArrayList<ArrayList<Integer>> elements,
      double phi,
      double deltaPhi,
      double psi,
      double deltaPsi,
      double quadraticCut,
      double torsionForceConstant,
      double distanceForceConstant) {
    if (!elements.isEmpty()) {
      for (int i = 0; i < elements.size(); i++) { // For every secondary element.
        ArrayList<Restraint> restraintList = new ArrayList<>();
        ArrayList<Integer> element = elements.get(i);
        int elementStart = element.get(0);
        int elementEnd = element.get(1);
        ArrayList<Residue> residues = (ArrayList<Residue>) molecularAssembly.getResidueList();

        // TORSION RESTRAINTS
        for (int j = elementStart + 1;
            j < elementEnd - 1;
            j++) { // For every residue in an element.
          Residue residueJminus1 = residues.get(j - 1);
          Residue residueJ = residues.get(j);
          Residue residueJplus1 = residues.get(j + 1);

          Atom carbon = (Atom) residueJ.getAtomNode("C");
          Atom alphaC = (Atom) residueJ.getAtomNode("CA");
          Atom nitrogen = (Atom) residueJ.getAtomNode("N");
          Atom lastCarbon = (Atom) residueJminus1.getAtomNode("C");
          Atom nextNitrogen = (Atom) residueJplus1.getAtomNode("N");

          // Phi torsion restraint.
          try {
            int lastCarbonIndex = lastCarbon.getArrayIndex();
            int nitrogenIndex = nitrogen.getArrayIndex();
            int alphaCIndex = alphaC.getArrayIndex();
            int carbonIndex = carbon.getArrayIndex();
            TorsionRestraint torsionRestraint =
                new TorsionRestraint(
                    scaler,
                    ramp,
                    lastCarbonIndex,
                    nitrogenIndex,
                    alphaCIndex,
                    carbonIndex,
                    phi,
                    deltaPhi,
                    torsionForceConstant);
            restraintList.add(torsionRestraint);
          } catch (Exception e) {
            logger.severe(" Meld phi torsion restraint cannot be created.\n");
          }

          // Psi torsion restraint.
          try {
            int nitrogenIndex = nitrogen.getArrayIndex();
            int alphaCIndex = alphaC.getArrayIndex();
            int carbonIndex = carbon.getArrayIndex();
            int nextNitrogenIndex = nextNitrogen.getArrayIndex();
            TorsionRestraint torsionRestraint =
                new TorsionRestraint(
                    scaler,
                    ramp,
                    nitrogenIndex,
                    alphaCIndex,
                    carbonIndex,
                    nextNitrogenIndex,
                    psi,
                    deltaPsi,
                    torsionForceConstant);
            restraintList.add(torsionRestraint);
          } catch (Exception e) {
            logger.severe(" Meld psi torsion restraint cannot be created.\n");
          }
        }

        // DISTANCE RESTRAINTS
        Residue elementStartRes = residues.get(elementStart);
        Residue elementStartResPlus1 = residues.get(elementStart + 1);
        Residue elementStartResPlus3 = residues.get(elementStart + 3);
        Residue elementStartResPlus4 = residues.get(elementStart + 4);

        Atom alphaC = (Atom) elementStartRes.getAtomNode("CA");
        Atom alphaCPlus1 = (Atom) elementStartResPlus1.getAtomNode("CA");
        Atom alphaCPlus3 = (Atom) elementStartResPlus3.getAtomNode("CA");
        Atom alphaCPlus4 = (Atom) elementStartResPlus4.getAtomNode("CA");

        // The four floats (r1-r4) are in nanometers and were provided by the Dill Research Group
        try {
          int alphaCIndex = alphaC.getArrayIndex();
          int alphaCPlus3Index = alphaCPlus3.getArrayIndex();
          double r1;
          double r2;
          double r3;
          double r4;
          if (elementType == 'H') {
            r1 = 0; // Units are nanometers for all r constants.
            r2 = 0.485;
            r3 = 0.561;
            r4 = 0.561 + quadraticCut;
          } else if (elementType == 'E') {
            r1 = 0; // Units are nanometers for all r constants.
            r2 = 0.785;
            r3 = 1.063;
            r4 = 1.063 + quadraticCut;
          } else {
            r1 = -1;
            r2 = -1;
            r3 = -1;
            r4 = -1;
            logger.severe("Cannot create restraints without alpha helix or beta sheet present.");
          }
          DistanceRestraint distanceRestraint =
              new DistanceRestraint(
                  scaler,
                  ramp,
                  alphaCIndex,
                  alphaCPlus3Index,
                  r1,
                  r2,
                  r3,
                  r4,
                  distanceForceConstant);
          restraintList.add(distanceRestraint);
        } catch (Exception e) {
          logger.severe(" Meld distance restraint cannot be created.\n");
        }

        try {
          int alphaCPlus1Index = alphaCPlus1.getArrayIndex();
          int alphaCPlus4Index = alphaCPlus4.getArrayIndex();
          double r1;
          double r2;
          double r3;
          double r4;
          if (elementType == 'H') {
            r1 = 0; // Units are nanometers for all r constants.
            r2 = 0.485;
            r3 = 0.561;
            r4 = 0.561 + quadraticCut;
          } else if (elementType == 'E') {
            r1 = 0; // Units are nanometers for all r constants.
            r2 = 0.785;
            r3 = 1.063;
            r4 = 1.063 + quadraticCut;
          } else {
            r1 = -1;
            r2 = -1;
            r3 = -1;
            r4 = -1;
            logger.severe("Cannot create restraints without alpha helix or beta sheet present.");
          }
          DistanceRestraint distanceRestraint =
              new DistanceRestraint(
                  scaler,
                  ramp,
                  alphaCPlus1Index,
                  alphaCPlus4Index,
                  r1,
                  r2,
                  r3,
                  r4,
                  distanceForceConstant);
          restraintList.add(distanceRestraint);
        } catch (Exception e) {
          logger.severe(" Meld distance restraint cannot be created.\n");
        }

        try {
          int alphaCIndex = alphaC.getArrayIndex();
          int alphaCPlus4Index = alphaCPlus4.getArrayIndex();
          double r1;
          double r2;
          double r3;
          double r4;
          if (elementType == 'H') {
            r1 = 0; // Units are nanometers for all r constants.
            r2 = 0.581;
            r3 = 0.684;
            r4 = 0.684 + quadraticCut;
          } else if (elementType == 'E') {
            r1 = 0; // Units are nanometers for all r constants.
            r2 = 1.086;
            r3 = 1.394;
            r4 = 1.394 + quadraticCut;
          } else {
            r1 = -1;
            r2 = -1;
            r3 = -1;
            r4 = -1;
            logger.severe("Cannot create restraints without alpha helix or beta sheet present.");
          }
          DistanceRestraint distanceRestraint =
              new DistanceRestraint(
                  scaler,
                  ramp,
                  alphaCIndex,
                  alphaCPlus4Index,
                  r1,
                  r2,
                  r3,
                  r4,
                  distanceForceConstant);
          restraintList.add(distanceRestraint);
        } catch (Exception e) {
          logger.severe(" Meld distance restraint cannot be created.\n");
        }
        RestraintGroup restraintGroup = new RestraintGroup(restraintList, restraintList.size());
        restraintGroupList.add(restraintGroup);
      }
    }
  }

  /**
   * This method adds a MELD restraint to the system.
   *
   * @param restraint A Restraint.
   * @param alpha A float between [0,1] indicating the lambda value (MC-OST) or temperature (Replica
   *     Exchange) of the system.
   * @param timestep An integer indicating the molecular dynamics timestep.
   * @return The index of the restraint.
   */
  private int addMeldRestraint(Restraint restraint, double alpha, double timestep) {
    int restIndex;
    if (restraint instanceof DistanceRestraint) {
      DistanceRestraint distanceRestraint = (DistanceRestraint) restraint;
      double scale = distanceRestraint.scaler.call(alpha) * distanceRestraint.ramp.call(timestep);
      double scaledForceConstant = ((DistanceRestraint) restraint).distanceForceConstant * scale;
      restIndex =
          OpenMMMeldLibrary.OpenMM_MeldForce_addDistanceRestraint(
              meldForce,
              distanceRestraint.alphaCIndex,
              distanceRestraint.alphaCPlus3Index,
              (float) distanceRestraint.r1,
              (float) distanceRestraint.r2,
              (float) distanceRestraint.r3,
              (float) distanceRestraint.r4,
              (float) scaledForceConstant);
    } else if (restraint instanceof TorsionRestraint) {
      TorsionRestraint torsionRestraint = (TorsionRestraint) restraint;
      double scale = torsionRestraint.scaler.call(alpha) * torsionRestraint.ramp.call(timestep);
      double scaledForceConstant = torsionRestraint.torsionForceConstant * scale;
      restIndex =
          OpenMMMeldLibrary.OpenMM_MeldForce_addTorsionRestraint(
              meldForce,
              torsionRestraint.atom1Index,
              torsionRestraint.atom2Index,
              torsionRestraint.atom3Index,
              torsionRestraint.atom4Index,
              (float) torsionRestraint.angle,
              (float) torsionRestraint.deltaAngle,
              (float) scaledForceConstant);
    } else {
      restIndex = 0;
      logger.severe("Restraint type cannot be handled.");
    }
    return restIndex;
  }

  /**
   * This method updates a MELD restraint in the system.
   *
   * @param restraint A Restraint.
   * @param alpha A float between [0,1] indicating the lambda value (MC-OST) or temperature (Replica
   *     Exchange) of the system.
   * @param timestep An integer indicating the molecular dynamics timestep.
   * @param distanceIndex The index of the distance restraint.
   * @param torsionIndex The index of the torsion restraint.
   * @return The new index.
   */
  private int updateMeldRestraint(
      Restraint restraint, double alpha, double timestep, int distanceIndex, int torsionIndex) {
    if (restraint instanceof DistanceRestraint) {
      DistanceRestraint distanceRestraint = (DistanceRestraint) restraint;
      double scale = distanceRestraint.scaler.call(alpha) * distanceRestraint.ramp.call(timestep);
      double scaledForceConstant = distanceRestraint.distanceForceConstant * scale;
      double smallScaledForceConstant = scaledForceConstant * 0.5;
      OpenMMMeldLibrary.OpenMM_MeldForce_modifyDistanceRestraint(
          meldForce,
          distanceIndex,
          distanceRestraint.alphaCIndex,
          distanceRestraint.alphaCPlus3Index,
          (float) distanceRestraint.r1,
          (float) distanceRestraint.r2,
          (float) distanceRestraint.r3,
          (float) distanceRestraint.r4,
          (float) smallScaledForceConstant);
      distanceIndex++;
      return distanceIndex;
    } else if (restraint instanceof TorsionRestraint) {
      TorsionRestraint torsionRestraint = (TorsionRestraint) restraint;
      double scale = torsionRestraint.scaler.call(alpha) * torsionRestraint.ramp.call(timestep);
      double scaledForceConstant = torsionRestraint.torsionForceConstant * scale;
      double smallScaledForceConstant = scaledForceConstant * 0.5;
      OpenMMMeldLibrary.OpenMM_MeldForce_modifyTorsionRestraint(
          meldForce,
          torsionIndex,
          torsionRestraint.atom1Index,
          torsionRestraint.atom2Index,
          torsionRestraint.atom3Index,
          torsionRestraint.atom4Index,
          (float) torsionRestraint.angle,
          (float) torsionRestraint.deltaAngle,
          (float) smallScaledForceConstant);
      torsionIndex++;
      return torsionIndex;
    } else {
      logger.severe("Restraint type cannot be handled.");
      return 0;
    }
  }

  /** Restraint objects are added to the meldForce. */
  private class Restraint {
    PointerByReference meldForce;
  }

  /** SelectableRestraints can be turned on or off throughout the simulation. */
  private class SelectableRestraint extends Restraint {}

  /** AlwaysOnRestraints cannot be turned on/off during the simulation. */
  private class AlwaysOnRestraint extends Restraint {}

  /** DistanceRestraints restrain the distance between two atoms. */
  private class DistanceRestraint extends SelectableRestraint {
    RestraintScaler scaler;
    LinearRamp ramp;
    int alphaCIndex;
    int alphaCPlus3Index;
    double r1; // nanometers
    double r2; // nanometers
    double r3; // nanometers
    double r4; // nanometers
    double distanceForceConstant; // kJ/mol/nm^2

    private DistanceRestraint(
        RestraintScaler scaler,
        LinearRamp ramp,
        int alphaCIndex,
        int alphaCPlus3Index,
        double r1,
        double r2,
        double r3,
        double r4,
        double distanceForceConstant) {
      this.scaler = scaler;
      this.ramp = ramp;
      this.alphaCIndex = alphaCIndex;
      this.alphaCPlus3Index = alphaCPlus3Index;
      this.r1 = r1;
      this.r2 = r2;
      this.r3 = r3;
      this.r4 = r4;
      this.distanceForceConstant = distanceForceConstant;
    }
  }

  /** TorsionRestraints restrain the angle created by selected atoms. */
  private class TorsionRestraint extends SelectableRestraint {
    RestraintScaler scaler;
    LinearRamp ramp;
    int atom1Index;
    int atom2Index;
    int atom3Index;
    int atom4Index;
    double angle; // Degrees
    double deltaAngle; // Degrees
    double torsionForceConstant; // kJ/mol/Degree^2

    private TorsionRestraint(
        RestraintScaler scaler,
        LinearRamp ramp,
        int atom1Index,
        int atom2Index,
        int atom3Index,
        int atom4Index,
        double angle,
        double deltaAngle,
        double torsionForceConstant) {
      this.scaler = scaler;
      this.ramp = ramp;
      this.atom1Index = atom1Index;
      this.atom2Index = atom2Index;
      this.atom3Index = atom3Index;
      this.atom4Index = atom4Index;
      this.angle = angle;
      this.deltaAngle = deltaAngle;
      this.torsionForceConstant = torsionForceConstant;
    }
  }

  /**
   * A RestraintGroup is a collection of Restraints that belong together (for example, all
   * DistanceRestraints and TorsionRestraints that are part of a single alpha helix).
   */
  private class RestraintGroup {
    ArrayList<Restraint> restraints;
    int numActive;
    int numRestraints;

    private RestraintGroup(ArrayList<Restraint> restraints, int numActive) {
      if (restraints.isEmpty()) {
        logger.severe(" Restraint list cannot be empty.");
      }
      if (numActive < 0) {
        logger.severe(" NumActive cannot be less than 0.");
      }
      this.restraints = restraints;
      this.numActive = numActive;
      this.numRestraints = restraints.size();
      if (numActive > numRestraints) {
        logger.severe(" Number of active restraints must be <= number of total restraints.");
      }
    }

    private int getNumRestraints() {
      return numRestraints;
    }

    private ArrayList<Restraint> getRestraints() {
      return restraints;
    }
  }

  /** A collection of multiple RestraintGroups. */
  private class SelectivelyActiveCollection {
    ArrayList<RestraintGroup> restraintGroups;
    int numActive;

    private SelectivelyActiveCollection() {}

    private SelectivelyActiveCollection(ArrayList<RestraintGroup> restraintGroups, int numActive) {
      this.restraintGroups = restraintGroups;
      this.numActive = numActive;

      if (restraintGroups.isEmpty()) {
        logger.severe("SelectivelyActiveCollection cannot have empty restraint list.");
      }
      int numGroups = restraintGroups.size();
      if (numActive > numGroups) {
        logger.severe("Number of active restraint groups must be less than number of groups.");
      }
    }

    private void setRestraintGroups(ArrayList<RestraintGroup> restraintGroups) {
      this.restraintGroups = restraintGroups;
    }

    private void setNumActive(int numActive) {
      this.numActive = numActive;
    }

    private void addRestraint(RestraintGroup restraintGroup) {
      restraintGroups.add(restraintGroup);
    }

    private void addRestraint(Restraint restraint) {
      ArrayList<Restraint> restraintList = new ArrayList<>();
      restraintList.add(restraint);
      RestraintGroup restraintGroup = new RestraintGroup(restraintList, 1);
      restraintGroups.add(restraintGroup);
    }

    private ArrayList<RestraintGroup> getRestraintGroups() {
      return restraintGroups;
    }

    private int getNumRestraintGroups() {
      return restraintGroups.size();
    }
  }

  /**
   * Contains the collections, groups and restraints for a system. Controls the addition/updating of
   * restraints during a simulation.
   */
  public class MeldRestraintTransformer {
    ArrayList<SelectivelyActiveCollection> selectivelyActiveCollections;
    ArrayList<AlwaysOnRestraint> alwaysActiveRestraints;

    private MeldRestraintTransformer(
        ArrayList<SelectivelyActiveCollection> selectivelyActiveCollections,
        ArrayList<AlwaysOnRestraint> alwaysActiveRestraints) {
      this.selectivelyActiveCollections = selectivelyActiveCollections;
      this.alwaysActiveRestraints = alwaysActiveRestraints;
    }

    private void addInteractions() {
      for (int collectionInd = 0;
          collectionInd < selectivelyActiveCollections.size();
          collectionInd++) {
        SelectivelyActiveCollection collection = selectivelyActiveCollections.get(collectionInd);
        PointerByReference meldGroupIndices = OpenMM_IntArray_create(0);
        ArrayList<RestraintGroup> restraintGroups = collection.getRestraintGroups();
        for (int groupInd = 0; groupInd < collection.getNumRestraintGroups(); groupInd++) {
          RestraintGroup group = restraintGroups.get(groupInd);
          PointerByReference meldRestraintIndices = OpenMM_IntArray_create(0);
          ArrayList<Restraint> restraints = group.getRestraints();
          for (int restraintInd = 0; restraintInd < group.getNumRestraints(); restraintInd++) {
            Restraint restraint = restraints.get(restraintInd);
            int meldRestraintIndex = (int) addMeldRestraint(restraint, 0.0, 0.0);
            OpenMM_IntArray_append(meldRestraintIndices, meldRestraintIndex);
          }
          int meldGroupIndex =
              OpenMMMeldLibrary.OpenMM_MeldForce_addGroup(
                  meldForce, meldRestraintIndices, group.numActive);
          OpenMM_IntArray_append(meldGroupIndices, meldGroupIndex);
          OpenMM_IntArray_destroy(meldRestraintIndices);
        }
        OpenMMMeldLibrary.OpenMM_MeldForce_addCollection(
            meldForce, meldGroupIndices, collection.numActive);
        OpenMM_IntArray_destroy(meldGroupIndices);
      }
    }

    public void update(double alpha, double timestep) {
      int counter = 0;
      int distanceIndex = 0;
      int torsionIndex = 0;
      for (int collectionInd = 0;
          collectionInd < selectivelyActiveCollections.size();
          collectionInd++) {
        SelectivelyActiveCollection collection = selectivelyActiveCollections.get(collectionInd);
        ArrayList<RestraintGroup> restraintGroups = collection.getRestraintGroups();
        for (int groupInd = 0; groupInd < collection.getNumRestraintGroups(); groupInd++) {
          RestraintGroup group = restraintGroups.get(groupInd);
          ArrayList<Restraint> restraints = group.getRestraints();
          for (int restraintInd = 0; restraintInd < group.getNumRestraints(); restraintInd++) {
            Restraint restraint = restraints.get(restraintInd);
            double time = System.nanoTime();
            if (restraint instanceof DistanceRestraint) {
              distanceIndex =
                  updateMeldRestraint(restraint, alpha, timestep, distanceIndex, torsionIndex);
            } else if (restraint instanceof TorsionRestraint) {
              torsionIndex =
                  updateMeldRestraint(restraint, alpha, timestep, distanceIndex, torsionIndex);
            }
            double elapsedTime = (System.nanoTime() - time) / (1E9);
            // logger.info(String.format("Elapsed time for updateMeldRestraintMethod: %f",
            // elapsedTime));
            counter++;
          }
        }
      }
      // logger.info(String.format("Number of restraints being updated: %d", counter));
    }
  }

  /** Abstract class that allows force constants on restraints to be scaled. */
  private abstract class RestraintScaler {
    double alphaMin = 0.0;
    double alphaMax = 1.0;

    void checkAlphaRange(double alpha) {
      if (alpha < 0 || alpha > 1.0) {
        logger.severe("Alpha must be in range [0,1].");
      }
    }

    abstract double call(double alpha);
  }

  /** Scaler that simply returns alpha, or the current lambda value. */
  private class NoScaler extends RestraintScaler {
    private NoScaler() {}

    double call(double alpha) {
      checkAlphaRange(alpha);
      return alpha;
    }
  }

  /** Scaler that returns 1.0, so restraints are not scaled. */
  private class ConstantScaler extends RestraintScaler {
    private ConstantScaler() {}

    double call(double alpha) {
      checkAlphaRange(alpha);
      return 1.0;
    }
  }

  /** Scaler that scales force constants non-linearly. */
  private class NonLinearScaler extends RestraintScaler {
    double factor;
    double strengthAtAlphaMax = 0.0;
    double strengthAtAlphaMin = 1.0;

    private NonLinearScaler(double factor) {
      this.factor = factor;
      if (factor < 1) {
        logger.severe("Factor must be greater than 1.");
      }
    }

    private double handleBoundaries(double alpha) {
      if (alpha <= alphaMin) {
        return 1.0;
      } else if (alpha >= alphaMax) {
        return 0.0;
      } else {
        return -1;
      }
    }

    double call(double alpha) {
      checkAlphaRange(alpha);
      double scale = handleBoundaries(alpha);
      if (scale == -1) {
        double delta = (alpha - alphaMin) / (alphaMax - alphaMin);
        double norm = 1.0 / (Math.exp(factor) - 1.0);
        scale = norm * (Math.exp(factor * (1.0 - delta)) - 1.0);
      }
      scale = (1.0 - scale) * (strengthAtAlphaMax - strengthAtAlphaMin) + strengthAtAlphaMin;
      return scale;
    }
  }

  /**
   * The linear ramp is used to slowly turn restraints on/off at particular point in a simulation.
   * It works by multiplying the force constant of a restraint by a linear interpolation of
   * user-specified start and end weights.
   */
  private class LinearRamp {
    double startTime;
    double endTime;
    double startWeight;
    double endWeight;

    private LinearRamp(double startTime, double endTime, double startWeight, double endWeight) {
      this.startTime = startTime;
      this.endTime = endTime;
      this.startWeight = startWeight;
      this.endWeight = endWeight;
    }

    double call(double timeStep) {
      if (timeStep < 0) {;
        logger.severe("Timestep cannot be less than 0.");
      }
      if (timeStep < startTime) {
        return startWeight;
      } else if (timeStep < endTime) {
        return startWeight
            + (endWeight - startWeight) * (timeStep - startTime) / (endTime - startTime);
      } else {
        return endWeight;
      }
    }
  }
}
