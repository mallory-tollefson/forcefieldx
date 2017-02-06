/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.potential.extended;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.reduction.SharedDouble;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtUtils.SB;
import ffx.potential.nonbonded.MultiplicativeSwitch;
import ffx.potential.parameters.MultipoleType;

import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.TitrationESV.TitrationUtils.isTitratableHydrogen;

/**
 * A generalized extended system variable.
 * Treatment of ESVs:
 *  a. Bonded terms interpolate linearly between end states.
 *      ESVs based on MultiResidue (e.g. TitrationESV) place a multiplier in the term objects themselves.
 *  b. PME and vdW scaling and derivatives are handled inside these classes' inner loops.
 *      Softcoring follows the same form used by OSRW lambda.
 * @author slucore
 */
public abstract class ExtendedVariable {

    // System handles
    private static final Logger logger = Logger.getLogger(ExtendedVariable.class.getName());
    private static int esvIndexer = 0;
    public final int index;
    private boolean ready = false;

    // Properties
    private static final boolean esvPropagation = prop("esv-propagation", false);
    private static final Double biasOverride = prop("esv-biasOverride", Double.NaN);
    private static final double thetaMass = prop("esv-thetaMass", 1.0e-18);            // from OSRW, reasonably 100 a.m.u.
    private static final double thetaFriction = prop("esv-thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps

    // Lambda and derivative variables
    private double lambda;                          // ESVs travel on {0,1}
    private double theta;                           // Propagates lambda particle via "lambda=sin(theta)^2"
    private double halfThetaVelocity = 0.0;         // from OSRW, start theta with zero velocity
    private final Random stochasticRandom = ThreadLocalRandom.current();
    /**
     * Magnitude of the discretization bias in kcal/mol.
     */
    private final double discrBiasMag;
    /**
     * Discretization bias and its (chain rule) derivative.
     */
    private double discrBias, dDiscrBiasdL;
    /**
     * Switched lambda value and its (chain rule) derivative.
     */
    private double lSwitch, dlSwitch;
    /**
     * Sigmoidal switching function. Maps lambda -> S(lambda) which has a flatter deriv near zero/unity.
     * Properties: {S(0)=0, S(1)=1, dS(0)=0, dS(1)=0, S(1-L)=1-S(L)}.
     */
    private final MultiplicativeSwitch switchingFunction;

    /* Atom lists and scaled terms                  */
    private final MultiResidue multiRes;
    private Residue residueForeground;              // resOne*lamedh + resZero*(1-lamedh)
    private Residue residueBackground;
    private final List<Atom> backbone;
    private final List<Atom> atomsForeground;       // (L); side-chain only
    private final List<Atom> atomsBackground;       // (1-L); side-chain only; permanently disconnected from assembly
    private final List<Atom> atomsShared;           // all foreground atoms except titrating hydrogens
    private final List<Atom> atomsUnshared;         // titrating (and thus foreground) atoms
//    private final List<Atom> masterAtomList;        // shared + unshared + background
    private final int moleculeNumber;

    /* Bonded energy and derivative handling        */
    private List<BondedTerm> bondedFg, bondedBg;  // valence terms for each side; mola won't see zro by default
    private MSNode termNode;                        // modified to contain all applicable bonded terms
    private final SharedDouble bondedDeriv = new SharedDouble();    // bonded dUdL reduction target
    private final HashMap<Class<? extends BondedTerm>,SharedDouble> fgBondedDerivDecomp;    // foreground dUdL by term
    private final HashMap<Class<? extends BondedTerm>,SharedDouble> bgBondedDerivDecomp;    // background dUdL by term
    private final HashMap<Atom,Atom> fg2bg = new HashMap<>();   // maps multipole end points of this ESV's lambda path

    public ExtendedVariable(MultiResidue multiRes, double biasMag, double initialLambda) {
        index = esvIndexer++;
        discrBiasMag = Double.isFinite(biasOverride) ? biasOverride : biasMag;
        this.switchingFunction = new MultiplicativeSwitch(0.0, 1.0);
        setLambda(initialLambda);

        this.multiRes = multiRes;
        residueForeground = multiRes.getActive();
        termNode = residueForeground.getTerms();
        residueBackground = multiRes.getInactive().get(0);
        moleculeNumber = residueForeground.getAtomList().get(0).getMoleculeNumber();

        backbone = new ArrayList<>();
        atomsForeground = new ArrayList<>();
        atomsBackground = new ArrayList<>();
        atomsShared = new ArrayList<>();
        atomsUnshared = new ArrayList<>();
        
        if (ExtendedSystem.esvDecomposeBonded) {
            fgBondedDerivDecomp = new HashMap<>();
            bgBondedDerivDecomp = new HashMap<>();
        } else {
            fgBondedDerivDecomp = null;
            bgBondedDerivDecomp = null;
        }
    }

    public ExtendedVariable(MultiResidue multiRes, double biasMag) {
        this(multiRes, biasMag, 1.0);
    }
    public ExtendedVariable(MultiResidue multiRes) {
        this(multiRes, 0.0, 1.0);
    }

    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     */
    public abstract double getTotalBias(double temperature, boolean print);
    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     */
    public abstract double getTotalBiasDeriv(double temperature, boolean print);

    /**
     * Propagate lambda using Langevin dynamics.
     * Check that temperature goes to the value used below (when set as a constant),
     * even when sim is decoupled. Be sure it call setLambda() rather than using
     * direct access for array resizing, etc.
     */
    public void propagate(double dEdEsv, double dt, double setTemperature) {
        if (!esvPropagation) {
            return;
        }
        double rt2 = 2.0 * ExtConstants.Boltzmann * setTemperature * thetaFriction / dt;
        double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / ExtConstants.forceToKcal;
        double dEdL = -dEdEsv * sin(2.0 * theta);
        halfThetaVelocity = (halfThetaVelocity * (2.0 * thetaMass - thetaFriction * dt)
                + ExtConstants.forceToKcalSquared * 2.0 * dt * (dEdL + randomForce))
                / (2.0 * thetaMass + thetaFriction * dt);
        theta = theta + dt * halfThetaVelocity;

        if (theta > PI) {
            theta -= 2.0 * PI;
        } else if (theta <= -PI) {
            theta += 2.0 * PI;
        }

        double sinTheta = sin(theta);
        setLambda(sinTheta * sinTheta);
    }

    /**
     * BEWARE calling this manually; a call to ExtendedSystem.updateListeners() is required 
     * before this change can properly take effect.
     */
    public final void setLambda(double lambda) {
        this.lambda = lambda;
        this.lSwitch = switchingFunction.taper(lambda);
        this.dlSwitch = switchingFunction.dtaper(lambda);
        theta = Math.asin(Math.sqrt(lambda));
        discrBias = discrBiasMag - 4*discrBiasMag*(lambda-0.5)*(lambda-0.5);
        dDiscrBiasdL = -8*discrBiasMag*(lambda-0.5);
        if (ready) {
            updateMultipoleTypes();
        }
    }
    
    private void updateMultipoleTypes() {
        for (Atom fg : atomsShared) {
            Atom bg = fg2bg.get(fg);
            MultipoleType Ptype = fg.getMultipoleType();
            MultipoleType Utype =  bg.getMultipoleType();
            MultipoleType types[] = new MultipoleType[]{Ptype, Utype};
            double mWeights[] = new double[]{lambda, 1.0 - lambda};
            double mdotWeights[] = new double[]{1.0, -1.0};
            MultipoleType multipoleM;
            MultipoleType multipoleMdot;
            if (Ptype == null) {
                // Multipoles not yet defined.
                continue;
            }
            
            if (Utype == null) {
                SB.logfn("Programming error @updateMultipoleTypes(): bgType was null.");
                SB.logfn("   fg: %s, %s", fg.toString(), fg.getName());
                SB.logfn(" Background atoms available for match: ");
                for (Atom debug : atomsBackground) {
                    SB.logfn("   bg: %s, %s", debug.toString(), debug.getName());
                }
                SB.warning();
                continue;
            }
            int frameTypes[] = Ptype.frameAtomTypes;
            multipoleM = MultipoleType.scale(types, mWeights, frameTypes);
            multipoleMdot = MultipoleType.scale(types, mdotWeights, frameTypes);
            if (multipoleM == null) {
                logger.warning("Programming error @updateMultipoleTypes(): combo was null.");
                continue;
            }
            if (multipoleMdot == null) {
                logger.warning("Programming error @updateMultipoleTypes(): combo was null.");
                continue;
            }
            fg.setEsvMultipoleM(multipoleM);
            fg.setEsvMultipoleMdot(multipoleMdot);
            SB.logf(" Assigning ESV MultipoleTypes for atom %s", fg.toNameNumberString());
            SB.nlogf("  U: %.2f*%s", lambda, Ptype.toCompactBohrString());
            SB.nlogf("  P: %.2f*%s", 1.0 - lambda, Utype.toCompactBohrString());
            SB.nlogf("  M:      %s", multipoleM.toCompactBohrString());
            SB.nlogf("  Mdot:   %s", multipoleMdot.toCompactBohrString());
            SB.print();
        }
    }

    /**
     * The unswitched lambda value, ie input to S(L).
     * This is probably not what you want.
     */
    public final double getLambda() {
        return lambda;      // L
    }

    public final double getLambdaSwitch() {
        return lSwitch;     // S(L)
    }

    public final double getSwitchDeriv() {
        return dlSwitch;    // dS(L)dL
    }

    public final int getIndex() {
        return index;
    }

    public String getName() {
        return String.format("ESV_%d", index);
    }

    @Override
    public String toString() {
        return String.format("ESV_%d (%4.2f->%4.2f)",
                index, getLambda(), getLambdaSwitch());
    }

    /**
     * From Shen&Huang 2016; drives ESVs to zero/unity.
     * bias = 4B*(L-0.5)^2
     */
    public double getDiscrBias() {
        return discrBias;
    }

    /**
     * dBiasdL = -8B*(L-0.5)
     */
    public double getDiscrBiasDeriv() {
        return dDiscrBiasdL;
    }

    /**
     * Fill the atom arrays; apply persistent indexing; set atom esv properties;
     * fill the bonded term arrays; set esv lambda on bonded terms.
     */
    public void readyup() {
        if (ready) {
            logger.warning("Program error @readyup(): should not be reinvoked.");
            return;
        }
        // Fill the atom lists.
        for (String bbName : ExtConstants.backboneNames) {
            Atom bb = (Atom) residueForeground.getAtomNode(bbName);
            if (bb != null) {
                backbone.add(bb);
//                bb.applyPersistentIndex();
            }
        }
        for (Atom atom : residueForeground.getAtomList()) {
            if (!backbone.contains(atom)) {
                atomsForeground.add(atom);
//                SB.logfn(" searching for fg mate: '%s'", atom.getName());
//                Atom companionAtom = residueBackground.getNodeAsAtom(atom.getName());
                Atom companionAtom = residueBackground.getAtomList().stream()
                        .filter((Atom bg) -> {
//                            SB.logfn(" querying bg atom: '%s'", bg.getName());
                            return bg.getName().equalsIgnoreCase(atom.getName());
                                })
                        .findFirst().orElse(null);
                SB.print();
                if (companionAtom == null) {
                    atomsUnshared.add(atom);
                    atom.setEsvState(1);
                    /* The following check ought to be safely removable if you've
                     * defined ExtendedVariables that are not TitrationESVs.      */
                    assert(isTitratableHydrogen(atom));
                    if (!isTitratableHydrogen(atom)) {
                        logger.warning(format("%s could not identify a companion "
                                + "for foreground atom %s.", this.toString(), atom));
                        throw new IllegalStateException();
                    }
                } else {
                    atomsShared.add(atom);
                    fg2bg.put(atom, companionAtom);
                }
            }
        }
        for (Atom a0 : residueBackground.getAtomList()) {
            if (!backbone.contains(a0)) {
                if (atomsForeground.contains(a0)) {
                    logger.severe("Yup.");
                }
                assert(!atomsForeground.contains(a0));
                assert(!isTitratableHydrogen(a0));
                atomsBackground.add(a0);
                a0.sendToBackground();
            }
        }
        
        for (Atom atom : ExtUtils.joinedListView(atomsForeground, atomsBackground)) {
            atom.setESV(this);
//            atom.applyPersistentIndex();
        }

//        TODO: Map Cartesian gradients assigned to background atoms to
//                their corresponding foreground atom if available.
//                Cancelled: background atoms are not true D.o.F. and as such do not
//                    appear in the atom array shown to VdW, PME.
//                    Instead, the input atom has pre-scaled multipoles.

        // Fill bonded term list and set all esvLambda values.
        bondedFg = residueForeground.getDescendants(BondedTerm.class);
        bondedBg = residueBackground.getDescendants(BondedTerm.class);
        MSNode extendedTermNode = new MSNode(format("Extended (%d)", bondedBg.size()));
        for (MSNode node : residueBackground.getTerms().getChildList()) {
            extendedTermNode.add(node);
        }
        multiRes.getActive().getTerms().add(extendedTermNode);

        ready = true;
        updateBondedLambdas();
        describe();
    }
    
    public Atom getBackground(Atom foreground) {
        return fg2bg.get(foreground);
    }

    /**
     * List all the atoms and bonded terms associated with each end state.
     */
    public void describe() {
        if (!ready) {
            return;
        }
        SB.logfn(" %s", this.toString());
        SB.logfn(" Switching on (%.2f,%.2f) with e0,de0,e1,de1: %.2f %.2f %.2f %.2f", 0.0, 1.0,
                switchingFunction.taper(0.0), switchingFunction.dtaper(0.0),
                switchingFunction.taper(1.0), switchingFunction.dtaper(1.0));
        SB.logfn("   Shared Atoms");
        for (Atom atom : atomsShared) {
            SB.logfn("%s", atom);
        }
        SB.logfn("   Unshared Atoms");
        for (Atom atom : atomsUnshared) {
            SB.logfn("%s", atom);
        }
        SB.logfn("   Background Atoms");
        for (Atom atom : atomsBackground) {
            SB.logfn("%s", atom);
        }
        SB.logfn("   Bonded Terms");
        for (MSNode term : residueForeground.getTerms().getChildList()) {
            SB.logfn("     %s", term);
            if (term.toString().trim().contains("Extended")) {
                for (MSNode ext : term.getChildList()) {
                    SB.logfn("       %s", ext);
                }
            }
        }
        SB.print();
    }

    public List<Atom> viewUnsharedAtoms() {
        return Collections.unmodifiableList(atomsUnshared);
    }
    
    public List<Atom> viewSharedAtoms() {
        return Collections.unmodifiableList(atomsShared);
    }
    
    public List<Atom> viewBackgroundAtoms() {
        return Collections.unmodifiableList(atomsBackground);
    }

    public boolean isReady() {
        return ready;
    }

    public void updateBondedLambdas() {
        if (!ExtendedSystem.esvScaleBonded || !ready) {
            return;
        }
        bondedDeriv.set(0.0);
        double Sl = getLambdaSwitch();
        double dSldL = getSwitchDeriv();
        for (BondedTerm bt1 : bondedFg) {
            if (ExtendedSystem.esvDecomposeBonded) {
                bt1.attachExtendedVariable(Sl, dSldL, bondedDeriv, fgBondedDerivDecomp);
                fgBondedDerivDecomp.clear();
            } else {
                bt1.attachExtendedVariable(Sl, dSldL, bondedDeriv);
            }
        }
        for (BondedTerm bt0 : bondedBg) {
            if (ExtendedSystem.esvDecomposeBonded) {
                bt0.attachExtendedVariable(1.0 - Sl, -dSldL, bondedDeriv, bgBondedDerivDecomp);
                bgBondedDerivDecomp.clear();
            } else {
                bt0.attachExtendedVariable(1.0 - Sl, -dSldL, bondedDeriv);
            }
        }
    }

    public double getBondedDeriv() {
        return bondedDeriv.get();
    }

    public HashMap<Class<? extends BondedTerm>,SharedDouble> getBondedDerivDecomp() {
        return fgBondedDerivDecomp;
    }

    public HashMap<Class<? extends BondedTerm>,SharedDouble> getBackgroundBondedDerivDecomp() {
        return bgBondedDerivDecomp;
    }

    public int getMoleculeNumber() {
        return moleculeNumber;
    }

}
