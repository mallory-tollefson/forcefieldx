// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.nonbonded.implicit;

import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Parallel computation of Born radii via the Grycuk method.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class BornRadiiRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(BornRadiiRegion.class.getName());
  private static final double oneThird = 1.0 / 3.0;
  private static final double PI4_3 = 4.0 / 3.0 * PI;
  private static final double PI_12 = PI / 12.0;
  private static final double DEFAULT_SNECK = 0.4058;
  private final BornRadiiLoop[] bornRadiiLoop;
  /** An ordered array of atoms in the system. */
  protected Atom[] atoms;
  /** Periodic boundary conditions and symmetry. */
  private Crystal crystal;
  /** Atomic coordinates for each symmetry operator. */
  private double[][][] sXYZ;
  /** Neighbor lists for each atom and symmetry operator. */
  private int[][][] neighborLists;
  /** Base radius of each atom. */
  private double[] baseRadius;
  /** Descreen radius of each atom. */
  private double[] descreenRadius;
  /**
   * Overlap scale factor for each atom, when using the Hawkins, Cramer & Truhlar pairwise
   * descreening algorithm.
   *
   * <p>G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized Models of Aqueous Free Energies
   * of Solvation Based on Pairwise Descreening of Solute Atomic Charges from a Dielectric Medium",
   * J. Phys. Chem., 100, 19824-19839 (1996).
   */
  private double[] overlapScale;
  /** Born radius of each atom. */
  private double[] born;
  /** Flag to indicate if an atom should be included. */
  private boolean[] use;
  /** GK cut-off distance squared. */
  private double cut2;
  /** Forces all atoms to be considered during Born radius updates. */
  private boolean nativeEnvironmentApproximation;
  /** If true, the descreening integral includes overlaps with the volume of the descreened atom */
  private final boolean perfectHCTScale;

  private SharedDoubleArray sharedBorn;
  private SharedDouble ecavTot;
  private boolean verboseRadii;
  private boolean neckCorrection;
  private boolean tanhCorrection;
  private double Sneck;
  private BornRescalingTanh bornRescalingTanh;

  public BornRadiiRegion(int nt, ForceField forceField, boolean perfectHCTScale) {
    bornRadiiLoop = new BornRadiiLoop[nt];
    for (int i = 0; i < nt; i++) {
      bornRadiiLoop[i] = new BornRadiiLoop();
    }
    ecavTot = new SharedDouble(0.0);
    verboseRadii = forceField.getBoolean("VERBOSE_BORN_RADII", false);
    neckCorrection = forceField.getBoolean("NECK_CORRECTION", false);
    tanhCorrection = forceField.getBoolean("TANH_CORRECTION",false);
    Sneck = forceField.getDouble("SNECK",DEFAULT_SNECK);
    this.perfectHCTScale = perfectHCTScale;
    if (verboseRadii) {
      logger.info(" Verbose Born radii.");
    }
  }

  @Override
  public void finish() {
    int nAtoms = atoms.length;
    double bigRadius = 50.0;
    for (int i = 0; i < nAtoms; i++) {
      final double baseRi = baseRadius[i];
      if (!use[i]) {
        born[i] = baseRi;
      } else {
        double sum = sharedBorn.get(i);
        if (sum <= 0.0) {
          born[i] = bigRadius;
          if (verboseRadii) {
            logger.info(
                format(" Born integral < 0 for atom %d; set Born radius to %12.6f (Base Radius: %2.6f)",
                        i+1, born[i],baseRadius[i]));
          }
        } else {
          // TODO: Fix tanh correction
          if(tanhCorrection){
            //sum = BornRescalingTanh.Tanh.rescale(sum, baseRi);
          }
          born[i] = 1.0 / pow(sum / PI4_3, oneThird);
          if (born[i] < baseRi) {
            born[i] = baseRi;
            if (verboseRadii) {
              logger.info(
                  format(" Born radius < Base Radius for atom %d: set Born radius to %12.6f", i+1,
                      baseRi));
            }
          } else if (born[i] > bigRadius) {
            born[i] = bigRadius;
            if (verboseRadii) {
              logger.info(
                  format(" Born radius > 50.0 Angstroms for atom %d: set Born radius to %12.6f", i+1,
                      baseRi));
            }
          } else if (isInfinite(born[i]) || isNaN(born[i])) {
            if (verboseRadii) {
              logger.info(
                  format(" Born radius NaN / Infinite for atom %d; set Born radius to %12.6f", i+1,
                      baseRi));
            }
            born[i] = baseRi;
          } else {
            if (verboseRadii){
              logger.info(
                      format(" Set Born radius for atom %d to %12.6f " +
                                      "(Base Radius: %2.6f)", i+1, born[i], baseRi)
              );
            }
          }
        }
      }
    }
    if (verboseRadii) {
      // This could get very verbose if printed at each step.
      logger.info(" Disabling verbose radii printing.");
      verboseRadii = false;
    }
  }

  public void init(
      Atom[] atoms,
      Crystal crystal,
      double[][][] sXYZ,
      int[][][] neighborLists,
      double[] baseRadius,
      double[] descreenRadius,
      double[] overlapScale,
      boolean[] use,
      double cut2,
      boolean nativeEnvironmentApproximation,
      double[] born) {
    this.atoms = atoms;
    this.crystal = crystal;
    this.sXYZ = sXYZ;
    this.neighborLists = neighborLists;
    this.baseRadius = baseRadius;
    this.descreenRadius = descreenRadius;
    this.overlapScale = overlapScale;
    this.use = use;
    this.cut2 = cut2;
    this.nativeEnvironmentApproximation = nativeEnvironmentApproximation;
    this.born = born;
  }

  @Override
  public void run() {
    try {
      int nAtoms = atoms.length;
      execute(0, nAtoms - 1, bornRadiiLoop[getThreadIndex()]);
    } catch (Exception e) {
      String message = "Fatal exception computing Born radii in thread " + getThreadIndex() + "\n";
      logger.log(Level.SEVERE, message, e);
    }
  }

  @Override
  public void start() {
    int nAtoms = atoms.length;
    if (sharedBorn == null || sharedBorn.length() < nAtoms) {
      sharedBorn = new SharedDoubleArray(nAtoms);
    }
    for (int i = 0; i < nAtoms; i++) {
      sharedBorn.set(i, 0.0);
    }
  }

  /**
   * Compute Born radii for a range of atoms via the Grycuk method.
   *
   * @since 1.0
   */
  private class BornRadiiLoop extends IntegerForLoop {

    private double[] localBorn;
    private double ecav;
    // Extra padding to avert cache interference.
    private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
    private long pad8, pad9, pada, padb, padc, padd, pade, padf;

    BornRadiiLoop() {
      ecav = 0.0;
    }

    @Override
    public void finish() {
      sharedBorn.reduce(localBorn, DoubleOp.SUM);
      ecavTot.addAndGet(ecav);
    }

    @Override
    public void run(int lb, int ub) {
      // The descreening integral is initialized to the limit of the atom alone in solvent.
      for (int i = lb; i <= ub; i++) {
        final double baseRi = baseRadius[i];
        localBorn[i] = PI4_3 / (baseRi * baseRi * baseRi);
      }
      int nSymm = crystal.spaceGroup.symOps.size();
      if (nSymm == 0) {
        nSymm = 1;
      }
      double[] x = sXYZ[0][0];
      double[] y = sXYZ[0][1];
      double[] z = sXYZ[0][2];
      for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
        double[][] xyz = sXYZ[iSymOp];
        for (int i = lb; i <= ub; i++) {
          if (!nativeEnvironmentApproximation && !use[i]) {
            continue;
          }
          final double baseRi = descreenRadius[i];
          final double descreenRi = descreenRadius[i];
          final double xi = x[i];
          final double yi = y[i];
          final double zi = z[i];
          int[] list = neighborLists[iSymOp][i];
          for (int k : list) {
            final double baseRk = descreenRadius[k];
            final double descreenRk = descreenRadius[k];
            assert (descreenRk > 0.0);
            if (!nativeEnvironmentApproximation && !use[k]) {
              continue;
            }
            if (i != k) {
              boolean is_12or13 = atoms[i].is_12_or_13(atoms[k]);
              final double xr = xyz[0][k] - xi;
              final double yr = xyz[1][k] - yi;
              final double zr = xyz[2][k] - zi;
              final double r2 = crystal.image(xr, yr, zr);
              if (r2 > cut2) {
                continue;
              }
              final double r = sqrt(r2);
              // Atom i being descreeened by atom k.
              double sk = overlapScale[k];
              if (sk > 0.0) {
                localBorn[i] += descreen(r, r2, baseRi, descreenRk, sk);
              }
              if(neckCorrection && !is_12or13 && !atoms[k].isHydrogen()) {
                // TODO: Add neck contribution to atom i being descreeened by atom k.
                /*if(neckDescreen(r,baseRi,baseRk) != 0.0) {
                  logger.info(
                          format("Modify local Born for atom %d with neck from %d by %2.8f",
                                  i + 1, k + 1, neckDescreen(r, baseRi, baseRk)));
                }*/
                localBorn[i] += neckDescreen(r, baseRi, baseRk);
              }

              // Atom k being descreeened by atom i.
              double si = overlapScale[i];
              if (si > 0.0) {
                localBorn[k] += descreen(r, r2, baseRk, descreenRi, si);
              }
              if(neckCorrection && !is_12or13 && !atoms[i].isHydrogen()) {
                // TODO: Add neck contribution to atom k being descreeened by atom i.
               /* if(neckDescreen(r,baseRk,baseRi) != 0.0) {
                  logger.info(format("Modify local Born for atom %d with neck from %d by %2.8f",
                          k + 1, i + 1, neckDescreen(r, baseRk, baseRi)));
                }*/
                localBorn[k] += neckDescreen(r, baseRk, baseRi);
              }

            } else if (iSymOp > 0) {
              final double xr = xyz[0][k] - xi;
              final double yr = xyz[1][k] - yi;
              final double zr = xyz[2][k] - zi;
              final double r2 = crystal.image(xr, yr, zr);
              if (r2 > cut2) {
                continue;
              }
              final double r = sqrt(r2);
              // Atom i being descreeened by atom k.
              double sk = overlapScale[k];
              if (sk > 0.0) {
                localBorn[i] += descreen(r, r2, baseRi, descreenRk, sk);
              }
              if(neckCorrection && !atoms[k].isHydrogen()) {
                // TODO: Add neck contribution to atom i being descreeened by atom k.
                localBorn[i] += neckDescreen(r, baseRi, baseRk);
              }

              // For symmetry mates, atom k is not descreeened by atom i.
            }
          }
        }
      }
    }

    @Override
    public void start() {
      int nAtoms = atoms.length;
      if (localBorn == null || localBorn.length < nAtoms) {
        localBorn = new double[nAtoms];
      }
      fill(localBorn, 0.0);
    }

    /**
     * Compute the integral of 1/r^6 over a neck region.
     *
     * @param r atomic separation.
     * @param radius base radius of the atom being descreened.
     * @param radiusK radius of the atom doing the descreening.
     * @return this contribution to the descreening integral.
     */
    private double neckDescreen(double r, double radius, double radiusK) {
      double neckIntegral;
      double radiusWater = 1.4;

      // Get Aij and Bij from Aguilar/Onufriev 2010 paper
      double[] constants = NeckIntegralOnufriev.NeckIntegralOnufrievConstants.run(radius, radiusK);

      double Aij = constants[0];
      double Bij = constants[1];

      logger.info(format("Aij: %2.8f Bij: %2.8f", Aij, Bij));
      if(r > radius + radiusK + 2*radiusWater || Bij > r){
          return 0.0;
      }

      // If a neck is formed, Aij can never be zero
      if(Aij <= 0.000000000){
        // Set to minimum in Aij matrix
        Aij = 0.0000161523;
      }

      /*logger.info(format("For atom rad. %2.3f descreened by atom rad. %2.3f, Aij: %2.8f Bij: %2.8f",
              radius,radiusK,Aij,Bij));*/
      double rMinusBij = r - Bij;
      double radiiMinusr = radius + radiusK + 2*radiusWater - r;
      double power1 = rMinusBij * rMinusBij * rMinusBij * rMinusBij;
      double power2 = radiiMinusr * radiiMinusr * radiiMinusr * radiiMinusr;

      // Use Aij and Bij to get neck integral using Equation 13 from Aguilar/Onufriev 2010 paper
      neckIntegral = Aij*power1*power2*Sneck;

      return -neckIntegral;
    }

    private double descreen(double r, double r2, double radius, double radiusK, double hctScale) {
      if (perfectHCTScale) {
        return perfectHCTIntegral(r, r2, radius, radiusK, hctScale);
      } else {
        return integral(r, r2, radius, radiusK * hctScale);
      }
    }

    /**
     * Use pairwise descreening to compute integral of 1/r^6.
     *
     * @param r atomic separation.
     * @param r2 atomic separation squared.
     * @param radius base radius of the atom being descreened.
     * @param scaledRadius scaled radius of the atom doing the descreening.
     * @return this contribution to the descreening integral.
     */
    private double integral(double r, double r2, double radius, double scaledRadius) {
      double integral = 0.0;
      // Descreen only if the scaledRadius is greater than zero.
      // and atom I does not engulf atom K.
      if (scaledRadius > 0.0 && (radius < r + scaledRadius)) {
        // Atom i is engulfed by atom k.
        if (radius + r < scaledRadius) {
          final double upper = scaledRadius - r;
          integral = (PI4_3 * (1.0 / (upper * upper * upper) - 1.0 / (radius * radius * radius)));
        }

        // Upper integration bound is always the same.
        double upper = r + scaledRadius;

        // Lower integration bound depends on atoms sizes and separation.
        double lower;
        if (radius + r < scaledRadius) {
          // Atom i is engulfed by atom k.
          lower = scaledRadius - r;
        } else if (r < radius + scaledRadius) {
          // Atoms are overlapped, begin integration from ri.
          lower = radius;
        } else {
          // No overlap between atoms.
          lower = r - scaledRadius;
        }

        double l2 = lower * lower;
        double l4 = l2 * l2;
        double lr = lower * r;
        double l4r = l4 * r;
        double u2 = upper * upper;
        double u4 = u2 * u2;
        double ur = upper * r;
        double u4r = u4 * r;
        double scaledRk2 = scaledRadius * scaledRadius;
        double term =
            (3.0 * (r2 - scaledRk2) + 6.0 * u2 - 8.0 * ur) / u4r
                - (3.0 * (r2 - scaledRk2) + 6.0 * l2 - 8.0 * lr) / l4r;
        integral -= PI_12 * term;
      }
      return integral;
    }

    private double perfectHCTIntegral(double r, double r2,
        double radius, double radiusK, double perfectHCT) {
      double integral = 0.0;
      // Descreen only if the scaledRadius is greater than zero.
      // and atom I does not engulf atom K.
      if (radiusK > 0.0 && (radius < r + radiusK)) {
        // Atom i is engulfed by atom k.
        // TODO: fix double counting of overlaps
        if (radius + r < radiusK) {
          final double upper = radiusK - r;
          integral = (PI4_3 * (1.0 / (upper * upper * upper) - 1.0 / (radius * radius * radius)));
        }

        // Upper integration bound is always the same.
        double upper = r + radiusK;

        // Lower integration bound depends on atoms sizes and separation.
        double lower;
        if (radius + r < radiusK) {
          // Atom i is engulfed by atom k.
          lower = radiusK - r;
        } else if (r < radius + radiusK) {
          // Atoms are overlapped, begin integration from ri.
          // TODO: fix double counting of the overlap
          lower = radius;
        } else {
          // No overlap between atoms.
          lower = r - radiusK;
        }

        double l2 = lower * lower;
        double l4 = l2 * l2;
        double lr = lower * r;
        double l4r = l4 * r;
        double u2 = upper * upper;
        double u4 = u2 * u2;
        double ur = upper * r;
        double u4r = u4 * r;
        double scaledK2 = radiusK * radiusK;
        double term =
            (3.0 * (r2 - scaledK2) + 6.0 * u2 - 8.0 * ur) / u4r
                - (3.0 * (r2 - scaledK2) + 6.0 * l2 - 8.0 * lr) / l4r;
        integral -= PI_12 * term;
      }

      return perfectHCT * integral;
    }
  }
}
