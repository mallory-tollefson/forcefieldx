/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.nonbonded;

import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;

/**
 * This class implements a spatial decomposition based on partitioning a
 * grid into octants.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class SpatialDensityRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(SpatialDensityRegion.class.getName());
    /**
     * The number of divisions along the A-axis.
     */
    private int nA;
    /**
     * The number of divisions along the B-axis.
     */
    private int nB;
    /**
     * The number of divisions along the C-Axis.
     */
    private int nC;
    /**
     * The number of cells in one plane (nDivisions^2).
     */
    private int nAB;
    /**
     * The number of cells (nDivisions^3).
     */
    private final int nCells;
    /**
     * Number of octant work cells.
     */
    private final int nWork;
    /**
     * Number of octant work cells with at least one atom 
     * (actualWork is less than or equal to nWork).
     */
    private int actualWork;
    /**
     * A temporary array that holds the index of the cell each atom is assigned
     * to.
     */
    private final int cellIndex[][];
    /**
     * The cell indices of each atom along a A-axis.
     */
    private final int cellA[];
    /**
     * The cell indices of each atom along a B-axis.
     */
    private final int cellB[];
    /**
     * The cell indices of each atom along a C-axis.
     */
    private final int cellC[];
    /**
     * The A index of each octant (0..nA - 1) that may not have any atoms.
     */
    protected final int workA[];
    /**
     * The B index of each octant (0..nB - 1) that may not have any atoms.
     */
    protected final int workB[];
    /**
     * The C index of each octant (0..nC - 1) that may not have any atoms.
     */
    protected final int workC[];
    /**
     * The A index of each octant (0..nA - 1) that may not have any atoms.
     */
    protected final int actualA[];
    /**
     * The B index of each octant (0..nB - 1) that may not have any atoms.
     */
    protected final int actualB[];
    /**
     * The C index of each octant (0..nC - 1) that may not have any atoms.
     */
    protected final int actualC[];
    /**
     * The list of atoms in each cell. [nsymm][natom] = atom index
     */
    protected final int cellList[][];
    /**
     * The offset of each atom from the start of the cell. The first atom atom
     * in the cell has 0 offset. [nsymm][natom] = offset of the atom
     */
    protected final int cellOffset[][];
    /**
     * The number of atoms in each cell. [nsymm][ncell]
     */
    protected final int cellCount[][];
    /**
     * The index of the first atom in each cell. [nsymm][ncell]
     */
    protected final int cellStart[][];
    public int nSymm;
    private final double coordinates[][][];
    private final Atom atoms[];
    private final double xf[];
    private final double yf[];
    private final double zf[];
    private final Crystal crystal;
    private final int nAtoms;
    private final int gridSize;
    private double grid[] = null;
    private float floatGrid[] = null;
    private double initValue = 0.0;
    private SpatialDensityLoop spatialDensityLoop[];
    private GridInitLoop gridInitLoop;

    public SpatialDensityRegion(int gX, int gY, int gZ, double grid[],
                                int basisSize, int nSymm, int minWork,
                                int threadCount, Crystal crystal,
                                Atom atoms[], double coordinates[][][]) {
        this(gX, gY, gZ, basisSize, nSymm, minWork, threadCount, crystal, atoms, coordinates);
        this.grid = grid;
    }

    public SpatialDensityRegion(int gX, int gY, int gZ, float grid[],
                                int basisSize, int nSymm, int minWork,
                                int threadCount, Crystal crystal,
                                Atom atoms[], double coordinates[][][]) {
        this(gX, gY, gZ, basisSize, nSymm, minWork, threadCount, crystal, atoms, coordinates);
        this.floatGrid = grid;
    }

    private SpatialDensityRegion(int gX, int gY, int gZ,
                                 int basisSize, int nSymm, int minWork,
                                 int threadCount, Crystal crystal,
                                 Atom atoms[], double coordinates[][][]) {
        /**
         * Chop up the 3D unit cell domain into fractional coordinate chunks to
         * allow multiple threads to put charge density onto the grid without
         * needing the same grid point. First, we partition the X-axis, then
         * the Y-axis, and finally the Z-axis if necesary.
         */
        this.crystal = crystal;
        this.coordinates = coordinates;
        this.nSymm = nSymm;
        this.atoms = atoms;
        this.nAtoms = atoms.length;

        gridInitLoop = new GridInitLoop();
        gridSize = gX * gY * gZ * 2;

        xf = new double[nAtoms];
        yf = new double[nAtoms];
        zf = new double[nAtoms];

        int nX = gX / basisSize;
        int nY = gY / basisSize;
        int nZ = gZ / basisSize;
        int div = 1;
        int currentWork = 0;
        if (threadCount > 1 && nZ > 1) {
            if (nZ % 2 != 0) {
                nZ--;
            }
            nC = nZ;
            div = 2;
            currentWork = nC / div / threadCount;
            // If we have enough work per thread, stop dividing the domain.
            if (currentWork > minWork || nY < 2) {
                // Reduce the number of divisions along the X-axis if possible
                while (currentWork >= minWork) {
                    nC -= 2;
                    currentWork = nC / div / threadCount;
                }
                nC += 2;
                nA = 1;
                nB = 1;
            } else {
                if (nY % 2 != 0) {
                    nY--;
                }
                nB = nY;
                div = 4;
                currentWork = nB * nC / div / threadCount;
                // If we have 4 * threadCount * minWork chunks, stop dividing the domain.
                if (currentWork > minWork || nX < 2) {
                    while (currentWork >= minWork) {
                        nB -= 2;
                        currentWork = nB * nC / div / threadCount;
                    }
                    nB += 2;
                    nA = 1;
                } else {
                    if (nX % 2 != 0) {
                        nX--;
                    }
                    nA = nX;
                    div = 8;
                    currentWork = nA * nB * nC / div / threadCount;
                    while (currentWork >= minWork) {
                        nA -= 2;
                        currentWork = nA * nB * nC / div / threadCount;
                    }
                    nA += 2;
                }
            }
            nAB = nA * nB;
            nCells = nAB * nC;
            nWork = nCells / div;
        } else {
            nA = 1;
            nB = 1;
            nC = 1;
            nAB = 1;
            nCells = 1;
            nWork = 1;
        }

        logger.info(String.format(" Grid chunks per thread:    %d / %d = %8.3f\n",
                                  nWork, threadCount, ((double) nWork) / threadCount));
        workA = new int[nWork];
        workB = new int[nWork];
        workC = new int[nWork];
        actualA = new int[nWork];
        actualB = new int[nWork];
        actualC = new int[nWork];
        int index = 0;
        for (int h = 0; h < nA; h += 2) {
            for (int k = 0; k < nB; k += 2) {
                for (int l = 0; l < nC; l += 2) {
                    workA[index] = h;
                    workB[index] = k;
                    workC[index++] = l;
                }
            }
        }
        cellList = new int[nSymm][nAtoms];
        cellIndex = new int[nSymm][nAtoms];
        cellOffset = new int[nSymm][nAtoms];
        cellStart = new int[nSymm][nCells];
        cellCount = new int[nSymm][nCells];
        cellA = new int[nAtoms];
        cellB = new int[nAtoms];
        cellC = new int[nAtoms];
    }

    public void setDensityLoop(SpatialDensityLoop loops[]) {
        spatialDensityLoop = loops;
    }

    private class GridInitLoop extends IntegerForLoop {

        private final IntegerSchedule schedule = IntegerSchedule.fixed();
        // Extra padding to avert cache interference.
        long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        long pad8, pad9, pada, padb, padc, padd, pade, padf;

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void run(int lb, int ub) {
            if (floatGrid != null) {
                for (int i = lb; i <= ub; i++) {
                    floatGrid[i] = (float) initValue;
                }
            }
            if (grid != null) {
                for (int i = lb; i <= ub; i++) {
                    grid[i] = initValue;
                }
            }
        }
    }

    public void setInitValue(double initValue) {
        this.initValue = initValue;
    }

    @Override
    public void run() {
        int ti = getThreadIndex();
        int work1 = actualWork - 1;
        SpatialDensityLoop loop = spatialDensityLoop[ti];
        try {
            execute(0, gridSize - 1, gridInitLoop);
            execute(0, work1, loop.setOctant(0));
            // Fractional chunks along the C-axis.
            if (nC > 1) {
                execute(0, work1, loop.setOctant(1));
                // Fractional chunks along the B-axis.
                if (nB > 1) {
                    execute(0, work1, loop.setOctant(2));
                    execute(0, work1, loop.setOctant(3));
                    // Fractional chunks along the A-axis.
                    if (nA > 1) {
                        execute(0, work1, loop.setOctant(4));
                        execute(0, work1, loop.setOctant(5));
                        execute(0, work1, loop.setOctant(6));
                        execute(0, work1, loop.setOctant(7));
                    }
                }
            }
        } catch (Exception e) {
            logger.severe(e.toString());
        }
    }

    /**
     * Assign asymmetric and symmetry mate atoms to cells. This is very fast;
     * there is little to be gained from parallelizing it at this point.
     */
    public void assignAtomsToCells() {
        // Zero out the cell counts.
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int cellIndexs[] = cellIndex[iSymm];
            final int cellCounts[] = cellCount[iSymm];
            final int cellStarts[] = cellStart[iSymm];
            final int cellLists[] = cellList[iSymm];
            final int cellOffsets[] = cellOffset[iSymm];
            for (int i = 0; i < nCells; i++) {
                cellCounts[i] = 0;
            }

            // Convert to fractional coordinates.
            final double redi[][] = coordinates[iSymm];
            final double x[] = redi[0];
            final double y[] = redi[1];
            final double z[] = redi[2];
            crystal.toFractionalCoordinates(nAtoms, x, y, z, xf, yf, zf);

            // Assign each atom to a cell using fractional coordinates.
            for (int i = 0; i < nAtoms; i++) {
                double xu = xf[i];
                double yu = yf[i];
                double zu = zf[i];

                // Move the atom into the range 0.0 <= x < 1.0
                while (xu >= 1.0) {
                    xu -= 1.0;
                }
                while (xu < 0.0) {
                    xu += 1.0;
                }
                while (yu >= 1.0) {
                    yu -= 1.0;
                }
                while (yu < 0.0) {
                    yu += 1.0;
                }
                while (zu >= 1.0) {
                    zu -= 1.0;
                }
                while (zu < 0.0) {
                    zu += 1.0;
                }

                // The cell indices of this atom.
                final int a = (int) Math.floor(xu * nA);
                final int b = (int) Math.floor(yu * nB);
                final int c = (int) Math.floor(zu * nC);
                if (iSymm == 0) {
                    cellA[i] = a;
                    cellB[i] = b;
                    cellC[i] = c;
                }
                // The cell index of this atom.
                final int index = a + b * nA + c * nAB;
                cellIndexs[i] = index;
                // The offset of this atom from the beginning of the cell.
                cellOffsets[i] = cellCounts[index]++;
            }

            // Define the starting indices.
            cellStarts[0] = 0;
            for (int i = 1; i < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }

            // Move atom locations into a list ordered by cell.
            for (int i = 0; i < nAtoms; i++) {
                final int index = cellIndexs[i];
                cellLists[cellStarts[index]++] = i;
            }

            // Redefine the starting indices again.
            cellStarts[0] = 0;
            for (int i = 1; i < nCells; i++) {
                final int i1 = i - 1;
                cellStarts[i] = cellStarts[i1] + cellCounts[i1];
            }

            // Loop over work chunks and get rid of empty chunks.
            actualWork = 0;
            for (int icell = 0; icell < nWork; icell++) {
                int ia = workA[icell];
                int ib = workB[icell];
                int ic = workC[icell];
                int empty = count(ia, ib, ic);
                // Fractional chunks along the C-axis.
                if (nC > 1 && empty == 0) {
                    empty += count(ia, ib, ic + 1);
                    // Fractional chunks along the B-axis.
                    if (nB > 1 && empty == 0) {
                        empty += count(ia, ib + 1, ic);
                        empty += count(ia, ib + 1, ic + 1);
                        // Fractional chunks along the A-axis.
                        if (nA > 1 && empty == 0) {
                            empty += count(ia + 1, ib, ic);
                            empty += count(ia + 1, ib, ic + 1);
                            empty += count(ia + 1, ib + 1, ic);
                            empty += count(ia + 1, ib + 1, ic + 1);
                        }
                    }
                }
                // If there is work in this chunk, include it.
                if (empty > 0) {
                    actualA[actualWork] = ia;
                    actualB[actualWork] = ib;
                    actualC[actualWork++] = ic;
                }
            }
        }
    }

    private int count(int ia, int ib, int ic) {
        int count = 0;
        for (int iSymm = 0; iSymm < nSymm; iSymm++) {
            final int index = index(ia, ib, ic);
            final int start = cellStart[iSymm][index];
            final int stop = start + cellCount[iSymm][index];
            count += (stop - start);
        }
        return count;
    }

    public int index(int ia, int ib, int ic) {
        return ia + ib * nA + ic * nAB;
    }
}
