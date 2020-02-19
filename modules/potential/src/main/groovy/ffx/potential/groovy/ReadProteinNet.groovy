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
package ffx.potential.groovy

import ffx.potential.cli.PotentialScript

import org.biojava.nbio.structure.StructureIO

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import org.biojava.nbio.structure.Structure
import org.biojava.nbio.structure.Chain
import org.biojava.nbio.structure.Group
import org.biojava.nbio.structure.Atom
import org.biojava.nbio.structure.GroupType


/**
 * The ReadProteinNet reads a text-based protein net and outputs B-Factors
 *
 * <br>
 * Usage:
 * <br>
 * ffxc ReadProteinNet [options] &lt;filename&gt;
 */
@Command(description = " ReadProteinNet reads a text-based protein net and outputs B-Factors.", name = "ffxc ReadProteinNet")
class ReadProteinNet extends PotentialScript {

    /**
     * Argument should be one text-based ProteinNet file.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'Text-based ProteinNet file.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    String nextLine
    BufferedReader br
    int totalProteins = 0
    int badProteins = 0
    
    /**
     * Execute the script.
     */
    @Override
    ReadProteinNet run() {

        if (!init()) {
            return this
        }

        String inputPath = filenames.get(0)

        logger.info("\n Opening ProteinNet " + inputPath)
        br = new BufferedReader(new FileReader(inputPath))

        while (true) {
            LinkedHashMap<String, String> dict = readRecord()
            if (dict != null) {
                logger.info(" " + dict.get("id"))
                logger.info(" " + dict.get("primary"))
                totalProteins += 1
                try {
                    String[] idArr = dict.get("id").split("_")
                    Structure s1 = StructureIO.getStructure(idArr[0])
                    Chain c1 = s1.getChain(idArr[idArr.length - 1])
                    if (c1 != null) {
                        for (Group group : c1.getAtomGroups(GroupType.AMINOACID)) {
                            println(group)
                            Atom a = group.getAtom("CA")
                            if (a != null) {
                                println(a.getTempFactor())
                            } else {
                                println("NO C-ALPHA FOUND")
                                badProteins += 1
                                break
                            }
                        }
                    } else {
                        println("INVALID CHAIN ID: " + dict.get("id"))
                        badProteins += 1
                    }
                } catch(FileNotFoundException e) {
                    println("PDB not found: " + dict.get("id").split("_")[0])
                    badProteins += 1
                }
            } else {
                break
            }
        }

        println("Total: " + totalProteins)
        println("Errors: " + badProteins)
        double percent = badProteins / totalProteins * 100
        println("Percent: " + percent)

        return this
    }

    /**
     *
     * @return: LinkedHashMap containing PDB ID and Primary Sequence
     */
    private LinkedHashMap<String, String> readRecord() {
        LinkedHashMap<String,String> dictionary = new LinkedHashMap<>()
        while ((nextLine = br.readLine()) != null) {
            switch(nextLine) {
                case "[ID]":
                    nextLine = br.readLine()
                    dictionary.put("id",nextLine)
                    break;
                case "[PRIMARY]":
                    nextLine = br.readLine()
                    dictionary.put("primary",nextLine)
                    break;  
                case "":
                    return dictionary                  
            }
        }
        return null
    }

}
