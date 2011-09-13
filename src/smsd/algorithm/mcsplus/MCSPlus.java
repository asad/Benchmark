/* Copyright (C) 2006-2011  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package smsd.algorithm.mcsplus;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import smsd.algorithm.mcgregor.McGregor;
import smsd.global.TimeOut;
import smsd.tools.TimeManager;

/**
 * This class handles MCS plus algorithm which is a combination of
 * c-clique algorithm and McGregor algorithm.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.SMSDBondSensitiveTest")
public class MCSPlus {

    private boolean shouldMatchRings;
    private boolean shouldMatchBonds;

    /**
     * Default constructor added 
     */
    public MCSPlus() {
    }
    private TimeManager timeManager = null;

    /**
     * @return the timeout
     */
    protected synchronized double getTimeout() {
        return TimeOut.getInstance().getMCSPlusTimeout();
    }

    /**
     * @return the timeManager
     */
    protected synchronized TimeManager getTimeManager() {
        return timeManager;
    }

    /**
     * @param aTimeManager the timeManager to set
     */
    public synchronized void setTimeManager(TimeManager aTimeManager) {
        TimeOut.getInstance().setTimeOutFlag(false);
        timeManager = aTimeManager;
    }

    /**
     * 
     * @param ac1
     * @param ac2
     * @param shouldMatchBonds 
     * @param shouldMatchRings 
     * @return
     * @throws CDKException
     */
    protected List<List<Integer>> getOverlaps(
            IAtomContainer ac1,
            IAtomContainer ac2,
            boolean shouldMatchBonds,
            boolean shouldMatchRings)
            throws CDKException {

        this.shouldMatchBonds = shouldMatchBonds;
        this.shouldMatchRings = shouldMatchRings;

        List<List<Integer>> extendMappings = null;

//        System.err.println("ac1 : " + ac1.getAtomCount());
//        System.err.println("ac2 : " + ac2.getAtomCount());
        setTimeManager(new TimeManager());
        try {
            GenerateCompatibilityGraph gcg = new GenerateCompatibilityGraph(ac1, ac2, isMatchBonds(), isMatchRings());
            List<Integer> comp_graph_nodes = gcg.getCompGraphNodes();

            List<Integer> cEdges = gcg.getCEgdes();
            List<Integer> dEdges = gcg.getDEgdes();
//
//            System.err.println("**************************************************");
//            System.err.println("C_edges: " + cEdges.size());
//            System.err.println("D_edges: " + dEdges.size());

            BKKCKCF init = new BKKCKCF(comp_graph_nodes, cEdges, dEdges);
//            Koch init = new Koch(comp_graph_nodes, cEdges, dEdges);
            Stack<List<Integer>> maxCliqueSet = null;
            maxCliqueSet = init.getMaxCliqueSet();

//            System.err.println("Max_Cliques_Set: " + maxCliqueSet);
//            System.err.println("Best Clique Size: " + init.getBestCliqueSize());
//            System.err.println("**************************************************");


            //clear all the compatibility graph content
            gcg.clear();
            List<Map<Integer, Integer>> mappings = new ArrayList<Map<Integer, Integer>>();

            while (!maxCliqueSet.empty()) {
                Map<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();
                indexindexMapping = ExactMapping.extractMapping(comp_graph_nodes, maxCliqueSet.peek());
                mappings.add(indexindexMapping);
                maxCliqueSet.pop();
            }
            extendMappings = searchMcGregorMapping(ac1, ac2, mappings);
        } catch (IOException ex) {
            Logger.getLogger(MCSPlus.class.getName()).log(Level.SEVERE, null, ex);
        }
//        System.err.println("extendMappings: " + extendMappings.iterator().next().size() / 2);

        return extendMappings;
    }

    private List<List<Integer>> searchMcGregorMapping(
            IAtomContainer ac1,
            IAtomContainer ac2,
            List<Map<Integer, Integer>> allMCSCopy) throws IOException {

        List<List<Integer>> extendMappings = new ArrayList<List<Integer>>();

        boolean ROPFlag = true;
        for (Map<Integer, Integer> firstPassMappings : allMCSCopy) {
            Map<Integer, Integer> extendMapping = new TreeMap<Integer, Integer>(firstPassMappings);
            McGregor mgit = null;

            if (ac1.getAtomCount() > ac2.getAtomCount()) {
                mgit = new McGregor(ac1, ac2, extendMappings, isMatchBonds(), isMatchRings());
            } else {
                extendMapping.clear();
                mgit = new McGregor(ac2, ac1, extendMappings, isMatchBonds(), isMatchRings());
                ROPFlag = false;
                for (Map.Entry<Integer, Integer> map : firstPassMappings.entrySet()) {
                    extendMapping.put(map.getValue(), map.getKey());
                }
            }
//            System.out.println("\nStart McGregor search");
            //Start McGregor search
            mgit.startMcGregorIteration(mgit.getMCSSize(), extendMapping);
            extendMappings = mgit.getMappings();
            mgit = null;

            if (isTimeOut()) {
                System.err.println("\nMCSPlus hit by timeout in McGregor\n");
                break;
            }

//            System.out.println("\nSol count after MG" + extendMappings.size());
        }
        List<List<Integer>> finalMappings = setMcGregorMappings(ROPFlag, extendMappings);
//        System.out.println("After set Sol count MG" + finalMappings.size());
//        System.out.println("MCSSize " + finalMappings.iterator().next().size() + "\n");

        return finalMappings;
    }

    private List<List<Integer>> setMcGregorMappings(
            boolean RONP,
            List<List<Integer>> mappings) {
        int counter = 0;
        int mcsSize = 0;
        List<List<Integer>> finalMappings = new ArrayList<List<Integer>>();
        for (List<Integer> mapping : mappings) {
            List<Integer> indexindexMapping = new ArrayList<Integer>();
            for (int index = 0; index < mapping.size(); index += 2) {
                Integer qIndex = 0;
                Integer tIndex = 0;

                if (RONP) {
                    qIndex = mapping.get(index);
                    tIndex = mapping.get(index + 1);
                } else {
                    qIndex = mapping.get(index + 1);
                    tIndex = mapping.get(index);
                }

                if (qIndex != null && tIndex != null) {
                    indexindexMapping.add(qIndex);
                    indexindexMapping.add(tIndex);
                }
            }
            if (!indexindexMapping.isEmpty() && indexindexMapping.size() > mcsSize) {
                mcsSize = indexindexMapping.size();
                finalMappings.clear();
                counter = 0;
            }
            if (!indexindexMapping.isEmpty() && !finalMappings.contains(indexindexMapping)
                    && (indexindexMapping.size()) == mcsSize) {
                finalMappings.add(counter, indexindexMapping);
                counter++;
            }
        }
        return finalMappings;
    }

    public synchronized boolean isTimeOut() {
        if (getTimeout() > -1 && getTimeManager().getElapsedTimeInMinutes() > getTimeout()) {
            TimeOut.getInstance().setTimeOutFlag(true);
            return true;
        }
        return false;
    }

    /**
     * @return the shouldMatchRings
     */
    public boolean isMatchRings() {
        return shouldMatchRings;
    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isMatchBonds() {
        return shouldMatchBonds;
    }
}
