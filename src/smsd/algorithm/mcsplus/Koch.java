/**
 *
 * Copyright (C) 2006-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
 * You should have received index copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package smsd.algorithm.mcsplus;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [Ina Koch: Enumerating all connected maximal common subgraphs
 * in two graphs; T.Comp. Sc. (2001); vol 250; pp.
 * 1-30]
 *
 * @cdk.githash
 * @cdk.module smsd
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Koch {

    private Stack<List<Integer>> maximumCliquesSet;
    /********************************************************************/
    /*
     *T: is a set of vertices which have already been used 
     * 
     */
    private List<Integer> T;
    /*
     * C: set of vertices belonging to the current clique
     */
    private List<Integer> C;
    /*
     * S: set of vertices which are not allowed to be added
     * to C
     */
    private List<Integer> S;
    /*
     *P: is a set of vertices which <b>can</b> be added to C, because they are
     * neighbours of vertex u via <i>c-edges</i>
     */
    private Stack<Integer> P;
    /*
     *D: is a set of vertices which <b>cannot</b> be added to C, because they are
     * neighbours of vertex u via <i>d-edges</i>
     */
    private Stack<Integer> D;
    /*
     *V: stored all the vertices for the Graph G
     * V[G]
     *nodes of vector comp_graph_nodes are stored in V
     */
    private Stack<Integer> V;
    /***********************************************************************/
    private List<Integer> C_edges;
    private List<Integer> D_edges;
    private int best_clique_size;
    private List<Integer> comp_graph_nodes;
    private List<Integer> C_copy;
    private Stack<Integer> P_copy;
    private Stack<Integer> D_copy;
    private List<Integer> S_copy;

    /**
     * Creates a new instance of BronKerboschKochCliqueFinder
     * @param compGraphNodes
     * @param cEdges
     * @param dEdges  
     */
    protected Koch(
            List<Integer> compGraphNodes,
            List<Integer> cEdges,
            List<Integer> dEdges) {

        this.comp_graph_nodes = compGraphNodes;
        this.C_edges = cEdges;
        this.D_edges = dEdges;
        best_clique_size = 0;

        //Initialization maximumCliquesSet

        maximumCliquesSet = new Stack<List<Integer>>();

        T = new ArrayList<Integer>(); //Initialize the T Vector

        C = new ArrayList<Integer>();
        P = new Stack<Integer>();
        D = new Stack<Integer>();
        S = new ArrayList<Integer>();



        V = new Stack<Integer>(); //Initialization of Vector V

        int V_set_size = comp_graph_nodes.size() / 3;
        for (int a = 0; a < V_set_size; a++) {
            V.add(comp_graph_nodes.get(a * 3 + 2));
        }
        // System.out.println();

        V.add(0);

        Init_Algorithm();

    }

    private void Init_Algorithm() {

        /*
         * N[u]: set of neighbours of vertex u in Graph G
         *
         */

        List<Integer> N = new ArrayList<Integer>();

        int b = 0;

        /*
         * Let T be the set of Nodes already been used in the initialization
         *
         */

        T.clear();

        while (V.get(b) != 0) {

            // V[b] is node u, v belogs to V[G]
            int central_node = V.get(b);
            P.clear();
            D.clear();
            S.clear();
            C.clear();

            //find the neighbors of the central node from V
            N = find_neighbors(V.get(b));
            for (int c = 0; c < N.size(); c = c + 2) {
                // N[c] is node v
                //Grouping of the neighbors in S,P and D

                /*
                 * u and v are adjacent via a C-edge
                 */
                if (N.get(c + 1) == 1) {
                    if (T.contains(N.get(c))) {
                        S.add(N.get(c));
                    } else {
                        P.push(N.get(c));
                    }

                } else if (N.get(c + 1) == 2) {
                    // u and v are adjacent via a D-edge
                    D.add(N.get(c));
                }
                //find respective neighbor position in P, which is needed for the deletion from V
                //int V_size = V.size();
                int neighbor_position = -1;

                int elementAtC = N.get(c);

                for (int d = 0; d < V.size(); d++) {
                    if (elementAtC == V.elementAt(d)) {
                        neighbor_position = d;
                    }
                }

                //delete neighbor from set V
                if (neighbor_position != -1) {
                    for (int e = neighbor_position; e < V.size() - 1; e++) {
                        V.set(e, V.get(e + 1));
                    }
                    V.pop();
                    if (neighbor_position < b) {
                        b = b - 1;
                    }
                }
            }
            P.add(0);
            C.add(central_node);
            Enumerate_Cliques(C, P, D, S);
            T.add(V.get(b));
            b++;
        }
    }

    private int Enumerate_Cliques(List<Integer> C, Stack<Integer> P, Stack<Integer> D, List<Integer> S) {

        List<Integer> N = new ArrayList<Integer>();
        Stack<Integer> ut_set = new Stack<Integer>();//Defined as P' in the paper

        C_copy = new ArrayList<Integer>();
        P_copy = new Stack<Integer>();
        D_copy = new Stack<Integer>();
        S_copy = new ArrayList<Integer>();


        for (Integer I : P) {
            ut_set.add(I);
        }

        if (P.size() == 1) {
            if (S.isEmpty()) {

                //store best solutions in stack maximumCliquesSet
                int clique_size = C.size();

                if (clique_size >= best_clique_size) {
                    if (clique_size > best_clique_size) {
                        while (!maximumCliquesSet.empty()) {
                            maximumCliquesSet.pop();
                        }
                        best_clique_size = clique_size;
                        //System.out.println("Best Cliques Size: " + best_clique_size + " " + clique_size );
                    }
                    if (clique_size == best_clique_size) {
                        maximumCliquesSet.push(C);
                    }
                }
                return 0;
            }
        }

        int a = 0;

        while (ut_set.elementAt(a) != 0) {
            // P[a] is node ut
            int ui = ut_set.get(a);
            //remove ut_set[a] from P
            //find position of ut_set node in P
            int P_size = P.size();
            Integer ut_node_pos = (P_size + 1);
            for (int counter = 0; counter < P_size - 1; counter++) {
                if (P.elementAt(counter) == ui) {
                    ut_node_pos = counter;
                }
            }
            if (ut_node_pos > P_size) {
                System.out.println("Index out of bound for deletion");
            }
            //delete ut_set node in P
            for (int counter = ut_node_pos; counter < P_size - 1; counter++) {
                P.set(counter, P.get(counter + 1));
            }

            P.pop(); //POP in CPP Asad

            C_copy.clear();
            P_copy.clear();
            D_copy.clear();
            S_copy.clear();
            N.clear();


            for (Integer obj : C) {
                C_copy.add(obj);
            }

            for (Integer obj : P) {
                P_copy.add(obj);
            }
            for (Integer obj : D) {
                D_copy.add(obj);
            }
            for (Integer obj : S) {
                S_copy.add(obj);
            }


            P_copy.pop();
            N = find_neighbors(ut_set.get(a));

            int N_size = N.size();

            for (int b = 0; b < N_size; b = b + 2) {

                int D_set_size = D.size();
                int Nelement_at_b = N.get(b);

                for (int c = 0; c < D_set_size; c++) {
                    if (Nelement_at_b == D.elementAt(c)) {
                        if (N.get(b + 1) == 1) {
                            if (T.contains(Nelement_at_b)) {
                                S_copy.add(N.get(b));
                            } else {
                                P_copy.push(N.get(b));
                            }
                            //delete N[b] bzw. D[c] from set D_copy
                            int D_copy_size = D_copy.size();
                            int Nb_position = 10000;
                            for (int e = 0; e < D_copy_size; e++) {
                                if (Nelement_at_b == D_copy.elementAt(e)) {
                                    Nb_position = e;
                                }
                            }
                            for (int e = Nb_position; e < D_copy_size - 1; e++) {
                                D_copy.set(e, D_copy.get(e + 1));
                            }

                            D_copy.pop();
                        }
                    }
                }
                //find respective neighbor position in ut_set, 
                //which is needed for the deletion from ut_set
                int ut_set_size = ut_set.size();
                int neighbor_position = -1;
                for (int e = 0; e < ut_set_size; e++) {
                    if (Nelement_at_b == ut_set.elementAt(e)) {
                        neighbor_position = e;
                    }
                }
                if (neighbor_position != -1) {
                    //delete neighbor from set P
                    for (int e = neighbor_position; e < ut_set_size - 1; e++) {
                        ut_set.set(e, ut_set.get(e + 1));
                    }
                    ut_set.pop();
                    if (neighbor_position < a) {
                        a = a - 1;
                    }
                }
            }
            Stack<Integer> P_copy_N_intersec = new Stack<Integer>();
            Stack<Integer> D_copy_N_intersec = new Stack<Integer>();
            List<Integer> S_copy_N_intersec = new ArrayList<Integer>();

            int nElement = -1;

            for (int sec = 0; sec < N_size; sec = sec + 2) {

                nElement = N.get(sec);

                if (P_copy.contains(nElement)) {
                    P_copy_N_intersec.push(nElement);
                }
                if (D_copy.contains(nElement)) {
                    D_copy_N_intersec.add(nElement);
                }
                if (S_copy.contains(nElement)) {
                    S_copy_N_intersec.add(nElement);
                }
            }
            P_copy_N_intersec.add(0);
            C_copy.add(ui);

            Enumerate_Cliques(C_copy, P_copy_N_intersec, D_copy_N_intersec, S_copy_N_intersec);
            S.add(ui);
            a++;
        }
        return 0;
    }

    private List<Integer> find_neighbors(int central_node) {

        List<Integer> neighbor_vec = new ArrayList<Integer>();
        int C_edge_number = C_edges.size() / 2;

        for (int a = 0; a < C_edge_number; a++) {
            if (C_edges.get(a * 2 + 0) == central_node) {
                //          System.out.println( C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 1));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }
            if (C_edges.get(a * 2 + 1) == central_node) {
                //           System.out.println(C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 0));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }
        }

        int D_edge_number = D_edges.size() / 2;

        //System.out.println("");
        //  System.out.println("D_edges Size: "+ D_edges.size());
        for (int a = 0; a < D_edge_number; a++) {
            if (D_edges.get(a * 2 + 0) == central_node) {
                //       System.out.println( D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 1));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }
            if (D_edges.get(a * 2 + 1) == central_node) {
                //        System.out.println(D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 0));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }
        }
        return neighbor_vec;
    }

    protected int getBestCliqueSize() {
        return best_clique_size;
    }

    /**
     * 
     * @return
     */
    protected Stack<List<Integer>> getMaxCliqueSet() {
        return maximumCliquesSet;
    }
}
