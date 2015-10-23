/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package profilehmm;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.function.Consumer;

/**
 *
 * @author Gunner
 */
public class ProfileHmm {

    public enum Type {

        Insertion, Deletion, Main, Start, End
    };
    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
//    static String[] proteins = new String[]{"A", "R", "N", "C", "D", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
    static int p = 20;
    final static double e2 = 0.11920292202211755;
    static Map<String, Double> amino_acid_maps = new HashMap<>();

    public static void main(String[] args) throws IOException {
        String substr = "123456789";
        System.out.println(substr.substring(0, 3) + substr.substring(4, substr.length()));
        amino_acid_maps.put("A", 8.26);
        amino_acid_maps.put("R", 5.53);
        amino_acid_maps.put("N", 4.06);
        amino_acid_maps.put("D", 5.46);
        amino_acid_maps.put("C", 1.37);
        amino_acid_maps.put("Q", 3.93);
        amino_acid_maps.put("E", 6.74);
        amino_acid_maps.put("G", 7.08);
        amino_acid_maps.put("H", 2.27);
        amino_acid_maps.put("I", 5.94);
        amino_acid_maps.put("L", 9.66);
        amino_acid_maps.put("K", 5.83);
        amino_acid_maps.put("M", 2.41);
        amino_acid_maps.put("F", 3.86);
        amino_acid_maps.put("P", 4.71);
        amino_acid_maps.put("S", 6.58);
        amino_acid_maps.put("T", 5.34);
        amino_acid_maps.put("W", 1.09);
        amino_acid_maps.put("Y", 2.92);
        amino_acid_maps.put("V", 6.87);
        File fl = new File("fasta_input");
        Scanner sc = new Scanner(fl);
        Map<String, String> multiple_sequence = new HashMap<>();
        String current_ms = "";
        while (sc.hasNext()) {
            String line = sc.nextLine();
            if (line.toCharArray()[0] == '>') {
//                String[] split = (line.split("\\|"));
                current_ms = (line.split("\\|")[1]);
//                System.out.println((line.split("\\|"))[1]);
                multiple_sequence.put(current_ms, "");
            } else {
                multiple_sequence.put(current_ms, multiple_sequence.get(current_ms) + line);
            }
        }
//        String[][] sequences = new String[multiple_sequence.values().iterator().next().length()][multiple_sequence.keySet().size()];
//        multiple_sequence.keySet().stream().forEach((name) -> {
//
//        });
        Map<Integer, ArrayList<AcidCount>> column_values = new HashMap<>();
        Map<Integer, ArrayList<AcidCount>> column_values_main = new HashMap<>();
        /**
         * refining the states, deleting states which have 50% gaps
         */
        for (int i = 0; i < multiple_sequence.values().iterator().next().length(); i++) {
            ArrayList<AcidCount> acidCounts = new ArrayList<>();
            for (final Map.Entry<String, String> entry : multiple_sequence.entrySet()) {
                AcidCount acidCount;
                final int tmp = i;
                acidCount = acidCounts.stream()
                        .filter(p -> p.name.equals(entry.getValue().charAt(tmp) + ""))
                        .findFirst().orElse(null);
                boolean present = acidCount != null;
                if (acidCount == null) {
                    acidCount = new AcidCount();
                }
                acidCount.name = (entry.getValue().charAt(tmp) + "").toUpperCase();
                acidCount.count++;
                if (!present) {
                    acidCounts.add(acidCount);
                }
            }
            column_values.put(i, acidCounts);
            //column_values_main.put(i, acidCounts);
        }
//        System.out.println(column_values.get(1));
//        for (Map.Entry<Integer, ArrayList<AcidCount>> cv : column_values.entrySet()) {
//            for (AcidCount ac : cv.getValue()) {
//                if (ac.name.equals("-") && ac.count >= (multiple_sequence.size() * .5)) {
//                    column_values_main.remove(cv.getKey());
//                }
//            }
//        }
        int del = 0;

        for (Map.Entry<Integer, ArrayList<AcidCount>> entrySet : column_values.entrySet()) {
            ArrayList<AcidCount> b = entrySet.getValue();
            boolean flg = false;
            AcidCount searchDash = b.stream().filter((AcidCount p) -> p.name.equals("-")).findFirst().orElse(null);
            if (searchDash != null) {
                if (searchDash.count >= (multiple_sequence.size() * .5)) {
                    del++;
                    flg = !flg;
                    //have to reconstruct multiple_sequence input here
                    for (Map.Entry<String, String> entrySet1 : multiple_sequence.entrySet()) {
                        //  multiple_sequence.forEach((name, sequence) -> {
                        multiple_sequence.put(entrySet1.getKey(), entrySet1.getValue().substring(0, entrySet.getKey() - del) + "" + entrySet1.getValue().substring(entrySet.getKey() - del + 1, entrySet1.getValue().length()));
//                    });
                    }
                }
            }
            if (!flg) {
                column_values_main.put(entrySet.getKey() - del, b);
            }

//            b.forEach((AcidCount p) -> {
//                if (p.name.equals("-") && p.count >= (multiple_sequence.size() * .5)) {
//                    column_values_main.remove(a);
//
//                }
//            });
        }

        column_values = null;
        System.gc();
        System.out.println("total columns " + column_values_main.keySet().size());
//        System.out.println("total columns " + column_values_main.get(column_values_main.keySet().size()));

        /*
         Main calculation
         */
        int stateCount = 0;
        final int totalSequenceCount = multiple_sequence.size();
        States startState = new States(0, Type.Start);
        List<States> AllStateList = new ArrayList<States>();
        List<States> parentListTmp = new ArrayList<>();
        parentListTmp.add(startState);
        AllStateList.add(startState);
        //insert state starting
        States insertState = new States(stateCount + 1, Type.Insertion);
        insertState.parent.put(insertState, e2);
        insertState.parent.put(startState, e2);
        startState.children.put(insertState, e2);//delete
        parentListTmp.add(insertState);
        AllStateList.add(insertState);
        for (Map.Entry<Integer, ArrayList<AcidCount>> entrySet : column_values_main.entrySet()) {
            Integer a = entrySet.getKey();
            ArrayList<AcidCount> b = entrySet.getValue();

//        column_values_main.forEach((a, b) -> {
            int prev = a - 1, next = a + 1;

            States state = new States(stateCount + a, Type.Main);

            AcidCount tmpAcidCount = column_values_main.get(a).stream().filter(p -> p.name.equals("-")).findFirst().orElse(null);
            int gap_in_self = (tmpAcidCount != null) ? tmpAcidCount.count : 0;
            //emission part

            for (Map.Entry<String, Double> entrySet1 : amino_acid_maps.entrySet()) {
                AcidCount search = b.stream().filter(ac -> ac.name.equals(entrySet1.getKey())).findFirst().orElse(null);
                state.emission.put(entrySet1.getKey(), (getPseudoCounts(totalSequenceCount - gap_in_self, (search == null) ? 0 : search.count, entrySet1.getValue())));
            }
//            amino_acid_maps.forEach((acid, w) -> {
//                int count = (int) b.stream().filter((AcidCount ac) -> ac.name.equals(acid)).count();
//                state.emission.put(acid, getPseudoCounts(count, totalSequenceCount, w));
//            });
//            b.forEach((AcidCount ac) -> {
//                if (!ac.name.equals("-")) {
//                    state.emission.put(ac.name, getPseudoCounts(ac.count, totalSequenceCount, amino_acid_maps.get(ac.name))); //System.out.println(" " + ac.name + " " + ac.count + ",");
//                }
//            });

            //assigning parent
            States tmpParent = null;
            ReturnVals retvalsParent = getLG_LL_GL(multiple_sequence, a, prev);
            for (States parent : parentListTmp) {
                switch (parent.type) {
                    case Insertion:
                        state.parent.put(parent, (1 - e2));//plz be sure about it,right now its original log -2 val
                        parent.children.put(state, 1 - e2);//delete later
                        break;
                    case Main:
                        tmpParent = parent;
                        //getting counts for transition probability
                        tmpAcidCount = column_values_main.get(prev).stream().filter(p -> p.name.equals("-")).findFirst().orElse(null);
                        int gap_in_parent = (tmpAcidCount != null) ? tmpAcidCount.count : 0;
                        if (gap_in_parent == gap_in_self && gap_in_self == 0) {//no computation if no sequence has gaps
                            state.parent.put(parent, getPseudoCounts(totalSequenceCount, totalSequenceCount));
                            parent.children.put(state, getPseudoCounts(totalSequenceCount, totalSequenceCount));//delete later
                        } else {

                            state.parent.put(parent, getPseudoCounts(totalSequenceCount, retvalsParent.ll));//chek for gaps to be negeted
                            parent.children.put(state, getPseudoCounts(totalSequenceCount, retvalsParent.ll));//delete later

                            //System.out.println("lg" + retvalsParent.lg);
                        }
                        break;
                    case Start:
                        tmpParent = parent;
//                        state.parent.put(parent, getPseudoCounts(totalSequenceCount, totalSequenceCount - gap_in_self));
                        state.parent.put(parent, 1 - e2 - e2);
                        parent.children.put(state, 1 - e2 - e2);//delete later
                        break;
                    case Deletion:
//                        state.parent.put(parent, 1 - parent.parent.values().iterator().next());
                        state.parent.put(parent, 1 - e2);
                        parent.children.put(state, 1 - e2);//delete later
                        break;
                    default:
                }
            }
//            if (a == 13 || a == 11) {
//                System.out.println("processing 11");
//            }
            parentListTmp = new ArrayList<>();
            parentListTmp.add(state);
            AllStateList.add(state);

            //delete state making
//            States deleteInsertState = null;
            ReturnVals returnValsNext = getLG_LL_GL(multiple_sequence, a, next);
//            if (a < column_values_main.keySet().size() ) {
//                System.out.println("ending" + a);

            States deleteInsertState = new States(stateCount + a, Type.Deletion);
            deleteInsertState.parent.put(deleteInsertState, e2);
            deleteInsertState.parent.put(tmpParent, (tmpParent.type == Type.Start) ? e2
                    : getPseudoCounts(totalSequenceCount, retvalsParent.lg));
            tmpParent.children.put(deleteInsertState,
                    (tmpParent.type == Type.Start) ? e2 : getPseudoCounts(totalSequenceCount, retvalsParent.lg));//delete later
            parentListTmp.add(deleteInsertState);
            AllStateList.add(deleteInsertState);

//            }
            //insert state next
            deleteInsertState = new States(stateCount + 1 + a, Type.Insertion);
            deleteInsertState.parent.put(deleteInsertState, e2);
            deleteInsertState.parent.put(state, getPseudoCounts(totalSequenceCount, returnValsNext.gl));

            state.children.put(deleteInsertState,
                    getPseudoCounts(totalSequenceCount, returnValsNext.gl));//delete later
            parentListTmp.add(deleteInsertState);
            AllStateList.add(deleteInsertState);

        }
        States endState = new States(column_values_main.keySet().size(), Type.End);

        for (States s : parentListTmp) {
            switch (s.type) {
                case Insertion:
                    endState.parent.put(s, (1 - e2));//plz be sure about it,right now its original log -2 val
                    break;
                case Main:
                    endState.parent.put(s, 1 - ((1 - e2) + s.parent.values().iterator().next()));
                    break;
                case Deletion:
                    endState.parent.put(s, 1 - s.parent.values().iterator().next());
                    break;
                default:
            }
            //System.out.println("name " + s.type + s.state_no);
        }
        System.out.println("total state count " + AllStateList.size());
//        printAll(startState);
        System.out.println(printTransmissionMatrix(AllStateList));

//        System.out.println("size of last parents " + parentListTmp.size());
//        printAll(parentListTmp.get(0));
//        printAll(parentListTmp.get(1));
//        printAll(parentListTmp.get(2));
//        double p1 = .4, p2 = .1;
//        System.out.println("prob of -2 " + getProbability(-2.0));
//        System.out.println("p1 in logodds" + getLogOdds(p1));
//        System.out.println("p2 in logodds" + getLogOdds(p2));
//        System.out.println("p3 in logodds" + getLogOdds(1 - p1));
//        System.out.println("sum " + (double) (getLogOdds(p1) + getLogOdds(p2) + getLogOdds(1 - p2 - p1)));
//        System.out.println("prob " + getProbability((double) (getLogOdds(p1) + getLogOdds(p2) + getLogOdds(1 - p2 - p1))));
//        System.out.println(column_values_main.keySet().size());
        //System.out.println(multiple_sequence.values());
    }

    public static String printEmissionMatrix(List<States> states) {
        String csv = "";
        for (States s : states) {
            if (s.type == Type.Main) {
                SortedSet<String> keys = new TreeSet<String>(s.emission.keySet());
                csv += (s.type + "-" + s.state_no + ",");
                for (String acid : keys) {
                    csv += (s.emission.get(acid) + ",");
                }
                csv += "\r\n";
//                System.out.println("");
            }
        }
        return csv;
    }

    public static String printTransmissionMatrix(List<States> states) {
        String csv = "";
        String tmp = "";
        for (States s1 : states) {
            tmp += s1.type + "-" + s1.state_no + ",";
            csv += s1.type + "-" + s1.state_no + ",";
            for (States s2 : states) {
                if ((s2.parent.keySet().stream().filter(a -> a.equals(s1)).count()) > 0) {
                    csv += s2.parent.get(s1) + ",";
                } else {
                    csv += "0,";
                }
            }
            csv += ("\r\n");
        }
        return "States,"+tmp+"\r\n"+csv;
    }

    public static void printAll(States state) {
        if (state.type != Type.Main && state.type != Type.Start) {
            return;
        }
        double i = 0;
        for (Map.Entry<States, Double> entrySet : state.children.entrySet()) {
            //System.out.println("Key " + entrySet.getKey().type + "" + entrySet.getKey().state_no + " val " + entrySet.getValue());
            i += entrySet.getValue();
        }
        System.out.println(state.type + "" + state.state_no + " Total = " + i);
        if (i > 1.0) {
//                || (state.state_no > 10 && state.state_no < 13)) {
            for (Map.Entry<States, Double> entrySet : state.children.entrySet()) {
                System.out.println("Key " + entrySet.getKey().type + "" + entrySet.getKey().state_no + " val " + entrySet.getValue());
//            i += entrySet.getValue();
            }
        }
        state.children.forEach((a, b) -> {
            printAll(a);
        });

    }

    public static ReturnVals getLG_LL_GL(Map<String, String> multiple_sequence, final int a, final int b) {
        final ReturnVals retVals = new ReturnVals();
        if (multiple_sequence.entrySet().iterator().next().getValue().length() == b || b < 0) {
            retVals.ll = multiple_sequence.size();
            return retVals;
        }
//        final BiFunction<Integer, Integer, Integer> mapper = (i, j) -> {
//            return i + j;
//        };
//        multiple_sequence.forEach((name, sequence) -> {
//            //System.out.println("a " + a + sequence.toCharArray()[a] + " " + sequence.toCharArray()[b]);
//            if (sequence.toCharArray()[a] == '-' && sequence.toCharArray()[b] != '-') {
//                mapper.apply(retVals.gl, 1);
//            } else if (sequence.toCharArray()[a] != '-' && sequence.toCharArray()[b] == '-') {
//                mapper.apply(retVals.lg, 1);
//            } else if (sequence.toCharArray()[a] != '-' && sequence.toCharArray()[b] != '-') {
//                mapper.apply(retVals.ll, 1);
//            }
//        });
        for (Map.Entry<String, String> entrySet : multiple_sequence.entrySet()) {
            String key = entrySet.getKey();
            String sequence = entrySet.getValue();
            try {
                if (sequence.toCharArray()[a] == '-' && sequence.toCharArray()[b] != '-') {
                    retVals.gl++;
                } else if (sequence.toCharArray()[a] != '-' && sequence.toCharArray()[b] == '-') {
                    retVals.lg++;
                } else if ((sequence.toCharArray()[a] != '-' && sequence.toCharArray()[b] != '-')
                        || (sequence.toCharArray()[a] == '-' && sequence.toCharArray()[b] == '-')) {
                    retVals.ll++;
                }
            } catch (ArrayIndexOutOfBoundsException ex) {
                System.out.println("");
            }
        }
        return retVals;
    }

    /**
     * this is for transition
     *
     * @param sum
     * @param count
     * @return
     */
    public static double getPseudoCounts(int sum, int count) {
        return getLogOdds((double) (1 + count) / (3 + sum));

    }

    /**
     * this is for emission matrix
     *
     * @param sum
     * @param count
     * @param acid_background_prob
     * @return
     */
    public static double getPseudoCounts(int sum, int count, double acid_background_prob) {
        return getLogOdds(((double) ((double) (acid_background_prob / 100) * p + count)) / (p + sum));
//        return (double) count / sum;
    }

    public static double getLogOdds(double p) {
        return Math.log((p) / (1 - p));
    }

    public static double getLogits(int probability) {

        return Math.log(probability / (1 - probability));
    }

    public static double getProbability(double logitValue) {
        double tmp = Math.exp(logitValue);
        return (tmp) / (tmp + 1);

    }

    public static class ReturnVals {

        public int ll = 0;
        public int lg = 0;
        public int gl = 0;
    }

    public static class AcidCount {

        public String name;
        public int count = 0;
    }
}
