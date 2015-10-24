/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package profilehmm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import static profilehmm.Config.e2;
import static profilehmm.Config.p;
import profilehmm.ProfileHmm.Type;

/**
 *
 * @author shad942
 */
public class HMMConstruct {

    Map<String, String> multiple_sequence = new HashMap<>();

    public HMMConstruct(Map<String, String> multiple_sequence) {
        this.multiple_sequence = multiple_sequence;
    }

    public Map<Integer, ArrayList<AcidCount>> getColumns() {
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
        return column_values_main;
    }

    public List<States> getStates() {

        Map<Integer, ArrayList<AcidCount>> column_values_main = getColumns();
        /*
         Main calculation
         */
        int stateCount = 0;
        final int totalSequenceCount = multiple_sequence.size();
        States startState = new States(0, Type.Start);
        List<States> AllStateList = new ArrayList<>();
        List<States> parentListTmp = new ArrayList<>();
        parentListTmp.add(startState);
        AllStateList.add(startState);
        //insert state starting
        States insertState = new States(stateCount + 1, Type.Insertion);
        insertState.parent.put(insertState, e2);
        insertState.parent.put(startState, e2);
        startState.children.put(insertState, e2);//delete
        setEmission(insertState, true);

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
            //emission part for main states
            for (Map.Entry<String, Double> entrySet1 : new Config().amino_acid_maps.entrySet()) {
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
            States tmpParent = null, tmpDeleteState = null;
            ReturnVals retvalsParent = getLG_LL_GL(multiple_sequence, prev, a);
            for (States parent : parentListTmp) {
                switch (parent.type) {
                    case Insertion:
                        state.parent.put(parent, (1 - e2));//plz be sure about it,right now its original log -2 val
                        parent.children.put(state, 1 - e2);//delete later
                        break;
                    case Main:
                        tmpParent = parent;
                        //getting counts for transition probability
//                        tmpAcidCount = column_values_main.get(prev).stream().filter(p -> p.name.equals("-")).findFirst().orElse(null);
//                        int gap_in_parent = (tmpAcidCount != null) ? tmpAcidCount.count : 0;
//                        if (gap_in_parent == gap_in_self && gap_in_self == 0) {//no computation if no sequence has gaps
//                            state.parent.put(parent, getPseudoCounts(totalSequenceCount, totalSequenceCount));
//                            parent.children.put(state, getPseudoCounts(totalSequenceCount, totalSequenceCount));//delete later
//                        } else {
                        state.parent.put(parent, getPseudoCounts(totalSequenceCount - retvalsParent.gg, retvalsParent.ll));//chek for gaps to be negeted
                        parent.children.put(state, getPseudoCounts(totalSequenceCount - retvalsParent.gg, retvalsParent.ll));//delete later
                        //System.out.println("lg" + retvalsParent.lg);
//                        }
                        break;
                    case Start:
                        tmpParent = parent;
//                        state.parent.put(parent, getPseudoCounts(totalSequenceCount, totalSequenceCount - gap_in_self));
                        state.parent.put(parent, 1 - e2 - e2);
                        parent.children.put(state, 1 - e2 - e2);//delete later
                        break;
                    case Deletion:

                        tmpDeleteState = parent;
//                        state.parent.put(parent, 1 - parent.parent.values().iterator().next());
                        state.parent.put(parent, 1 - e2);
                        parent.children.put(state, 1 - e2);//delete later
                        break;
                    default:
                }
            }
//            if (a == 31 ) {
//                System.out.println("processing 11");
//            }
            parentListTmp = new ArrayList<>();
            parentListTmp.add(state);
            AllStateList.add(state);

            ReturnVals returnValsNext = getLG_LL_GL(multiple_sequence, a, next);
            States deleteInsertState = new States(stateCount + a, Type.Deletion);
            if (tmpDeleteState != null) {
                deleteInsertState.parent.put(tmpDeleteState, e2);
            }
            deleteInsertState.parent.put(tmpParent, (tmpParent.type == Type.Start) ? e2
                    : getPseudoCounts(totalSequenceCount, retvalsParent.lg));
            tmpParent.children.put(deleteInsertState,
                    (tmpParent.type == Type.Start) ? e2
                            : getPseudoCounts(totalSequenceCount - retvalsParent.gg, retvalsParent.lg));//delete later
            parentListTmp.add(deleteInsertState);
            AllStateList.add(deleteInsertState);
//            }
            //insert state next
            deleteInsertState = new States(stateCount + 1 + a, Type.Insertion);
            deleteInsertState.parent.put(deleteInsertState, e2);
            deleteInsertState.parent.put(state, getPseudoCounts(totalSequenceCount - returnValsNext.gg, returnValsNext.gl));
            state.children.put(deleteInsertState,
                    getPseudoCounts(totalSequenceCount - returnValsNext.gg, returnValsNext.gl));//delete later
            setEmission(deleteInsertState, true);
            parentListTmp.add(deleteInsertState);
            AllStateList.add(deleteInsertState);
        }
        States endState = new States(column_values_main.keySet().size(), Type.End);
//        States tmp = null;
//        double priorProb = 0.0;
        for (States s : parentListTmp) {
            switch (s.type) {
                case Insertion:
//                    priorProb += 1 - e2;
                    s.children.put(endState, 1 - e2);
                    endState.parent.put(s, (1 - e2));//plz be sure about it,right now its original log -2 val
                    break;
                case Main:
                    double tmp = 0.0;
                    for (Map.Entry<States, Double> entrySet : s.children.entrySet()) {
                        tmp += entrySet.getValue();

                    }
                    s.children.put(endState, 1 - tmp);
                    endState.parent.put(s, 1 - tmp);
                    break;
                case Deletion:

                    s.children.put(endState, 1 - e2);
                    endState.parent.put(s, 1 - e2);
                    break;
                default:
            }
            //System.out.println("name " + s.type + s.state_no);
        }
        AllStateList.add(endState);
        return AllStateList;
    }

    private void setEmission(States state, boolean isInsert) {
//emission part for insert states
        for (Map.Entry<String, Double> entrySet1 : new Config().amino_acid_maps.entrySet()) {
//                AcidCount search = b.stream().filter(ac -> ac.name.equals(entrySet1.getKey())).findFirst().orElse(null);
            if (isInsert) {
                state.emission.put(entrySet1.getKey(), getLogOdds((double) entrySet1.getValue() / 100));
            } else {

            }
        }
    }

    private static ReturnVals getLG_LL_GL(Map<String, String> multiple_sequence, final int a, final int b) {
        final ReturnVals retVals = new ReturnVals();
        if (multiple_sequence.entrySet().iterator().next().getValue().length() == b || b < 0 || a < 0) {
            retVals.ll = multiple_sequence.size();
            return retVals;
        } else {
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
                //String key = entrySet.getKey();
                String sequence = entrySet.getValue();
                try {
                    if (sequence.toCharArray()[a] == '-' && sequence.toCharArray()[b] != '-') {
                        retVals.gl++;
                    } else if (sequence.toCharArray()[a] != '-' && sequence.toCharArray()[b] == '-') {
                        retVals.lg++;
                    } else if ((sequence.toCharArray()[a] != '-' && sequence.toCharArray()[b] != '-')) {
                        retVals.ll++;
                    } else if (sequence.toCharArray()[a] == '-' && sequence.toCharArray()[b] == '-') {
                        retVals.gg++;
                    }
                } catch (ArrayIndexOutOfBoundsException ex) {
                    System.out.println("");
                }
            }
        }
        return retVals;
    }

    /**
     * this is for transition matrix
     *
     * @param sum
     * @param count
     * @return
     */
    private double getPseudoCounts(int sum, int count) {
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
    private double getPseudoCounts(int sum, int count, double acid_background_prob) {
        return getLogOdds(((double) ((double) (acid_background_prob / 100) * p + count)) / (p + sum));
//        return (double) count / sum;
    }

    private double getLogOdds(double p) {
//        return Math.log((p) / (1 - p));
        return Math.log(p);
    }

//    private double getLogits(int probability) {
//        return Math.log(probability);
////        return Math.log(probability / (1 - probability));
//    }
//
//    private double getProbability(double logitValue) {
//        double tmp = Math.exp(logitValue);
//        return (tmp) / (tmp + 1);
//
//    }
    public static class ReturnVals {

        public int ll = 0;
        public int lg = 0;
        public int gl = 0;
        public int gg = 0;
    }

    public static class AcidCount {

        public String name;
        public int count = 0;
    }
}
