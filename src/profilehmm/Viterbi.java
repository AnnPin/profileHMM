/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package profilehmm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import profilehmm.ProfileHmm.Type;

/**
 *
 * @author shad942
 */
public class Viterbi {

    public static void main(String[] args) throws IOException {
        Config config = new Config();
        Map<String, String> input = config.readFasta("fasta_simple");
        Map<String, String> input2 = config.readFasta("viterbi_input");
        String viterbiInput = "";
        for (char ch : input2.entrySet().iterator().next().getValue().toCharArray()) {
            if (ch != '-') {
                viterbiInput += ch;
            }
        }
        HMMConstruct hmmc = new HMMConstruct(input);
        List<States> allStates = hmmc.getStates();
        System.out.println("total states " + allStates.size());
        config.amino_acid_maps.put("-", 0.0);
        Stack<Path> stack = new Stack<>();
        States[][] path = new States[viterbiInput.length() + 1][allStates.size()];
        Double[][] values = new Double[viterbiInput.length() + 1][allStates.size()];
        for (int i = 0; i < viterbiInput.length() + 1; i++) {
            Arrays.fill(values[i], -1.0);//init with zeros
        }
        values[0][0] = 1.0;
        Map<String, Double> valuesMap = new HashMap<>();
        valuesMap.put(allStates.get(0).toString(), values[0][0]);
        String valuesViterbi = ",";
        for (int j = 0; j < allStates.size(); j++) {
            valuesViterbi += allStates.get(j).toString() + ",";
        }
        valuesViterbi += "\r\n" + "-,1.0,";
        Map<Integer, Double> deleteProbabilities = new HashMap<>();
        for (int j = 1; j < allStates.size(); j++) {

            if (allStates.get(j).type == Type.Deletion) {//as deletion are silent states

                double priorProbability = 0.0;
//                States parentDelete = (allStates.get(j).parent.keySet().stream().filter(a -> a.type == Type.Deletion || a.type == Type.Start).findFirst().orElse(null));
                States parentDelete = null;
                for (Map.Entry<States, Double> entrySet : allStates.get(j).parent.entrySet()) {
                    if (entrySet.getKey().type == Type.Deletion) {
                        parentDelete = entrySet.getKey();
                        priorProbability = getInverseLog(entrySet.getValue());
                        break;
                    } else if (entrySet.getKey().type == Type.Start) {
                        priorProbability = getInverseLog(entrySet.getValue());
                        break;
                    }
                }
                values[0][j] = priorProbability * ((parentDelete == null) ? 1 : deleteProbabilities.get(parentDelete.state_no));
                deleteProbabilities.put(allStates.get(j).state_no, values[0][j]);
            } else {
                values[0][j] = 0.0;
            }
            valuesViterbi += values[0][j] + ",";
            valuesMap.put(allStates.get(j).toString(), values[0][j]);
            if (values[0][j] > 0) {
                stack.push(new Path(allStates.get(j), allStates.get(0)));
            }
        }
        valuesViterbi += "\r\n";
        int column = 1;
//        String csv = "";
//        Arrays.asList(values[0]).forEach(n -> System.out.println(n));

        /**
         * main block
         */
        for (char ch : viterbiInput.toCharArray()) {
            valuesViterbi += ch + ",";
            for (int j = 0; j < allStates.size(); j++) {

                double probability = (values[column][j] != null) ? values[column][j] : 0.0;
                if (probability > 0) {
                    System.out.println("found less than 0");
                }
                States currentStates = allStates.get(j);
//                if (currentStates.emission.get(ch + "") == null) {
////                        && (currentStates.type != Type.Deletion
////                        && currentStates.type != Type.End)) {
//                    values[column][j] = 0.0;
//                }
//                else {
                double emissionProbability = (currentStates.type == Type.Deletion || currentStates.type == Type.End) ? 1
                        : ((currentStates.emission.get(ch + "") != null) ? getInverseLog(currentStates.emission.get(ch + "")) : 0);
                for (Map.Entry<States, Double> p : currentStates.parent.entrySet()) {
                    double tmp = ((valuesMap.get((p.getKey().toString())) != null) ? (valuesMap.get((p.getKey().toString()))) : 0);
                    if (probability
                            < (tmp * getInverseLog(p.getValue()) * emissionProbability)) {
                        probability = tmp * getInverseLog(p.getValue()) * emissionProbability;
                        path[column][j] = p.getKey();
                        stack.push(new Path(currentStates, p.getKey()));
                        if (probability > 1) {
                            System.out.println("pro " + probability);
                        }
                    }
//                    }
                    values[column][j] = probability;
                }
                if (path[column][j] != null) {
                    valuesViterbi += probability + path[column][j].toString() + ",";
                } else {
                    valuesViterbi += ",";
                }
                if (probability > 0) {
                    valuesMap.put(currentStates.toString(), probability);
                }
            }
            valuesViterbi += "\r\n";
//            allStates.forEach(n -> {
//                System.out.print(n.toString() + ",");
//            });
//            System.out.println("");
//            Arrays.asList(values[column]).forEach(n -> {
//                System.out.print(n + ",");
//            });
//            System.out.println("");
            column++;
        }
        Config.printCSVFile("viterbiValues.csv", valuesViterbi);
        getPath(stack);

//        System.out.println(allStates.get(allStates.size()-1).state_no+""+
//                allStates.get(allStates.size()-1).type);
//        for (int i = 0; i < allStates.size(); i++) {
////            System.out.println(allStates.get(i).type+""+allStates.get(i).state_no+"__"+values[viterbiInput.length()][i]);
//            System.out.println((path[viterbiInput.length()][i] != null) ? path[viterbiInput.length()][i].toString() : 0
//                    + "__" + values[viterbiInput.length()][i]);
//
//        }
//        Arrays.asList(values[viterbiInput.length()]).forEach(n -> System.out.println(n));
    }

    private static void getPath(Stack<Path> list) {
        String query = "";
        while (!list.empty()) {
            Path p = list.pop();
            if (query.equals("")) {
                query = p.parent.toString();
                System.out.println("State " + p.current.toString());
            } else if (p.current.toString().equals(query)) {
                query = p.parent.toString();
                System.out.println("State " + p.current.toString());
            }

        }
    }

    private static double getInverseLog(double probability) {
        return Math.pow(Math.E, probability);
//        return probability;
    }

    public static class Path {

        States current;
        States parent;

        public Path(States current, States parent) {
            this.current = current;
            this.parent = parent;
        }

    }
}
