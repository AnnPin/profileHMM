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

    static List<States> allStates;
    static Map<String, String> input2;

    public static void main(String[] args) throws IOException {
        Config config = new Config();
        Map<String, String> input = config.readFasta("fasta_simple");
        input2 = config.readFasta("viterbi_input");
        String viterbiInput = "";
        for (char ch : input2.entrySet().iterator().next().getValue().toCharArray()) {
            if (ch != '-') {
                viterbiInput += ch;
            }
        }
        HMMConstruct hmmc = new HMMConstruct(input);
        allStates = hmmc.getStates();
        System.out.println("total states " + allStates.size());
        config.amino_acid_maps.put("-", 0.0);
//        Stack<Path> stack = new Stack<>();
        Map<String, States> pathmap = new HashMap<>();
        Map<String, Double> valuesMap = new HashMap<>();
        List<Path> stack = new ArrayList<>();
        States[][] path = new States[viterbiInput.length() + 1][allStates.size()];
        Double[][] values = new Double[viterbiInput.length() + 1][allStates.size()];
        for (int i = 0; i < viterbiInput.length() + 1; i++) {
            Arrays.fill(values[i], -1.0);//init with zeros
        }
        values[0][0] = 1.0;

        valuesMap.put(allStates.get(0).toString(), values[0][0]);
        int ind = 1;
        for (char ch : viterbiInput.toCharArray()) {
            values[0][ind++] = 0.0;
        }
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
                pathmap.put(allStates.get(j) + "-", allStates.get(0));
                stack.add(new Path(allStates.get(j), allStates.get(0), "-"));
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

                States currentStates = allStates.get(j);

//                if (currentStates.emission.get(ch + "") == null) {
////                        && (currentStates.type != Type.Deletion
////                        && currentStates.type != Type.End)) {
//                    values[column][j] = 0.0;
//                }
//                else {
//                if (currentStates.state_no == 5 && currentStates.type) {
//                    System.out.println("State 4");
//                }
                if ((currentStates.state_no == 1 || currentStates.state_no == 0) && (currentStates.type == Type.Main)) {
                    System.out.println("found less than 0");
                }
                double emissionProbability = (currentStates.type == Type.End) ? 1
                        : ((currentStates.emission.get(ch + "") != null)
                                ? getInverseLog(currentStates.emission.get(ch + "")) : 0);
                for (Map.Entry<States, Double> p : currentStates.parent.entrySet()) {

                    if (currentStates.type == Type.Deletion) {
                        emissionProbability = (p.getKey().emission.get(ch + "") != null)
                                ? getInverseLog(p.getKey().emission.get(ch + "")) : 0;
                    }
                    if (currentStates.state_no > column) {
                        emissionProbability = 0;
                    }
                    double tmp = ((valuesMap.get((p.getKey().toString())) != null) ? (valuesMap.get((p.getKey().toString()))) : 0)
                            * (getInverseLog(p.getValue()) * emissionProbability);
                    if (probability < tmp) {
                        probability = tmp;// * getInverseLog(p.getValue()) * emissionProbability;
                        path[column][j] = p.getKey();
                        pathmap.put(currentStates.toString() + ch, p.getKey());

                        Stack<Path> tmpStack = new Stack<>();
                        for (Path tmppath : stack) {

                            if (tmppath.current != currentStates) {
                                tmpStack.add(tmppath);
                            }
                        }
                        stack = tmpStack;
                        stack.add(new Path(currentStates, p.getKey(), ch + ""));

//                        if (probability > 1) {
//                            System.out.println("pro " + probability);
//                        }
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
        //Config.printCSVFile("viterbiValues", valuesViterbi);
        getPath(pathmap);
//        getPath(stack);
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

    private static void getPath(Map<String, States> list) {
        String input = input2.values().iterator().next();
        String query = "", ch = input.charAt(input.length() - 1) + "";
        int found = input.length() - 1;
        States tmp = allStates.get((allStates.size() - 1));
        while (tmp != null) {
            System.out.println(tmp);
            tmp = list.get(tmp.toString() + ch);
//            if (found == 0) {
//                break;
//            }
            --found;
            if (found < 0) {
                ch = "-";
            } else {
                ch = input.charAt(found) + "";
            }
            while (ch == "-") {
                ch = input.charAt(--found) + "";
                if (found == 0) {
                    break;
                }
            }
        }
    }

    private static void getPath(List<Path> list) {
        String query = "";
//        , ch = "";
        String sequence_to_align = input2.values().iterator().next();
        int j = sequence_to_align.length() - 1;
//        int i = list.size() - 1;
        for (int i = list.size() - 1; i >= 0; i--) {
//        while(true){

            Path p = list.get(i);

//            System.out.println("General " + p.current.toString() + " parent " + p.parent.toString());
            if (query.equals("")) {
                //ch = p.charac;
                query = p.parent.toString();
                System.out.println("State " + p.current.toString());
            } else if (p.current.toString().equals(query)) {//(sequence_to_align.toCharArray()[j]+"").equals(p.charac) &&
//                ch = p.charac;
                j--;
                query = p.parent.toString();
                System.out.println("State " + p.current.toString());
            }
//            else if (p.current.toString().equals(query)) {
//                System.out.println(p.current);
//            }
//            System.out.println(i);
        }
    }

    private static double getInverseLog(double probability) {
        return Math.pow(Math.E, probability);
//        return probability;
    }

    public static class Path {

        States current;
        States parent;
        String charac;

        public Path(States current, States parent, String character) {
            this.current = current;
            this.parent = parent;
            this.charac = character;
        }

    }
}
