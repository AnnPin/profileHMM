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
    static String valuesViterbi = ",";
    static Map<String, Integer> mapIndexStates;
    static Map<String, States> pathmap = new HashMap<>();
    static Map<String, String> pathmapSecond = new HashMap<>();

    public static void main(String[] args) throws IOException {
        Config config = new Config();
        Map<String, String> input = config.readFasta("fasta_simple");
        input2 = config.readFasta("viterbi_input_1");
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

        Map<String, Double> valuesMap = new HashMap<>();
        mapIndexStates = new HashMap<>();
        int ind = 0;
        for (States s : allStates) {
            mapIndexStates.put(s.toString(), ind++);
        }
        List<Path> stack = new ArrayList<>();
        States[][] path = new States[viterbiInput.length() + 1][allStates.size()];
        Double[][] values = new Double[viterbiInput.length() + 1][allStates.size()];
        for (int i = 0; i < viterbiInput.length() + 1; i++) {
            Arrays.fill(values[i], -1.0);//init with zeros
        }
        values[0][0] = 1.0;

        valuesMap.put(allStates.get(0).toString(), values[0][0]);
        ind = 1;
        for (char ch : viterbiInput.toCharArray()) {
            values[0][ind++] = 0.0;
        }

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
                pathmapSecond.put("0," + j, "0,0");
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
        getRestOfTheThingsDone(viterbiInput, values);
//        for (char ch : viterbiInput.toCharArray()) {
//            valuesViterbi += ch + ",";
//            for (int j = 0; j < allStates.size(); j++) {
//
//                double probability = (values[column][j] != null) ? values[column][j] : 0.0;
//
//                States currentStates = allStates.get(j);
//                if (currentStates.state_no > column + 1) {
//                    continue;
//                }
////                if (currentStates.emission.get(ch + "") == null) {
//////                        && (currentStates.type != Type.Deletion
//////                        && currentStates.type != Type.End)) {
////                    values[column][j] = 0.0;
////                }
////                else {
////                if (currentStates.state_no == 5 && currentStates.type) {
////                    System.out.println("State 4");
////                }
//                if (currentStates.type == Type.End && column == viterbiInput.length()) {
//                    System.out.println("");
//                }
//                double emissionProbability = (currentStates.type == Type.End) ? 1
//                        : ((currentStates.emission.get(ch + "") != null)
//                                ? getInverseLog(currentStates.emission.get(ch + "")) : 0);
//                for (Map.Entry<States, Double> p : currentStates.parent.entrySet()) {
//
//                    if (currentStates.type == Type.Deletion) {
//                        emissionProbability = (p.getKey().emission.get(ch + "") != null)
//                                ? getInverseLog(p.getKey().emission.get(ch + "")) : 0;
//                    }
//                    if (currentStates.state_no > column && currentStates.type != Type.End) {
//                        emissionProbability = 0;
//                    }
//                    String characFind = "" + ((column == 1) ? "" : viterbiInput.toCharArray()[column - 1]);
//                    double prevValue = ((valuesMap.get((p.getKey().toString())) != null) ? (valuesMap.get((p.getKey().toString()))) : 0);
//                    double tmp
//                            = prevValue * (getInverseLog(p.getValue()) * emissionProbability);
//                    if (probability < tmp) {
//                        probability = tmp;// * getInverseLog(p.getValue()) * emissionProbability;
//                        path[column][j] = p.getKey();
//                        pathmap.put(currentStates.toString() + ch, p.getKey());
//
//                        Stack<Path> tmpStack = new Stack<>();
//                        for (Path tmppath : stack) {
//
//                            if (tmppath.current != currentStates) {
//                                tmpStack.add(tmppath);
//                            }
//                        }
//                        stack = tmpStack;
//                        stack.add(new Path(currentStates, p.getKey(), ch + ""));
//
////                        if (probability > 1) {
////                            System.out.println("pro " + probability);
////                        }
//                    }
////                    }
//                    values[column][j] = probability;
//                }
//                if (probability > 0) {
//                    valuesViterbi += probability + path[column][j].toString() + ",";
//                } else {
//                    valuesViterbi += ",";
//                }
//                if (probability > 0) {
//                    valuesMapTmp.put(currentStates.toString() + viterbiInput.toCharArray()[column - 1], probability);
//                }
//            }
//            for (Map.Entry<String, Double> entrySet : valuesMapTmp.entrySet()) {
//                valuesMap.put(entrySet.getKey(), entrySet.getValue());
//            }
//            valuesMapTmp = new HashMap<>();
//            valuesViterbi += "\r\n";
////            allStates.forEach(n -> {
////                System.out.print(n.toString() + ",");
////            });
////            System.out.println("");
////            Arrays.asList(values[column]).forEach(n -> {
////                System.out.print(n + ",");
////            });
////            System.out.println("");
//            column++;
//        }
        Config.printCSVFile("viterbiValues", valuesViterbi);
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

        input = input.replace("-", "");
        String query = "", ch = input.charAt(input.length() - 1) + "";
        int found = input.length() - 1;
//        States tmp = allStates.get((allStates.size() - 1));
//        while (tmp != null) {
//            System.out.println(tmp);
//            tmp = list.get(tmp.toString() + ch);
////            if (found == 0) {
////                break;
////            }
//            --found;
//            if (found < 0) {
//                ch = "-";
//            } else {
//                ch = input.charAt(found) + "";
//            }
//            
//        }
//        System.out.println(allStates.get(allStates.size() - 1));
        States lastState = null;
        query = (allStates.get(allStates.size() - 1)).toString() + input.toCharArray()[input.length() - 1] + "";
        for (int i = input.length() - 1; i > -1; i--) {
            for (Map.Entry<String, States> entrySet : list.entrySet()) {
                if (entrySet.getKey().equals(query)) {
                    System.out.println(entrySet.getKey());
                    String add = "";
                    add = ((i < 0) ? "-" : input.toCharArray()[i] + "");

//                    if (i != input.length() - 2) {
//                        add = ((i < 0) ? "-" : input.toCharArray()[i] + "");
//                    } else {
//                        add = input.toCharArray()[input.length()-1] + "";
//                    }
                    lastState = entrySet.getValue();
                    query = entrySet.getValue().toString() + add;
                    break;
                }

            }

        }

//        int i = input.length() - 1;
//        System.out.println(query);
        while (lastState != null) {
            boolean found2 = false;
            for (Map.Entry<String, States> entrySet : list.entrySet()) {
                if (entrySet.getKey().equals(query) || entrySet.getKey().contains(query)) {
                    System.out.println(query);
                    query = entrySet.getValue().toString();
                    lastState = list.get(query);
                    found2 = true;
                }
            }
            if (!found2) {
                lastState = null;
            }
        }
        System.out.println(allStates.get(0));
//        while (i > 0) {
//            if (list.containsKey(query)) {
//                System.out.println(list.get(query));
//                query = list.get(query).toString() + input.toCharArray()[--i];
//            }
//        }
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

    private static void getRestOfTheThingsDone(String viterbiInput, Double[][] values) {
        double probabilityMain = 0.0;
        boolean changed = false;
        for (int i = 0; i < viterbiInput.length(); i++) {

            for (int j = 0; j < allStates.size(); j++) {
                values[i + 1][j] = 0.0;
            }
            String lookingFor = viterbiInput.toCharArray()[i] + "";
            valuesViterbi += lookingFor + ",";
            for (int j = 0; j < allStates.size(); j++) {
                States currentState = allStates.get(j);
                if (currentState.type == Type.Deletion && currentState.state_no == 3) {
                    System.out.print("");
                }

                if (currentState.type != Type.End) {
//                    if (currentState.state_no != i && (currentState.type != Type.Insertion)) {//big logic
//                        continue;
//                    }
//                    } else if (currentState.state_no > i + 1 && (currentState.type == Type.Insertion)) {
//                        continue;
//                    }
                }
//                if (currentState.state_no - 1 != i && currentState.type == Type.End) {
//                    continue;
//                }
                double emissionProbability = 0.0;
                if (currentState.type == Type.Main && currentState.state_no == 0) {
                    System.out.print("");
                }
                States parentUpper = null;
                emissionProbability = (currentState.emission.get(lookingFor) != null) ? currentState.emission.get(lookingFor) : 0;
                for (Map.Entry<States, Double> entrySet : currentState.parent.entrySet()) {
                    States parent = entrySet.getKey();
                    if (parent.state_no == 0 && parent.type == Type.Main) {
                        System.out.println("");
                    }
                    int parentIndex = mapIndexStates.get(parent.toString());
                    int currentIndex = mapIndexStates.get(currentState.toString());

                    double prevProbability = (parent == currentState) ? values[i][currentIndex] : values[i][parentIndex];
                    double transitionProbability = entrySet.getValue();
                    if (currentState.type == Type.Deletion) {
//                        State lastParent
//                        while(emissionProbability==0){
//                        
//                        }
                        emissionProbability = (parent.emission.get(lookingFor) == null) ? 0 : parent.emission.get(lookingFor);
                    } else if (currentState.type == Type.End) {
                        prevProbability = values[i + 1][parentIndex];
                        emissionProbability = (parent.emission.get(lookingFor) == null) ? 0 : parent.emission.get(lookingFor);
                        if (emissionProbability == 0) {
                            emissionProbability = 1;
                        }
                    } else if (currentState.state_no < i && currentState.type==Type.Main) {
                        emissionProbability = 0;

                    }
//                    else if (currentState.type == Type.End) {
//                        emissionProbability = 1;
//                    }

                    probabilityMain = ((emissionProbability == 0) ? 0 : getInverseLog(emissionProbability)) * prevProbability * getInverseLog(transitionProbability);
                    if (currentState.state_no == 2 && currentState.type == Type.Main) {
                        System.out.println("");
                    }
                    if (values[i + 1][currentIndex] < probabilityMain) {
                        if(i+1==1 && currentIndex==2)
                            System.out.println("");
                        values[i + 1][currentIndex] = probabilityMain;
                        pathmap.put(currentState.toString() + lookingFor, parent);
                        pathmapSecond.put(((currentState.type == Type.End) ? i : (i + 1)) + "," + currentIndex, i + "," + currentIndex);
                        changed = true;
                        parentUpper = parent;
                    }

                }
                if (probabilityMain > 0 && changed) {
                    valuesViterbi += probabilityMain + parentUpper.toString() + ",";
                } else {
                    valuesViterbi += "0,";
                }
                changed = false;
            }
            valuesViterbi += "\r\n";
        }
        System.out.println(pathmap.get(viterbiInput.length() + "," + allStates.size()));
        return;

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
