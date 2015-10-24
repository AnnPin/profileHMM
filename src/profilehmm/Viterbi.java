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
        States[][] path = new States[viterbiInput.length() + 1][allStates.size()];
        Double[][] values = new Double[viterbiInput.length() + 1][allStates.size()];
//        Arrays.fill(values, 0.0);//init with zeros
        values[0][0] = 1.0;
        Map<String, Double> valuesMap = new HashMap<>();
        valuesMap.put(allStates.get(0).state_no + "" + allStates.get(0).type, values[0][0]);

        Map<Integer, Double> deleteProbabilities = new HashMap<>();
        for (int j = 1; j < allStates.size(); j++) {
            if (allStates.get(j).type == Type.Deletion) {//as deletion are silent states

                double priorProbability = 0.0;
//                States parentDelete = (allStates.get(j).parent.keySet().stream().filter(a -> a.type == Type.Deletion || a.type == Type.Start).findFirst().orElse(null));
                States parentDelete = null;
                for (Map.Entry<States, Double> entrySet : allStates.get(j).parent.entrySet()) {
                    if (entrySet.getKey().type == Type.Deletion) {
                        parentDelete = entrySet.getKey();
                        priorProbability = entrySet.getValue();
                        break;
                    } else if (entrySet.getKey().type == Type.Start) {
                        priorProbability = entrySet.getValue();
                        break;
                    }
                }
                values[0][j] = priorProbability * ((parentDelete == null) ? 1 : deleteProbabilities.get(parentDelete.state_no));
                deleteProbabilities.put(allStates.get(j).state_no, values[0][j]);
            } else {
                values[0][j] = 0.0;
            }
            valuesMap.put(allStates.get(j).state_no + "" + allStates.get(j).type, values[0][j]);

        }
        int column = 1;

//        Arrays.asList(values[0]).forEach(n -> System.out.println(n));
        for (char ch : viterbiInput.toCharArray()) {
            for (int j = 0; j < allStates.size(); j++) {
                if (j + 1 == allStates.size() && column == viterbiInput.length()) {

                    System.out.println("test");
                }
                double probability = (values[column][j] != null) ? values[column][j] : 0.0;
                States currentStates = allStates.get(j);
                if (currentStates.emission.get(ch + "") == null) {
//                        && (currentStates.type != Type.Deletion
//                        && currentStates.type != Type.End)) {
                    values[column][j] = 0.0;
                } else {
                    double emissionProbability = getInverseLog(currentStates.emission.get(ch + ""));
                    for (Map.Entry<States, Double> p : currentStates.parent.entrySet()) {
                        double tmp = ((valuesMap.get((p.getKey().state_no + "" + p.getKey().type)) != null) ? (valuesMap.get((p.getKey().state_no + "" + p.getKey().type))) : 0);
                        if (probability < tmp
                                * p.getValue() * emissionProbability) {
                            probability = ((valuesMap.get((p.getKey().state_no + "" + p.getKey().type)) != null) ? (valuesMap.get((p.getKey().state_no + "" + p.getKey().type))) : 0)
                                    * p.getValue() * emissionProbability;
                            path[column][j] = currentStates;
                        }
                    }
                    values[column][j] = probability;
                }
                if (probability > 0) {
                    valuesMap.put(currentStates.state_no + "" + currentStates.type, probability);
                }
            }
            column++;
        }
//        System.out.println(allStates.get(allStates.size()-1).state_no+""+
//                allStates.get(allStates.size()-1).type);
        for (int i = 0; i < allStates.size(); i++) {
//            System.out.println(allStates.get(i).type+""+allStates.get(i).state_no+"__"+values[viterbiInput.length()][i]);
            System.out.println((path[viterbiInput.length()][i]!=null)?path[viterbiInput.length()][i].toString():0
                    + "__" + values[viterbiInput.length()][i]);

        }
//        Arrays.asList(values[viterbiInput.length()]).forEach(n -> System.out.println(n));
    }

    private static double getInverseLog(double probability) {
        return Math.pow(Math.E, probability);
    }
}
