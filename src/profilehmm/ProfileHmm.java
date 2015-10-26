/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package profilehmm;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

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
//    static int p = 20;
//    final static double e2 = 0.13533528323;
    public static void main(String[] args) throws IOException {
        if (args.length < 3) {
            System.out.println("Please specify the input in arguments(input fasta, emission filename, transition filename)");
            return;
        }

        Map<String, String> multiple_sequence = new Config().readFasta(args[0]);

//        String[][] sequences = new String[multiple_sequence.values().iterator().next().length()][multiple_sequence.keySet().size()];
//        multiple_sequence.keySet().stream().forEach((name) -> {
//
//        });
        HMMConstruct hmmc = new HMMConstruct(multiple_sequence);
        List<States> allStates = hmmc.getStates();
        System.out.println("total states " + allStates.size());
        Config.printCSVFile(args[1], printEmissionMatrix(allStates));

        Config.printCSVFile(args[2], printTransmissionMatrix(allStates));
    }

    

    private static String printEmissionMatrix(List<States> states) {
        String csv = ",";
        for (Map.Entry<String, Double> entrySet : new Config().amino_acid_maps.entrySet()) {
            csv += entrySet.getKey() + ",";
        }
        csv += "\r\n";
//        List<String> printedStates = new ArrayList<>();
        for (States s : states) {
            if (s.type == Type.Main || s.type == Type.Insertion) {
                if (s.emission.size() == new Config().amino_acid_maps.size()) {
//                printedStates.add(s.state_no + "" + s.type);
                    SortedSet<String> keys = new TreeSet<String>(s.emission.keySet());
                    csv += (s.type + "-" + s.state_no + ",");
                    for (String acid : keys) {
                        csv += (s.emission.get(acid) + ",");
                    }
                    csv += "\r\n";
                }
            }
        }
        return csv;
    }

    private static String printTransmissionMatrix(List<States> states) {
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
        return "States," + tmp + "\r\n" + csv;
    }

    private static void printAll(States state) {
        if (state.type != Type.Main && state.type != Type.Start) {
            return;
        }
        double i = 0;
        for (Map.Entry<States, Double> entrySet : state.children.entrySet()) {
            //System.out.println("Key " + entrySet.getKey().type + "" + entrySet.getKey().state_no + " val " + entrySet.getValue());
            i += entrySet.getValue();
        }
//        System.out.println(state.type + "" + state.state_no + " Total = " + i);
        if (i > 1.0 || i < 1.0) {
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

}
