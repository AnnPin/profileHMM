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
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author shad942
 */
public class Config {

    final static int p = 20;
    final static double le2 = 0.13533528323;
    final static double e2 = -2;

    Map<String, Double> amino_acid_maps;

    public static void printCSVFile(String fileName, String str) {
        PrintWriter writer;
        try {
            writer = new PrintWriter(fileName + ".csv", "UTF-8");
            writer.print(str);
            writer.close();
        } catch (FileNotFoundException | UnsupportedEncodingException ex) {
            Logger.getLogger(ProfileHmm.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public Map<String, String> readFasta(String FileName) throws IOException {
        File fl = new File(FileName);
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
        return multiple_sequence;
    }

    public Config() {
        amino_acid_maps = new HashMap<>();
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
    }

}
