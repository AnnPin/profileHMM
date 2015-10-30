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
import profilehmm.ProfileHmm.Type;

/**
 *
 * @author Gunner
 */
public class States {

    int state_no;
    Type type;
    Map<States, Double> parent;
    Map<States, Double> children;//no need, just checking
    Map<String, Double> emission;
    boolean silent = false;

    public Map<String, Double> getEmission() {
        return emission;
    }

    public void setEmission(Map<String, Double> emission) {
        this.emission = emission;
    }

    public States(int state_no, Type type) {
        this.state_no = state_no;
        this.type = type;
        parent = new HashMap<>();
        children = new HashMap<>();
        emission = new HashMap<String, Double>();
        if (type == Type.Deletion || type == Type.End) {
            silent = true;
        }
    }

    public int getState_no() {
        return state_no;
    }

    public void setState_no(int state_no) {
        this.state_no = state_no;
    }

    public Map<States, Double> getParent() {
        return parent;
    }

    public void setParent(HashMap<States, Double> parent) {
        this.parent = parent;
    }

    @Override
    public String toString() {
        return "s=" + state_no + "-t=" + type;
    }

    @Override
    public boolean equals(Object obj) {
        return ((States) obj).state_no == this.state_no && ((States) obj).type == this.type;//To change body of generated methods, choose Tools | Templates.
    }

}
