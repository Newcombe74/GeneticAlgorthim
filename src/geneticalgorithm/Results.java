/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geneticalgorithm;
import java.util.ArrayList;

/**
 *
 * @author c2-newcombe
 */
public class Results {
    Results[] results;
    ArrayList<Double> bestFitness = new ArrayList<>();
    ArrayList<Double> avgFitness = new ArrayList<>();
    ArrayList<Double> sumFitness = new ArrayList<>();

    public Results(int r) {
        this.results = new Results[r];
        
        for(int i = 0; i < r; i++){
            this.results[i] = new Results(0);
        }
    }

    public Results[] getResults() {
        return results;
    }

    public void setResults(Results[] results) {
        this.results = results;
    }
    
    public ArrayList<Double> getBestFitness() {
        return bestFitness;
    }

    public void setBestFitness(ArrayList<Double> bestFitness) {
        this.bestFitness = bestFitness;
    }
    
    public void addBestFitness(Double d){
        this.bestFitness.add(d);
    }
    
    public double getAvgBestFitness() {
        return calcAvg(this.bestFitness);
    }

    public ArrayList<Double> getAvgFitness() {
        return avgFitness;
    }

    public void setAvgFitness(ArrayList<Double> avgFitness) {
        this.avgFitness = avgFitness;
    }

    public void addAvgFitness(Double d){
        this.avgFitness.add(d);
    }
    
    public double getAvgAvgFitness() {
        return calcAvg(this.avgFitness);
    }
    
    public ArrayList<Double> getSumFitness() {
        return sumFitness;
    }

    public void setSumFitness(ArrayList<Double> sumFitness) {
        this.sumFitness = sumFitness;
    }
    
    public void addSumFitness(Double d){
        this.sumFitness.add(d);
    }
    
    public double getAvgSumFitness() {
        return calcAvg(this.sumFitness);
    }

    private static double calcAvg(ArrayList<Double> list) {
        if (list.isEmpty()) {
            return 0;
        }

        double ret = 0;

        for (Double d : list) {
            ret += d;
        }
        return ret / list.size();
    }
}
