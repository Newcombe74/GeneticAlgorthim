/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package geneticalgorithm;

/**
 *
 * @author c2-newcombe
 */
public class Individual {
    private int[] chromosome;
    private int fitness;

    public Individual(int[] chromosome) {
        this.chromosome = chromosome;
    }
    
    public Individual(int[] chromosome, int fitness) {
        this.chromosome = chromosome;
        this.fitness = fitness;
    }
    
    public int[] getChromosome() {
        return chromosome;
    }

    public void setChromosome(int[] gene) {
        this.chromosome = gene;
    }

    public int getFitness() {
        return fitness;
    }

    public void setFitness(int fitness) {
        this.fitness = fitness;
    }
    
    public void incFitness(){
        this.fitness++;
    }
}
