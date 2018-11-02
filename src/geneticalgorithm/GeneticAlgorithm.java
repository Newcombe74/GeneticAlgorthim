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
public class GeneticAlgorithm {

    //Selection Type
    public final int SELECTION_TOURNEMENT = 1,
            SELECTION_ROULETTE = 2;
    
    //Hyperparameters
    private int populationSize,
            numberOfGenerations,
            chromosomeSize;
    private double probabilityOfMutation;

    //Population
    private Individual population[];
    private Individual offspring[];
    
    //Results
    private Results genResults;

    public GeneticAlgorithm() {
    }

    public GeneticAlgorithm(int populationSize, int numberOfGenerations, int chromosomeSize) {
        this.populationSize = populationSize;
        this.numberOfGenerations = numberOfGenerations;
        this.chromosomeSize = chromosomeSize;
        this.probabilityOfMutation = calcRandMutationRate();
    }
    
    public GeneticAlgorithm(int populationSize, int numberOfGenerations, int chromosomeSize, double probabilityOfMutation) {
        this.populationSize = populationSize;
        this.numberOfGenerations = numberOfGenerations;
        this.chromosomeSize = chromosomeSize;
        this.probabilityOfMutation = probabilityOfMutation;
    }
    
    public void run(int selectionType) {
        //SET each individuals genes to be 1 or 0 at random
        for (int i = 0; i < population.length; i++) {
            int[] genes = new int[chromosomeSize];

            for (int j = 0; j < genes.length; j++) {
                genes[j] = (int) ((Math.random() * 2) % 2);
            }
            population[i] = new Individual(genes);
        }

        for (int g = 0; g < numberOfGenerations; g++) {
            population = calcFitness(population);

            genResults.results[g].addBestFitness(new Double(bestFitness(population)));
            genResults.results[g].addAvgFitness(new Double(avgFitness(population)));
            genResults.results[g].addSumFitness(new Double(sumFitness(population)));

            offspring = crossover();
            
            offspring = mutate();

            offspring = calcFitness(offspring);

            population = selection(selectionType);

        }
    }
    
    //START_Selection
    private Individual[] selection(int selectionType) {
        switch (selectionType) {
            case SELECTION_TOURNEMENT:
                return tournementSelection();
            case SELECTION_ROULETTE:
                return rouletteWheelSelection();
            default:
                System.err.println("Selection type not found: " + selectionType);
                return null;
        }
    }

    private Individual[] tournementSelection() {
        Individual[] nextGen = new Individual[populationSize];

        if (offspring.length > 0) {
            for (int i = 0; i < populationSize; i++) {
                int parent = (int) ((Math.random() * population.length) % population.length);
                int child = (int) ((Math.random() * offspring.length) % offspring.length);

                if (population[parent].getFitness() >= offspring[child].getFitness()) {
                    nextGen[i] = new Individual(population[parent].getChromosome());
                } else {
                    nextGen[i] = new Individual(offspring[child].getChromosome());
                }
            }

            return nextGen;
        } else {
            return population;
        }
    }

    private Individual[] rouletteWheelSelection() {
        Individual[] nextGen = new Individual[populationSize];

        //Put parents and children into a single population
        Individual[] currentGen = new Individual[population.length + offspring.length];
        int n = 0;
        for (Individual inividual : population) {
            currentGen[n] = inividual;
            n++;
        }
        for (Individual inividual : offspring) {
            currentGen[n] = inividual;
            n++;
        }

        int totalFitness = sumFitness(currentGen);
        for (int i = 0; i < populationSize; i++) {
            int runningTotal = 0, j = 0;

            int selectionPoint = (int) ((Math.random() * totalFitness) % totalFitness);

            while (runningTotal <= selectionPoint) {
                runningTotal += currentGen[j].getFitness();
                j++;
            }
            nextGen[i] = currentGen[j - 1];
        }

        return nextGen;
    }
    //END_Selection

    //START_Crossover
    private Individual[] crossover() {
        ArrayList<Individual> children = new ArrayList<>();

        for (int i = 0; i < populationSize - 1; i++) {
                children.addAll(singlePointCrossover(
                        population[i].getChromosome(),
                        population[(i + 1)].getChromosome()));
        }

        Individual[] ret = new Individual[children.size()];

        for (int i = 0; i < ret.length; i++) {
            ret[i] = children.get(i);
        }
        return ret;
    }

    private ArrayList<Individual> singlePointCrossover(int[] parent1, int[] parent2) {
        ArrayList<Individual> children = new ArrayList<>();
        int[][] crossoverGenes = new int[2][chromosomeSize];
        int child1 = 0, child2 = 1;
        int crossoverPoint = (int) (Math.random() * chromosomeSize) - 1;

        for (int i = 0; i < chromosomeSize; i++) {
            if (i < crossoverPoint) {
                crossoverGenes[child1][i] = parent1[i];
            } else {
                crossoverGenes[child1][i] = parent2[i];
            }
            if (i < crossoverPoint) {
                crossoverGenes[child2][i] = parent2[i];
            } else {
                crossoverGenes[child2][i] = parent1[i];
            }
        }

        children.add(new Individual(crossoverGenes[child1]));
        children.add(new Individual(crossoverGenes[child2]));

        return children;
    }
    //END_Crossover

    //START_Mutation
    private Individual[] mutate() {
        for (Individual child : offspring) {
            child.setChromosome(mutateChromosome(child.getChromosome()));
        }
        return offspring;
    }

    private int calcRandMutationRate() {
        int min = 1 / populationSize;
        int max = 1 / chromosomeSize;

        if (min >= max) {
            throw new IllegalArgumentException(
                    "Population size must be greater than chromosome size");
        }
        return (int) (Math.random() * ((max - min) + 1)) + min;
    }

    private int[] mutateChromosome(int[] chrom) {
        int[] mutatedGenes = chrom;

        for (int i = 0; i < chromosomeSize; i++) {
            double m = Math.random();
            if (probabilityOfMutation >= m) {
                mutatedGenes[i] = invert(chrom[i]);
            }
        }

        return mutatedGenes;
    }

    private int invert(int gene) {
        if (gene == 1) {
            return 0;
        } else {
            return 1;
        }
    }
    //END_Mutation

    //START_Fitness
    private Individual[] calcFitness(Individual[] pop) {
        return null;
    }
    
    private int avgFitness(Individual[] pop) {
        return sumFitness(pop) / pop.length;
    }

    private int bestFitness(Individual[] pop) {
        int ret = 0;
        for (Individual i : pop) {
            if (i.getFitness() > ret) {
                ret = i.getFitness();
            }
        }
        return ret;
    }

    private int sumFitness(Individual[] pop) {
        int ret = 0;
        for (Individual i : pop) {
            ret += i.getFitness();
        }
        return ret;
    }
    //END_Fitness

    //START_Utils
    private double calcAvg(double[] arr) {
        if (arr.length == 0) {
            return 0;
        }

        double ret = 0;

        for (int i = 0; i < arr.length; i++) {
            ret += arr[i];
        }
        return ret / arr.length;
    }

    private double calcPerc(double a, double b) {
        return (100 / b) * a;
    }
    //END_Utils
}
