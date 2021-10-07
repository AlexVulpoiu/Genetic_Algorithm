#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

ofstream g("evolution.txt");

struct Chromosome {

    string representation;
    double value{}, fitness{}, probability{}, rightBoundProbability{};
};

int step, populationDimension, precision, stages, chromosomeLength;
double leftBound, rightBound, a, b, c, crossoverProbability, mutationProbability, totalFitness;
Chromosome elitist;
vector<int> crossoverParticipants;
vector<Chromosome> chromosomes, intermediatePopulation;

void read() {

    ifstream f("input.in");

    f >> populationDimension;
    f >> leftBound >> rightBound;
    f >> a >> b >> c;
    f >> precision;
    f >> crossoverProbability;
    f >> mutationProbability;
    f >> stages;
    f.close();
}

double f(double value) {    /// f(v) = a * v^2 + b * v + c
    return a * pow(value, 2) + b * value + c;
}

string putSpaces(int value, int maxValue) {     /// put additional spaces for numbers alignment

    int digits, maxDigits;
    string spaces;

    /// number of digits for value
    if(!value) {
        digits = 1;
    } else {
        digits = 0;
        while(value) {
            digits++;
            value /= 10;
        }
    }

    /// number of digits for maxValue
    if(!maxValue) {
        maxDigits = 1;
    } else {
        maxDigits = 0;
        while(maxValue) {
            maxDigits++;
            maxValue /= 10;
        }
    }

    spaces = "";    /// a string containing additional spaces that need to be added
    for(int i = 0; i < maxDigits - digits; i++)
        spaces += ' ';

    return spaces;
}

double computeValue(double value) {     /// codified value = the value in [leftBound, rightBound] interval, represented by a chromosome
    return (1.0 * (rightBound - leftBound)) / (pow(2, chromosomeLength) - 1.0) * value + leftBound;
}

void generateChromosomes() {

    totalFitness = 0;
    chromosomeLength = ceil(log2(1.0 * (rightBound - leftBound) * pow(10, precision)));
    g << "Initial population\n";
    for(int i = 0; i < populationDimension; i++) {

        Chromosome auxChromosome;

        for(int j = 0; j < chromosomeLength; j++) {     /// generating random chromosome representation
            int bit = (int)(random() % 2);
            auxChromosome.representation += to_string(bit);
        }

        auxChromosome.value = 0.0;
        for(int j = chromosomeLength - 1; j >= 0; j--) {
            if(auxChromosome.representation[j] == '1') {    /// computing the value represented by the chromosome
                auxChromosome.value += pow(2, chromosomeLength - 1 - j);
            }
        }
        auxChromosome.value = computeValue(auxChromosome.value);
        auxChromosome.fitness = f(auxChromosome.value);  /// fitness function for chromosome value is the function applied to the value
        totalFitness += auxChromosome.fitness;
        chromosomes.push_back(auxChromosome);

        /// print chromosome values
        g << "  " << putSpaces(i + 1, populationDimension) << i + 1 << ". " << auxChromosome.representation << ",  x = ";
        if(auxChromosome.value >= 0.0) {
            g << ' ';
        }
        g << fixed << setprecision(6) << auxChromosome.value << ",  f = " << f(auxChromosome.value) << '\n';
    }
    g << '\n';
}

void selectionProbabilities() {

    if(step == 1) {
        g << "Selection probabilities\n";
    }
    for(int i = 0; i < chromosomes.size(); i++) {

        chromosomes[i].probability = chromosomes[i].fitness / totalFitness;
        if(step == 1) {
            g << "  chromosome " << putSpaces(i + 1, populationDimension) << i + 1 << ", probability: " << chromosomes[i].probability << '\n';
        }
        if(!i) {
            chromosomes[i].rightBoundProbability = chromosomes[i].probability;
        } else {
            chromosomes[i].rightBoundProbability = chromosomes[i - 1].rightBoundProbability + chromosomes[i].probability;
        }
    }

    if(step == 1) {
        g << "\nSelection probabilities intervals\n  0 ";
        for (auto & chromosome : chromosomes) {
            g << chromosome.rightBoundProbability << ' ';
        }
        g << "\n\n";
    }
}

void selectChromosomes() {

    if(step == 1) {
        g << "Selection\n";
    }
    intermediatePopulation.clear();
    for(int i = 0; i < populationDimension - 1; i++) {

        double u = (double)random() / RAND_MAX;     /// generate a random probability
        if(step == 1) {
            g << "  " << putSpaces(i + 1, populationDimension - 1) << i + 1 << ". u = " << u << ", select chromosome ";
        }

        int left, right, middle;
        left = 0;
        right = populationDimension - 1;
        while(left <= right) {      /// binary search for the first chromosome whose rightBoundProbability is greater than u
            middle = left + (right - left) / 2;

            if(chromosomes[middle].rightBoundProbability > u && (!middle || chromosomes[middle - 1].rightBoundProbability <= u)) {
                intermediatePopulation.push_back(chromosomes[middle]);
                break;
            } else if(chromosomes[middle].rightBoundProbability > u) {
                right = middle - 1;
            } else {
                left = middle + 1;
            }
        }
        if(step == 1) {
            g << putSpaces(middle + 1, populationDimension) << middle + 1 << '\n';
        }
    }

    if(step == 1) {
        g << "\nAfter selection\n";
        for (int i = 0; i < intermediatePopulation.size(); i++) {

            g << "  " << putSpaces(i + 1, intermediatePopulation.size()) << i + 1 << ". " << intermediatePopulation[i].representation << ",  x = ";
            if (intermediatePopulation[i].value > 0) {
                g << ' ';
            }
            g << intermediatePopulation[i].value << ",  f = " << intermediatePopulation[i].fitness << '\n';
        }
        g << '\n';
    }
}

void selectForCrossover() {

    if(step == 1) {
        g << "Crossover probability " << crossoverProbability << '\n';
    }
    crossoverParticipants.clear();
    for(int i = 0; i < populationDimension - 1; i++) {

        double u = (double)random() / RAND_MAX;     /// generate a random probability
        if(step == 1) {
            g << "  " << putSpaces(i + 1, populationDimension - 1) << i + 1 << ". " << intermediatePopulation[i].representation << ",  u = " << u;
        }
        if(u < crossoverProbability) {
            if(step == 1) {
                g << " < " << crossoverProbability << " participates\n";
            }
            crossoverParticipants.push_back(i);
        } else {
            if(step == 1) {
                g << '\n';
            }
        }

    }
    if(step == 1) {
        g << '\n';
    }

    shuffle(crossoverParticipants.begin(), crossoverParticipants.end(), mt19937(random_device()()));    /// shuffle the crossover participants
    if(crossoverParticipants.size() > 1) {
        for (int i = 0; i < crossoverParticipants.size(); i += 2) {  /// combine groups of 2 chromosomes; the last group may have 3 elements

            int breakPoint = (int)(random() % chromosomeLength);
            if (crossoverParticipants.size() % 2 == 1 && i == crossoverParticipants.size() - 3) { /// in this case, we have 3 chromosomes in the group

                int index1 = crossoverParticipants[i], index2 = crossoverParticipants[i + 1], index3 = crossoverParticipants[i + 2];

                if (step == 1) {
                    g << "Crossover between chromosome " << index1 + 1 << ", chromosome " << index2 + 1 << " and chromosome " << index3 + 1 << '\n';
                    g << "  " << intermediatePopulation[index1].representation << ' ' << intermediatePopulation[index2].representation << ' ' << intermediatePopulation[index3].representation << " point " << breakPoint << '\n';
                }
                if (breakPoint) {    /// if break point is 0, the chromosomes will not be modified

                    string auxRepresentation = intermediatePopulation[index1].representation;

                    /// chromosome 1
                    intermediatePopulation[index1].representation = intermediatePopulation[index1].representation.substr(0, breakPoint) + intermediatePopulation[index2].representation.substr(breakPoint);
                    intermediatePopulation[index1].value = 0.0;     /// changing the representation and computing the new value and fitness
                    for (int j = chromosomeLength - 1; j >= 0; j--) {
                        if (intermediatePopulation[index1].representation[j] == '1') {    /// computing the value represented by the chromosome
                            intermediatePopulation[index1].value += pow(2, chromosomeLength - 1 - j);
                        }
                    }
                    intermediatePopulation[index1].value = computeValue(intermediatePopulation[index1].value);
                    intermediatePopulation[index1].fitness = f(intermediatePopulation[index1].value);

                    /// chromosome 2
                    intermediatePopulation[index2].representation = intermediatePopulation[index2].representation.substr(0, breakPoint) + intermediatePopulation[index3].representation.substr(breakPoint);
                    intermediatePopulation[index2].value = 0.0;     /// changing the representation and computing the new value and fitness
                    for (int j = chromosomeLength - 1; j >= 0; j--) {
                        if (intermediatePopulation[index2].representation[j] == '1') {    /// computing the value represented by the chromosome
                            intermediatePopulation[index2].value += pow(2, chromosomeLength - 1 - j);
                        }
                    }
                    intermediatePopulation[index2].value = computeValue(intermediatePopulation[index2].value);
                    intermediatePopulation[index2].fitness = f(intermediatePopulation[index2].value);

                    /// chromosome 3
                    intermediatePopulation[index3].representation = intermediatePopulation[index3].representation.substr(0, breakPoint) + auxRepresentation.substr(breakPoint);
                    intermediatePopulation[index3].value = 0.0;     /// changing the representation and computing the new value and fitness
                    for (int j = chromosomeLength - 1; j >= 0; j--) {
                        if (intermediatePopulation[index3].representation[j] == '1') {    /// computing the value represented by the chromosome
                            intermediatePopulation[index3].value += pow(2, chromosomeLength - 1 - j);
                        }
                    }
                    intermediatePopulation[index3].value = computeValue(intermediatePopulation[index3].value);
                    intermediatePopulation[index3].fitness = f(intermediatePopulation[index3].value);

                    if (step == 1) {
                        g << "  Result  " << intermediatePopulation[index1].representation << ' ' << intermediatePopulation[index2].representation << ' ' << intermediatePopulation[index3].representation << '\n';
                    }
                }

                if (step == 1) {
                    g << "  Result  " << intermediatePopulation[index1].representation << ' ' << intermediatePopulation[index2].representation << ' ' << intermediatePopulation[index3].representation << '\n';
                }

                break;  /// the 3-group is the last one

            } else {    /// in this case, we combine 2 chromosomes

                int index1 = crossoverParticipants[i], index2 = crossoverParticipants[i + 1];

                if (step == 1) {
                    g << "Crossover between chromosome " << index1 + 1 << " and chromosome " << index2 + 1 << '\n';
                    g << "  " << intermediatePopulation[index1].representation << ' ' << intermediatePopulation[index2].representation << " point " << breakPoint << '\n';
                }
                if (breakPoint) {    /// if break point is 0, the chromosomes will not be modified

                    string auxRepresentation = intermediatePopulation[index1].representation;

                    /// chromosome 1
                    intermediatePopulation[index1].representation = intermediatePopulation[index1].representation.substr(0, breakPoint) + intermediatePopulation[index2].representation.substr(breakPoint);
                    intermediatePopulation[index1].value = 0.0;     /// changing the representation and computing the new value and fitness
                    for (int j = chromosomeLength - 1; j >= 0; j--) {
                        if (intermediatePopulation[index1].representation[j] == '1') {    /// computing the value represented by the chromosome
                            intermediatePopulation[index1].value += pow(2, chromosomeLength - 1 - j);
                        }
                    }
                    intermediatePopulation[index1].value = computeValue(intermediatePopulation[index1].value);
                    intermediatePopulation[index1].fitness = f(intermediatePopulation[index1].value);

                    /// chromosome 2
                    intermediatePopulation[index2].representation = intermediatePopulation[index2].representation.substr(0, breakPoint) + auxRepresentation.substr(breakPoint);
                    intermediatePopulation[index2].value = 0.0;     /// changing the representation and computing the new value and fitness
                    for (int j = chromosomeLength - 1; j >= 0; j--) {
                        if (intermediatePopulation[index2].representation[j] == '1') {    /// computing the value represented by the chromosome
                            intermediatePopulation[index2].value += pow(2, chromosomeLength - 1 - j);
                        }
                    }
                    intermediatePopulation[index2].value = computeValue(intermediatePopulation[index2].value);
                    intermediatePopulation[index2].fitness = f(intermediatePopulation[index2].value);
                }
                if (step == 1) {
                    g << "  Result  " << intermediatePopulation[index1].representation << ' ' << intermediatePopulation[index2].representation << '\n';
                }
            }
        }
    }

    if(step == 1) {
        g << "\nAfter recombination\n";
        for (int i = 0; i < intermediatePopulation.size(); i++) {
            g << "  " << putSpaces(i + 1, intermediatePopulation.size()) << i + 1 << ". " << intermediatePopulation[i].representation << ",  x = ";
            if (intermediatePopulation[i].value > 0) {
                g << ' ';
            }
            g << intermediatePopulation[i].value << ",  f = " << intermediatePopulation[i].fitness << '\n';
        }
        g << '\n';
    }
}

void applyMutation() {

    if(step == 1) {
        g << "Mutation probability for each gene " << mutationProbability << "\nThe modified chromosomes are\n";
    }
    for(int i = 0; i < intermediatePopulation.size(); i++) {

        double u = (double)random() / RAND_MAX;     /// generate a random probability for the chromosome to suffer a mutation
        /// rare mutation

        if(u < mutationProbability) {

            if(step == 1) {
                g << "  " << putSpaces(i + 1, intermediatePopulation.size()) << i + 1 << '\n';
            }
            int position = (int)(random() % chromosomeLength);  /// generate a random position from the representation to apply a mutation

            intermediatePopulation[i].representation = intermediatePopulation[i].representation.substr(0, position) + to_string(1 - (intermediatePopulation[i].representation[position] - '0')) + intermediatePopulation[i].representation.substr(position + 1);
            intermediatePopulation[i].value = 0.0;
            for(int j = chromosomeLength - 1; j >= 0; j--) {
                if(intermediatePopulation[i].representation[j] == '1') {    /// computing the value represented by the chromosome
                    intermediatePopulation[i].value += pow(2, chromosomeLength - 1 - j);
                }
            }
            intermediatePopulation[i].value = computeValue(intermediatePopulation[i].value);
            intermediatePopulation[i].fitness = f(intermediatePopulation[i].value);
        }
    }

    intermediatePopulation.push_back(elitist);  /// adding the elitist chromosome in the list
    if(step == 1) {
        g << "\nAfter mutation\n";
        for (int i = 0; i < intermediatePopulation.size(); i++) {

            g << "  " << putSpaces(i + 1, intermediatePopulation.size()) << i + 1 << ". " << intermediatePopulation[i].representation << ",  x = ";
            if (intermediatePopulation[i].value > 0) {
                g << ' ';
            }
            g << intermediatePopulation[i].value << ",  f = " << intermediatePopulation[i].fitness << '\n';
        }
        g << "\nMaximum & average performance evolution\n";
    }
}

int main() {

    read();
    srandom(time(nullptr));

    generateChromosomes();
    vector<Chromosome> elitists;
    for(step = 1; step <= stages; step++) {

        elitist = chromosomes[0];   /// finding the elitist chromosome, which will automatically pass in the next generation
        for(int i = 0; i < populationDimension; i++) {
            if(chromosomes[i].fitness > elitist.fitness) {
                elitist = chromosomes[i];
            }
        }
        selectionProbabilities();
        selectChromosomes();
        selectForCrossover();
        applyMutation();

        /// compute maximum fitness
        double maxValue = 0.0;
        double x;
        totalFitness = 0.0;
        for(auto & i : intermediatePopulation) {

            if(i.fitness > maxValue) {
                maxValue = i.fitness;
                x = i.value;    /// the value where maximum fitness is obtained
            }
            totalFitness += i.fitness;
        }

        elitists.push_back(elitist);
        int size = elitists.size();
        if(size >= stages / 10) {
            if(elitists[size - 1].fitness == elitists[size - 1 - stages / 10].fitness) {    /// stop the algorithm if the elitist is the same for 10% of stages
                g << "\nSolution " << maxValue << '\n';
                break;
            }
        }

        g << "  " << putSpaces(step, stages) << step << ". maximum = " << maxValue << ", x = " << x << ",  average performance = " << totalFitness / (1.0 * populationDimension) << '\n';
        chromosomes = intermediatePopulation;   /// the current intermediatePopulation will be the initial set of chromosomes for the next generation

        if(step == stages) {
            g << "\nSolution " << maxValue << '\n';
        }
    }
    g.close();

    return 0;
}
