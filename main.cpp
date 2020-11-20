#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <fstream>

typedef std::vector<int> Transition;
typedef std::pair<Transition, double> Pattern;
typedef std::vector<std::pair<Transition, double>> PatternVector;

bool contain_pattern(const Transition &, const Transition &);
void init_patterns(const std::vector<std::string> &, PatternVector &);
bool pattern_exist(const std::string &, const PatternVector &);
void display_freq_pattern(const PatternVector &);
PatternVector merge_freq_pattern(const PatternVector &);
std::vector<std::string> generate_transitions(const std::string &, int);
PatternVector calc_weight_transitions(const std::vector<Transition> &);
Transition choose_transition(PatternVector);
Transition choose_subset_pattern(Transition);
PatternVector sampling_frequency(const std::vector<Transition> &);
void count_frequency_sample(PatternVector &, std::vector<Transition>);
void count_area_sample(PatternVector &, std::vector<Transition>);
PatternVector calc_weight_transitions_2(const std::vector<Transition> &);
int draw_k (std::string);
PatternVector sampling_area(const std::vector<Transition> &);
std::vector<Transition> convert_file_in_dataset(std::string);
int rng(int);

int main() {
    int minsup = 2;
    std::vector<Transition> transitions;
    PatternVector freq_pattern;
    PatternVector sample_frequency;
    PatternVector sample_area;

    /*Transition t1 = {1, 2, 45};
    Transition t2 = {1, 2, 3};
    Transition t3 = {1, 2, 5};
    Transition t4 = {1, 2, 7};
    Transition t5 = {1, 3, 6, 8};

    transitions.push_back(t1);
    transitions.push_back(t2);
    transitions.push_back(t3);
    transitions.push_back(t4);
    transitions.push_back(t5);*/

    bool VERBOSE = true;

    //transitions = generate_transitions("ABCDEFGHIJ", 5000);

    std::chrono::time_point<std::chrono::system_clock> start_t, end_t; // chrono
    start_t = std::chrono::system_clock::now();

    end_t = std::chrono::system_clock::now();
    int elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>
                             (end_t - start_t).count();

    if (VERBOSE) {
        std::sort(freq_pattern.begin(),
                  freq_pattern.end(),
                  [](const Pattern & a, const Pattern & b) {
                      return a.first.size() == b.first.size()?
                      a.first < b.first: a.first.size() < b.first.size();
                  });

            //Test algo 1
            /*std::cout << "Frequence :\n";
            sample_frequency = sampling_frequency(transitions);
            count_frequency_sample(sample_frequency, transitions);
            display_freq_pattern(sample_frequency);
            std::cout << std::endl;*/

            //Test algo 2
            /*std::cout << "Aire :\n";
            sample_area = sampling_area(transitions);
            count_area_sample(sample_area, transitions);
            display_freq_pattern(sample_area);*/

            //Q4
            transitions = convert_file_in_dataset("mushrooms.txt");

            /*sample_frequency = sampling_frequency(transitions);
            count_frequency_sample(sample_frequency, transitions);
            display_freq_pattern(sample_frequency);*/

            sample_area = sampling_area(transitions);
            count_area_sample(sample_area, transitions);
            display_freq_pattern(sample_area);
    }
    std::cout << "Execution time: " << elapsed_ms << "ms" << std::endl;


}

//Algorithm 1 : Sampling by frequency, step 1 : calc weight
PatternVector calc_weight_transitions(const std::vector<Transition> & transitions) {
    PatternVector weights;

    for (int i = 0; i < transitions.size(); i++) {
        int tmp = pow(2, transitions[i].size());
        std::pair<Transition, int> p (transitions[i], tmp);
        weights.push_back(p);
    }

    return weights;
}

//Algorithm 1 & 2: Sampling by frequency, step 2 : Choose a transition
Transition choose_transition(PatternVector w) {
    int total_weight = 0;

    for (int i = 0; i < w.size(); i++) {
        total_weight += w[i].second;
    }

    int cumul_proba = 0;
    int p = rng(total_weight);

    for (int i = 0; i < w.size(); i++) {
        cumul_proba += w[i].second;
        if (p <= cumul_proba) {
            return w[i].first;
        }
    }
}

//Algorithm 1 : Sampling by frequency, step 3 : make the final set
PatternVector sampling_frequency(const std::vector<Transition> & transitions) {
    PatternVector weights;
    weights = calc_weight_transitions(transitions);

    PatternVector sample;
    bool exist = false;

    for (int i = 0; i < transitions.size(); i++) {
        Transition choosen_t = choose_transition(weights);
        Transition subset_choosen_t = choose_subset_pattern(choosen_t);
        for (int j = 0; j < sample.size(); j++) {
            if (sample[j].first == subset_choosen_t) {
                exist = true;
                break;
            }
        }
        if (exist == false) {
            Pattern p (subset_choosen_t, 0);
            sample.push_back(p);
        } else {
            exist = false;
        }
    }

    return sample;
}

//Algorithm 1 : Choose uniformaly a subset pattern of a transition
Transition choose_subset_pattern(Transition pattern) {
    Transition choosen_subset = {};

    while (choosen_subset.size() == 0) {
        for (int i = 0; i < pattern.size(); i++) {
            int r = rng(2);
            Transition tmp = {pattern[i]};
            if (r == 0 && !contain_pattern(choosen_subset, tmp)) {
                choosen_subset.push_back(pattern[i]);
            }
        }
    }

    std::sort(choosen_subset.begin(), choosen_subset.end());

    return choosen_subset;
}

//Algorithm 2 : Sampling by area, step 1 : calc weight
PatternVector calc_weight_transitions_2(const std::vector<Transition> & transitions) {
    PatternVector weights;

    for (int i = 0; i < transitions.size(); i++) {
        int tmp = transitions[i].size() * pow(2, transitions[i].size() - 1);
        std::pair<Transition, int> p (transitions[i], tmp);
        weights.push_back(p);
    }

    return weights;
}

//Algorithm 2 : Sampling by area, step 3 : draw k
int draw_k (Transition choosen_t) {
    std::vector<int> id;
    int total = 0;

    for (int i = 1; i < choosen_t.size()+1; i++) {
        id.push_back(i);
        total += i;
    }

    int r = rng(total) + 1;
    int p = 0;

    for (int i =0 ; i < id.size(); i++) {
        p += id[i];
        if (p >= r) {
            return id[i];
        }
    }
}

//Algorithm 2 : Choose uniformaly a subset pattern of size k into the transition
Transition choose_subset_pattern_2(Transition pattern, int k) {
    Transition choosen_subset = {};

    while (choosen_subset.size() < k) {
        for (int i = 0; i < pattern.size(); i++) {
            int r = rng(2);
            Transition tmp = {pattern[i]};
            if (r == 0 && !contain_pattern(choosen_subset, tmp)) {
                choosen_subset.push_back(pattern[i]);
                break;
            }
        }
    }

    std::sort(choosen_subset.begin(), choosen_subset.end());

    return choosen_subset;
}

//Algorithm 2 : Sampling by area, step 3 : make the final set
PatternVector sampling_area(const std::vector<Transition> & transitions) {
    PatternVector weights;
    weights = calc_weight_transitions_2(transitions);

    PatternVector sample;
    bool exist = false;

    for (int i = 0; i < transitions.size(); i++) {
        Transition choosen_t = choose_transition(weights);
        int k = draw_k(choosen_t);
        Transition subset_choosen_t = choose_subset_pattern_2(choosen_t, k);
        for (int j = 0; j < sample.size(); j++) {
            if (sample[j].first == subset_choosen_t) {
                exist = true;
                break;
            }
        }
        if (exist == false) {
            Pattern p (subset_choosen_t, 0);
            sample.push_back(p);
        } else {
            exist = false;
        }
    }

    return sample;
}

//Q3 : Count frequency of a sample
void count_frequency_sample(PatternVector & sample, std::vector<Transition> transitions) {
    for (int i = 0; i < transitions.size(); i++) {
        for (int j = 0; j < sample.size(); j++) {
            if (contain_pattern(transitions[i], sample[j].first)) {
                sample[j].second++;
            }
        }
    }
}

//Q3 : Count area of a sample
void count_area_sample(PatternVector & sample, std::vector<Transition> transitions) {
    for (int i = 0; i < transitions.size(); i++) {
        for (int j = 0; j < sample.size(); j++) {
            if (contain_pattern(transitions[i], sample[j].first)) {
                sample[j].second += sample[j].first.size();
            }
        }
    }
}

//Q4 : Read a txt file and create Dataset
std::vector<Transition> convert_file_in_dataset(std::string nameFile) {
    std::ifstream is(nameFile);  //Ouverture d'un fichier en lecture
    std::vector<Transition> transitions;

    if(is)
    {
        std::string data_record;
        std::string delimiter = " ";
        std::string e;
        size_t pos = 0;
        while (getline(is, data_record)) {
            pos = 0;
            Transition transition;
            while ((pos = data_record.find(delimiter)) != std::string::npos) {
                e = data_record.substr(0, pos);
                transition.push_back(std::stoi(e));
                data_record.erase(0, pos + delimiter.length());
            }
            transitions.push_back(transition);
        }
        return transitions;
    }
    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << std::endl;
    }
}

std::vector<std::string> generate_transitions(const std::string & symbols,
                                              int nb_transitions) {
    std::vector<std::string> transitions;
    std::string transition;
    std::string tmp_symbols(symbols);
    bool exist;
    for (int i(0); i < nb_transitions; ++i) {
        int nb_item = 1 + rand() % symbols.size();
        std::random_shuffle(tmp_symbols.begin(), tmp_symbols.end());
        transition = tmp_symbols.substr(0, nb_item);
        std::sort(transition.begin(), transition.end());
        exist = false;
        for (const std::string & t: transitions) {
            if (t == transition) {
                exist = true;
                break;
            }
        }
        if (!exist) {
            transitions.push_back(transition);
        }
    }
    return transitions;
}

//Function to know if a pattern is into a transition
bool contain_pattern(const Transition & transition,
                     const Transition & pattern) {

    if (transition.size() == 0) {
        return false;
    }

    for (unsigned i(0); i < pattern.size(); i++) {
        bool is_on = false;
        for (int j = 0; j < transition.size(); j++) {
            if (transition[j] == pattern[i]) {
               is_on = true;
            }
        }

        if (!is_on){
            return false;
        }
    }

    return true;
}

//Display the std::vector freq_pattern
void display_freq_pattern(const PatternVector & freq_pattern){
    for (unsigned i(0); i < freq_pattern.size(); ++i) {
        for (int j = 0; j < freq_pattern[i].first.size(); j++) {
            std::cout << freq_pattern[i].first[j] << " ";
        }

        std::cout << ": " << freq_pattern.at(i).second << std::endl;
    }
}

int rng(int mod) {
    unsigned int seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
    std::minstd_rand0 generator(seed + rand() % 1000);
    return (int) generator() % mod;
}
