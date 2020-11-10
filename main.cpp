#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>

typedef std::pair<std::string, int> Pattern;
typedef std::vector<std::pair<std::string, int>> PatternVector;

bool contain_pattern(const std::string &, const std::string &);
void init_patterns(const std::vector<std::string> &, PatternVector &);
bool pattern_exist(const std::string &, const PatternVector &);
void display_freq_pattern(const std::vector<std::pair<std::string, int>> &);
int count_pattern_in_transitions(const std::vector<std::string> &,
                                 const std::string &);
PatternVector find_frequent_patterns(const std::vector<std::string> &,
                                     const PatternVector &, int);
PatternVector merge_freq_pattern(const PatternVector &);
std::vector<std::string> generate_transitions(const std::string &, int);
std::vector<std::pair<std::string, int>> calc_weight_transitions(const std::vector<std::string> &);
std::string choose_transition(std::vector<std::pair<std::string, int>>);
std::string choose_subset_pattern(std::string);
PatternVector sampling_frequency(const std::vector<std::string> &);
int rng(int);

int main() {
    int minsup = 2;
    std::vector<std::string> transitions;
    PatternVector freq_pattern;
    PatternVector sample_frequency;

    std::string t1 = "DEF";
    std::string t2 = "BCE";
    std::string t3 = "ABCE";
    std::string t4 = "ABC";
    std::string t5 = "ABCE";
    std::string t6 = "BCE";

    transitions.push_back(t1);
    transitions.push_back(t2);
    transitions.push_back(t3);
    transitions.push_back(t4);
    transitions.push_back(t5);
    transitions.push_back(t6);

    bool VERBOSE = false;

    //transitions = generate_transitions("ABCDEFGHIJ", 5000);

    std::chrono::time_point<std::chrono::system_clock> start_t, end_t; // chrono
    start_t = std::chrono::system_clock::now();

    freq_pattern = find_frequent_patterns(transitions, freq_pattern, minsup);
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
        display_freq_pattern(freq_pattern);
    }
    std::cout << "Execution time: " << elapsed_ms << "ms" << std::endl;

    //Test algo 1
    sample_frequency = sampling_frequency(transitions);
    display_freq_pattern(sample_frequency);
}


//Algorithm 1 : Sampling by frequency, step 1 : calc weight
std::vector<std::pair<std::string, int>> calc_weight_transitions(const std::vector<std::string> & transitions) {
    std::vector<std::pair<std::string, int>> weights;

    for (int i = 0; i < transitions.size(); i++) {
        int tmp = pow(2, transitions[i].size());
        std::pair<std::string, int> p (transitions[i], tmp);
        weights.push_back(p);
    }

    return weights;
}

//Algorithm 1 : Sampling by frequency, step 2 : Choose a transition
std::string choose_transition(std::vector<std::pair<std::string, int>> w) {
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
PatternVector sampling_frequency(const std::vector<std::string> & transitions) {
    std::vector<std::pair<std::string, int>> weights;
    weights = calc_weight_transitions(transitions);

    PatternVector sample;
    bool exist = false;

    for (int i = 0; i < transitions.size(); i++) {
        std::string choosen_t = choose_transition(weights);
        std::string subset_choosen_t = choose_subset_pattern(choosen_t);
        for (int j = 0; j < sample.size(); j++) {
            if (sample[j].first == subset_choosen_t) {
                sample[j].second++;
                exist = true;
                break;
            }
        }
        if (exist == false) {
            std::pair<std::string, int> p (subset_choosen_t, 1);
            sample.push_back(p);
        } else {
            exist = false;
        }
    }

    return sample;
}

//Algorithm 1 : Choose uniformaly a subset pattern of a transition
std::string choose_subset_pattern(std::string pattern) {
    std::vector<std::string> possible_patterns;

    for (int i = 0; i < pattern.size(); i++) {
        std::string s(1,pattern[i]);
        possible_patterns.push_back(s);
    }

    std::sort(possible_patterns.begin(),
                possible_patterns.end());

    for (int i = 0; i < possible_patterns.size(); i++) {
        for (int j = i+1; j < possible_patterns.size(); j++) {
            std::string new_tmp = "";
            if (possible_patterns[i].size() == possible_patterns[j].size()
                    && possible_patterns[i].size() == 1) {
                new_tmp = possible_patterns[i] + possible_patterns[j];
            } else if (possible_patterns[i].size() == possible_patterns[j].size()
                        && possible_patterns[i].size() > 1) {
                std::string tmp1 = possible_patterns[i];
                tmp1.pop_back();

                std::string tmp2 = possible_patterns[j];
                tmp2.pop_back();

                if (tmp1 == tmp2) {
                    new_tmp = possible_patterns[i] + possible_patterns[j][possible_patterns[j].size() - 1];
                }
            }

            if (new_tmp != "") {
                possible_patterns.push_back(new_tmp);
            }

            if (new_tmp == pattern) {
                break;
            }
        }
    }

    int r = rng(possible_patterns.size() - 1);
    return possible_patterns[r];
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
bool contain_pattern(const std::string & transition,
                     const std::string & pattern) {
    if (transition.size() < pattern.size()) {
        return false;
    }

    for (unsigned i(0); i < pattern.size(); ++i) {
        if (transition.find(pattern[i]) == std::string::npos) {
            return false;
        }
    }

    return true;
}

//Function to init patterns with size 1 (A,B,C,D....)
void init_patterns(const std::vector<std::string> & transitions,
                   PatternVector & freq_pattern) {
    for (const std::string & transition: transitions) {
        for (char c: transition) {
            std::string p(1, c);
            if (!pattern_exist(p, freq_pattern)) {
                freq_pattern.push_back(Pattern(p, 0));
            }
        }
    }

    std::sort(freq_pattern.begin(),
                  freq_pattern.end(),
                  [](const Pattern & a, const Pattern & b) {
                      return a.first.size() == b.first.size()?
                      a.first < b.first: a.first.size() < b.first.size();
                  });
}

//Check if a pattern is already in the std::vector freq_pattern
bool pattern_exist(const std::string & pattern,
                   const PatternVector & freq_pattern) {
    for (Pattern p: freq_pattern) {
        if(p.first == pattern) {
            return true;
        }
    }

    return false;
}

//Display the std::vector freq_pattern
void display_freq_pattern(const PatternVector & freq_pattern){
    for (unsigned i(0); i < freq_pattern.size(); ++i) {
        std::cout << freq_pattern.at(i).first
        << " : " << freq_pattern.at(i).second << std::endl;
    }
}

//Count frequent of a pattern
int count_pattern_in_transitions(const std::vector<std::string> & transitions,
                                 const std::string & pattern) {
    int cpt = 0;
    for (const std::string & transition: transitions) {
        cpt += contain_pattern(transition, pattern);
    }
    return cpt;
}

//Algorithm to find frequent patterns
PatternVector find_frequent_patterns(const std::vector<std::string> & transitions,
                                     const PatternVector & freq_pattern,
                                     int minsup) {
    PatternVector tmp(freq_pattern);
    PatternVector result(freq_pattern);
    //First init the std::vector
    init_patterns(transitions, tmp);
    int cpt;
    while (!tmp.empty()) {
        std::vector<bool> ignore(tmp.size(), false);
        for (unsigned i(0); i < tmp.size(); ++i) {
            cpt = count_pattern_in_transitions(transitions, tmp[i].first);
            if (cpt >= minsup) {
                tmp[i].second = cpt;
            } else {
                ignore[i] = true;
            }
        }
        PatternVector new_tmp;
        for (unsigned i(0); i < tmp.size(); ++i) {
            if (!ignore[i]) {
                new_tmp.push_back(tmp[i]);
            }
        }
        tmp = new_tmp;
        //At the end of the first round we add frequent pattern in the result
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(result));
        tmp = merge_freq_pattern(tmp);
    }
    return result;
}

//Merge pattern
PatternVector merge_freq_pattern(const PatternVector & freq_pattern_tmp) {
    int size_v(freq_pattern_tmp.size());
    PatternVector result(freq_pattern_tmp);
    if (result.empty()) {
        return result;
    }
    if (result[0].first.size() == 1) { //Pattern with one char (first round)
        for (int i(0); i < size_v; ++i) {
            for (int j(i + 1); j < size_v; ++j) {
                result.push_back(Pattern(std::string(result[i].first
                                                     + result[j].first), 0));
            }
        }
    } else {
        std::string tmp1, tmp2;
        for (int i(0); i < size_v; ++i) {
            //Prefix of first pattern
            tmp1 = result[i].first;
            tmp1.pop_back();
            for (int j(i + 1); j < size_v; ++j) {
                //Prefix of second pattern
                tmp2 = result[j].first;
                tmp2.pop_back();
                if (tmp1 == tmp2) {
                    result.push_back(Pattern(std::string(result[i].first
                                                         + result[j].first.back()), 0));
                }
            }
        }
    }

    return PatternVector(result.begin() + size_v, result.end());
}

int rng(int mod) {
    unsigned int seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
    std::minstd_rand0 generator(seed + rand() % 1000);
    return (int) generator() % mod;
}
