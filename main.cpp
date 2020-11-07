#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>

typedef std::pair<std::string, int> Pattern;
typedef std::vector<std::pair<std::string, int>> PatternVector;

bool contain_pattern(const std::string &, const std::string &);
void init_patterns(std::vector<std::string>, PatternVector &);
bool pattern_exist(const std::string &, const PatternVector &);
void display_freq_pattern(const std::vector<std::pair<std::string, int>> &);
int count_pattern_in_transitions(const std::vector<std::string> &,
                                 const std::string &);
PatternVector find_frequent_patterns(const std::vector<std::string> &,
                                     const PatternVector &, int);
PatternVector merge_freq_pattern(const PatternVector &);
std::vector<std::string> generate_transitions(const std::string &, int);

int main() {
    int minsup = 2;
    std::vector<std::string> transitions;
    PatternVector freq_pattern;

    /*std::string t1 = "DEF";
    std::string t2 = "BCE";
    std::string t3 = "ABCE";
    std::string t4 = "BCE";
    std::string t5 = "ABCE";
    std::string t6 = "BCE";

    transitions.push_back(t1);
    transitions.push_back(t2);
    transitions.push_back(t3);
    transitions.push_back(t4);
    transitions.push_back(t5);
    transitions.push_back(t6);*/

    bool VERBOSE = true;

    transitions = generate_transitions("ABCDEFGHIJ", 5000);

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

    for (int i = 0; i < pattern.size(); i++) {
        if (transition.find(pattern.at(i)) == std::string::npos) {
            return false;
        }
    }

    return true;
}

//Function to init patterns with size 1 (A,B,C,D....)
void init_patterns(std::vector<std::string> transitions,
                   PatternVector & freq_pattern) {
    for (int i(0); i < transitions.size(); ++i) {
        for (int j(0); j < transitions[i].size(); ++j) {
            std::string p;
            std::stringstream ss;
            ss << transitions[i][j];
            ss >> p;
            if (!pattern_exist(p, freq_pattern)) {
                Pattern tmp (p, 0);
                freq_pattern.push_back(tmp);
            }
        }
    }

    std::sort(freq_pattern.begin(), freq_pattern.end());
}

//Check if a pattern is already in the std::vector freq_pattern
bool pattern_exist(const std::string & pattern,
                   const PatternVector & freq_pattern) {
    for (int i = 0; i < freq_pattern.size(); i++) {
        if(freq_pattern.at(i).first == pattern) {
            return true;
        }
    }

    return false;
}

//Display the std::vector freq_pattern
void display_freq_pattern(const PatternVector & freq_pattern){
    for (int i = 0; i < freq_pattern.size(); ++i) {
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
    PatternVector tmp = freq_pattern;
    PatternVector result = freq_pattern;
    //First init the std::vector
    init_patterns(transitions, tmp);
    int cpt;
    while (tmp.size() >= 1) {
        std::vector<bool> ignore(tmp.size(), false);
        for (int i(0); i < tmp.size(); ++i) {
            cpt = count_pattern_in_transitions(transitions, tmp[i].first);
            if (cpt >= minsup) {
                tmp[i].second = cpt;
            } else {
                ignore[i] = true;
            }
        }
        PatternVector new_tmp;
        for (int i(0); i < tmp.size(); ++i) {
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
    int size_v = freq_pattern_tmp.size();
    PatternVector result(freq_pattern_tmp);
    if (result.size() == 0) {
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
