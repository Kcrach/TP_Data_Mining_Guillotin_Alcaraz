#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

bool contain_pattern(string, string);
void init_patterns(vector<string>, vector<pair<string, int>> &);
bool pattern_exist(string , vector<pair<string, int>>);
void display_freq_pattern(vector<pair<string, int>>);
int count_pattern_in_transitions(vector<string>, string);
void find_frequent_patterns(vector<string>, vector<pair<string, int>> &, int);
void merge_freq_pattern(vector<pair<string, int>> &);

int main()
{
    int minsup = 2;
    vector<string> transitions;
    vector<pair<string, int>> freq_pattern;

    string t1 = "DEF";
    string t2 = "BCE";
    string t3 = "ABCE";
    string t4 = "BCE";
    string t5 = "ABCE";
    string t6 = "BCE";

    transitions.push_back(t1);
    transitions.push_back(t2);
    transitions.push_back(t3);
    transitions.push_back(t4);
    transitions.push_back(t5);
    transitions.push_back(t6);

    find_frequent_patterns(transitions, freq_pattern, minsup);
    display_freq_pattern(freq_pattern);
}

//Function to know if a pattern is into a transition
bool contain_pattern(string transition, string pattern) {
    if (transition.size() < pattern.size()) {
        return false;
    }

    for (int i = 0; i < pattern.size(); i++) {
        if (transition.find(pattern.at(i)) == string::npos) {
            return false;
        }
    }

    return true;
}

//Function to init patterns with size 1 (A,B,C,D....)
void init_patterns(vector<string> transitions, vector<pair<string, int>> &freq_pattern) {
    for (int i = 0; i < transitions.size(); i++) {
        for (int j = 0; j < transitions.at(i).size(); j++) {
            string p;
            stringstream ss;
            ss << transitions.at(i).at(j);
            ss >> p;
            bool exist = pattern_exist(p, freq_pattern);
            if (exist == false) {
                pair<string, int> tmp (p, 0);
                freq_pattern.push_back(tmp);
            }
        }
    }

    sort(freq_pattern.begin(), freq_pattern.end());
}

//Check if a pattern is already in the vector freq_pattern
bool pattern_exist(string pattern, vector<pair<string, int>> freq_pattern) {
    for (int i = 0; i < freq_pattern.size(); i++) {
        if(freq_pattern.at(i).first == pattern) {
            return true;
        }
    }

    return false;
}

//Display the vector freq_pattern
void display_freq_pattern(vector<pair<string, int>> freq_pattern){
    for (int i = 0; i < freq_pattern.size(); i++) {
        cout << freq_pattern.at(i).first << " : " << freq_pattern.at(i).second << endl;
    }
}

//Count frequent of a pattern
int count_pattern_in_transitions(vector<string> transitions, string pattern) {
    int cpt = 0;

    for (int i = 0; i < transitions.size(); i++) {
        if (contain_pattern(transitions.at(i),pattern)) {
            cpt++;
        }
    }

    return cpt;
}

//Algorithm to find frequent patterns
void find_frequent_patterns(vector<string> transitions, vector<pair<string, int>> &freq_pattern, int minsup) {
    vector<pair<string, int>> tmp = freq_pattern;

    //First init the vector
    init_patterns(transitions, tmp);
    while (tmp.size() >= 1) {
        for (int i = 0; i < tmp.size(); i++) {
            int cpt = count_pattern_in_transitions(transitions, tmp.at(i).first);
            if (cpt >= minsup) {
                tmp.at(i).second = cpt;
            } else {
                tmp.erase(tmp.begin()+i);
                i--;
            }
        }
        //At the end of the first round we add frequent pattern in the final list
        freq_pattern.insert(freq_pattern.begin() + freq_pattern.size(), tmp.begin(), tmp.end());
        merge_freq_pattern(tmp);
    }
}

//Merge pattern
void merge_freq_pattern(vector<pair<string, int>> &freq_pattern_tmp) {
    int size_v = freq_pattern_tmp.size();

    if (freq_pattern_tmp.at(0).first.size() == 1) { //Pattern with one char (first round)
        for (int i = 0; i < size_v; i++) {
            for (int j = i+1; j < size_v; j++) {
                string new_pattern = freq_pattern_tmp.at(i).first + freq_pattern_tmp.at(j).first;
                pair<string, int> tmp (new_pattern, 0);
                freq_pattern_tmp.push_back(tmp);
            }
        }

        for (int i = 0; i < size_v; i++) {
            freq_pattern_tmp.erase(freq_pattern_tmp.begin());
        }
        return;
    } else {
        for (int i = 0; i < size_v; i++) {
            //Prefix of first pattern
            string tmp1 = freq_pattern_tmp.at(i).first;
            tmp1.pop_back();
            for (int j = i+1; j < size_v; j++) {
                //Prefix of second pattern
                string tmp2 = freq_pattern_tmp.at(j).first;
                tmp2.pop_back();
                if (tmp1 == tmp2) {
                    string p1 = freq_pattern_tmp.at(i).first;
                    string p2 = freq_pattern_tmp.at(j).first;
                    string new_pattern = p1 + p2.at(p2.size()-1);
                    pair<string, int> tmp (new_pattern, 0);
                    freq_pattern_tmp.push_back(tmp);
                }
            }
        }

        for (int i = 0; i < size_v; i++) {
            freq_pattern_tmp.erase(freq_pattern_tmp.begin());
        }
    }
}
