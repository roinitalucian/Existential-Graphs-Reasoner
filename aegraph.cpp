// Copyright 2019 Roinita Lucian, Calin Dragos
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    // returns all paths in the three that lead to a node followed
    // by two subgraphs in a row (and no atoms)
    std::vector<std::vector<int>> paths;

    for (int i = 0; i < num_subgraphs(); i++) {
        if (subgraphs[i].num_subgraphs() == 1 &&
            subgraphs[i].num_atoms() == 0) {
            paths.push_back({i});
        }
    }

    for (int i = 0; i < num_subgraphs(); i++){
        auto r = subgraphs[i].possible_double_cuts();
        for (auto& v : r)
            v.insert(v.begin(), i);
        copy(r.begin(), r.end(), back_inserter(paths));
    }

    return paths;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    // removes the double negation and moves the atoms below up
    AEGraph result(repr());
    AEGraph *curr = &result;
    for (unsigned int i = 0; i < where.size() - 1; i++) {
        curr = &(curr->subgraphs[where[i]]);
    }
    int pos = where[where.size() - 1];
    curr->subgraphs[pos] = curr->subgraphs[pos].subgraphs[0];
    for (auto a : curr->subgraphs[pos].atoms) {
        curr->atoms.push_back(a);
    }
    if (curr->subgraphs[pos].num_subgraphs() > 0) {
        curr = &(curr->subgraphs[pos]);
    }
    curr->subgraphs[pos].atoms.clear();
    curr->subgraphs.erase(curr->subgraphs.begin() + pos);
    return result;
}

std::vector<std::vector<int>> AEGraph::possible_erasures_recursive
(AEGraph root, int level, int first){
    // iterates through all the subraphs and atoms and pushes to
    // the result vector the paths that have odd lengths
    std::vector<std::vector<int>> result;
    for (int i = 0; i < num_subgraphs(); i++) {
        std::vector<std::vector<int>> paths = root.get_paths_to(subgraphs[i]);
        if (first) {
            result.push_back({i});
        }
        for (auto path : paths) {
            if (path.size() % 2 == 1 && path.size() != 1) {
                result.push_back(path);
            }
        }

        paths = subgraphs[i].possible_erasures_recursive(root, level, 0);
        copy(paths.begin(), paths.end(), back_inserter(result));
    }

    for (int i = 0; i < num_atoms(); i++) {
        std::vector<std::vector<int>> paths = root.get_paths_to(atoms[i]);
        if (first) {
            result.push_back({i + num_subgraphs()});
        }
        for (auto path : paths) {
            if (path.size() % 2 == 1 && path.size() != 1) {
                result.push_back(path);
            }
        }
    }

    return result;
}

std::vector<std::vector<int>> AEGraph::possible_erasures(int level) {
    // initialize the root
    AEGraph root(repr());
    return possible_erasures_recursive(root, level, 1);
}


AEGraph AEGraph::erase(std::vector<int> where) {
    // identical to deiterate
    return deiterate(where);
}


std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    // calls get_path_to for each subgraph and atom and saves
    // the paths that can be deiterated
    std::vector<std::vector<int>> result;
    AEGraph current(repr());
    for (int i = 0; i < current.num_subgraphs(); i++) {
        std::vector<std::vector<int>> paths = current.get_paths_to(current[i]);
        for (auto path : paths) {
            if (path[0] != i) {
                result.push_back(path);
            }
        }
    }
    for (int i = num_subgraphs(); i < current.size(); i++) {
        std::string repr = current[i].repr();
        std::vector<std::vector<int>> paths =
        current.get_paths_to(current[i].repr().substr(1, 1));
        for (auto path : paths) {
            if (path[0] != i) {
                result.push_back(path);
            }
        }
    }
    return result;
}

void AEGraph::delete_down(AEGraph* g, int pos, bool first) {
    // deletes everything starting at graph g
    if (first) {
        if (g->num_subgraphs() > pos) {
            if (g->subgraphs[pos].size() > 0) {
                delete_down(&(g->subgraphs[pos]), pos, 0);
                g->subgraphs.erase(g->subgraphs.begin()+pos);
            }
        } else {
            g->atoms.erase(g->atoms.begin() + pos - g->num_subgraphs());
        }
    } else {
        for (int i = 0; i < g->num_subgraphs(); i++) {
            delete_down(&(g->subgraphs[i]), i, 0);
            g->subgraphs.erase(g->subgraphs.begin()+i);
        }
        for (int i = num_subgraphs(); i < g->size(); i++) {
            g->atoms.clear();
        }
    }
}

AEGraph AEGraph::deiterate(std::vector<int> where)  {
    // goes to the graph that where leads to and apply delete_down to it
    AEGraph result(repr());
    AEGraph *curr = &result;
    for (unsigned int i = 0; i < where.size() - 1; i++) {
        curr = &(curr->subgraphs[where[i]]);
    }
    delete_down(curr, where[where.size() - 1], 1);
    return result;
}

