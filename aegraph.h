// Copyright 2019 Roinita Lucian, Calin Dragos
#ifndef AEGRAPH_H_
#define AEGRAPH_H_

#include <vector>
#include <string>

class AEGraph {
 public:
    explicit AEGraph(std::string representation);

    std::string repr() const;

    void sort();

    bool operator<(const AEGraph& other) const;
    bool operator==(const AEGraph& other) const;
    bool operator!=(const AEGraph& other) const;
    AEGraph operator[](const int index) const;

    bool contains(const AEGraph& other) const;
    bool contains(const std::string other) const;

    int num_subgraphs() const;
    int num_atoms() const;
    int size() const;

    std::vector<std::vector<int>> possible_double_cuts() const;
    AEGraph double_cut(std::vector<int> where) const;

    std::vector<std::vector<int>> possible_erasures(int level = -1);
    AEGraph erase(std::vector<int>);

    std::vector<std::vector<int>> possible_deiterations() const;
    AEGraph deiterate(std::vector<int> where);
    std::vector<std::vector<int>> get_paths_to(const std::string other) const;
    std::vector<std::vector<int>> get_paths_to(const AEGraph& other) const;
    void delete_down(AEGraph* g, int pos, bool first);
    std::vector<std::vector<int>>
    possible_erasures_recursive(AEGraph root, int level, int first);

    std::vector<std::string> atoms;
    std::vector<AEGraph> subgraphs;

    friend std::ostream& operator<<(std::ostream &out, const AEGraph &g);

    bool is_SA;
};

#endif  // AEGRAPH_H_
