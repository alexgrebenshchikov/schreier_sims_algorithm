//
// Created by alex on 04.05.2020.
//

#include <iostream>
#include <algorithm>
#include "schreier.h"


Perm::Perm(std::vector<int> p) {
    _data = std::move(p);
    _size = _data.size();
}

Perm::Perm(size_t n) {
    _data.resize(n);
    _size = n;
    for (int i = 0; i < _size; i++)
        _data[i] = i;
}


bool Perm::operator<(const Perm &p) const {
    for (int i = 0; i < _data.size(); i++) {
        if (_data[i] < p._data[i])
            return true;
        if (_data[i] > p._data[i])
            return false;
    }
    return false;
}

bool Perm::operator==(const Perm &p) const {
    return !(*this < p || p < *this);
}


bool Perm::operator!=(const Perm &p) const {
    return !(*this == p);
}

Perm &Perm::operator*=(const Perm &p) {
    for (int &i : _data) {
        i = p._data[i];
    }
    return *this;
}

void Perm::print() const {
    for (auto &i:_data)
        std::cout << i << ' ';
    std::cout << '\n';
}

Perm operator*(Perm p1, const Perm &p2) {
    p1 *= p2;
    return p1;
}

int &Perm::operator[](size_t pos) {
    return _data[pos];
}

const int &Perm::operator[](size_t pos) const {
    return _data[pos];
}

Perm inverse(const Perm &p) {
    Perm res = p;
    for (int i = 0; i < p._data.size(); i++)
        res[p[i]] = i;
    return res;
}

Perm::operator bool() const {
    return !(_data.empty());
}

size_t Perm::get_size() const {
    return _size;
}


void build_schreier_tree(int w, std::set<Perm> &S, std::map<int, Perm> &orbit, Full_stabs_chain &ans) {
    for (auto &g : S) {
        if (orbit.find(g[w]) == orbit.end()) {
            orbit[g[w]] = orbit[w] * g;


            if (ans.perm_n.find(orbit[g[w]]) == ans.perm_n.end()) {
                if (orbit[w] != Perm(ans.n))
                    ans.perm_n[orbit[g[w]]] = ans.perm_n[orbit[w]];
                if (g != Perm(ans.n))
                    ans.perm_n[orbit[g[w]]].insert(ans.perm_n[orbit[g[w]]].end(), ans.perm_n[g].begin(),
                                                   ans.perm_n[g].end());
            }

            build_schreier_tree(g[w], S, orbit, ans);
        }
    }
}


std::vector<int> inverse_gen(std::vector<int> p) {
    std::reverse(p.begin(), p.end());
    for (auto &i : p)
        i = -i;
    return p;
}


std::set<Perm> make_gen(std::set<Perm> &S, std::map<int, Perm> &orbit, Full_stabs_chain &ans) {
    std::set<Perm> NewS;
    for (auto &s: S) {
        for (auto &u : orbit) {
            Perm p = u.second * s;
            p *= inverse(orbit[s[u.first]]);

            if (ans.perm_n.find(p) == ans.perm_n.end()) {
                if (u.second != Perm(ans.n))
                    ans.perm_n[p].insert(ans.perm_n[p].end(), ans.perm_n[u.second].begin(), ans.perm_n[u.second].end());
                if (s != Perm(ans.n))
                    ans.perm_n[p].insert(ans.perm_n[p].end(), ans.perm_n[s].begin(), ans.perm_n[s].end());
                auto tmp = inverse_gen(ans.perm_n[orbit[s[u.first]]]);
                if (orbit[s[u.first]] != Perm(ans.n))
                    ans.perm_n[p].insert(ans.perm_n[p].end(), tmp.begin(), tmp.end());
            }


            NewS.insert(p);
        }
    }
    return NewS;
}


void schreier_sims(std::set<Perm> S, Full_stabs_chain &ans, int n) {
    std::map<int, Perm> orbit;
    Perm id(n);
    ans.n = n;
    ans.init_gen_set = S;
    make_map(S, ans);
    int w = 0;
    std::set<Perm> J = {id};
    while (S != J) {
        orbit.clear();
        orbit[w] = id;

        build_schreier_tree(w, S, orbit, ans);

        ans.schreier_trees.push_back(orbit);

        S = make_gen(S, orbit, ans);
        ans.stabs.push_back(S);
        for (const auto &i : orbit)
            ans.strong_set.insert(i.second);
        w++;

    }

    for (int i = 0; i < w; i++)
        ans.base.push_back(i);
}


bool in_group(Full_stabs_chain &chain, Perm g, std::vector<Perm> &ans) {
    int w = 0;
    Perm id(g.get_size());
    for (auto &i:chain.schreier_trees) {
        if (i.find(g[w]) == i.end())
            return false;

        if (i[g[w]] != id) {
            ans.push_back(i[g[w]]);
            g *= inverse(i[g[w]]);
        }

        w++;
    }
    return g == id;
}

void make_map(const std::set<Perm> &S, Full_stabs_chain &ans) {
    ans.perm_n[Perm(ans.n)].push_back(0);
    ans.n_perm[0] = Perm(ans.n);
    int c = 1;
    for (const auto &i : S) {
        ans.perm_n[i].push_back(c);
        ans.perm_n[inverse(i)].push_back(-c);
        ans.n_perm[c] = i;
        ans.n_perm[-c] = inverse(i);
        c++;
    }
}

size_t Full_stabs_chain::get_group_size() {
    size_t ans = 1;
    for (auto &i : schreier_trees)
        ans *= i.size();
    return ans;
}

std::set<Perm> Full_stabs_chain::get_gen_set_stab(int i) {
    return stabs[i];
}

std::vector<std::set<Perm>> Full_stabs_chain::get_stab_chain(int i) {
    std::vector<std::set<Perm>> ans;
    ans.insert(ans.end(), stabs.begin() + i, stabs.end());
    return ans;
}

std::vector<int> Full_stabs_chain::get_orbit(int i) {
    std::vector<int> ans;
    for (const auto &j: schreier_trees[i])
        ans.push_back(j.first);
    return ans;
}

bool check_full_stab_chain(Full_stabs_chain &chain) {
    std::map<int, Perm> orbit;
    Full_stabs_chain tmp;
    Perm id(chain.n);
    tmp.n = chain.n;
    make_map(chain.init_gen_set, tmp);

    for (int w = 0; w < chain.base.size(); w++) {
        orbit.clear();
        orbit[w] = id;
        if (w == 0) {
            build_schreier_tree(w, chain.init_gen_set, orbit, tmp);
        } else {
            build_schreier_tree(w, chain.stabs[w - 1], orbit, tmp);
        }
        if (orbit != chain.schreier_trees[w])
            return false;
    }
    return true;
}
