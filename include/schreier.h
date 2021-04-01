//
// Created by alex on 04.05.2020.
//

#ifndef SS_SCHREIER_H
#define SS_SCHREIER_H

#include <vector>
#include <map>
#include <set>


class Perm {   /* Класс перестановок*/
public:
    Perm() = default;

    explicit Perm(std::vector<int> p);

    explicit Perm(size_t n);

    bool operator<(const Perm &p) const;

    bool operator==(const Perm &p) const;

    bool operator!=(const Perm &p) const;

    Perm &operator*=(const Perm &p);

    int &operator[](size_t pos);

    const int &operator[](size_t pos) const;

    explicit operator bool() const;

    friend Perm operator*(Perm p1, const Perm &p2);

    friend Perm inverse(const Perm &p);

    void print() const;

    [[nodiscard]] size_t get_size() const;

private:
    std::vector<int> _data;
    std::size_t _size{0};
};

struct Full_stabs_chain { /* Класс полная цепочка стабилизаторов */
    size_t get_group_size(); /* Возвращает порядок группы */

    std::set<Perm> get_gen_set_stab(int i); /* Возвращает множество образующих Stab_i */

    std::vector<std::set<Perm>> get_stab_chain(int i); /* Возвращает цепочку стабилизаторов для Stab_i */

    std::vector<int> get_orbit(int i); /* Возвращает орбиту элемента базы как множество */

    int n;
    std::vector<int> base; /* Элементы базы */
    std::vector<std::map<int, Perm>> schreier_trees; /* Деревья Шраера */
    std::vector<std::set<Perm>> stabs; /* Образзующие G_i*/
    std::set<Perm> strong_set; /* Сильное порождающее множество */
    std::map<Perm, std::vector<int>> perm_n; /* Выражение перестановок через исходные */
    std::map<int, Perm> n_perm; /* Пронумерованные исходные перествновки и обратные к ним*/
    std::set<Perm> init_gen_set; /* Изначальное порождающее множество */
};

std::vector<int> inverse_gen(std::vector<int> p);


void build_schreier_tree(int w, std::set<Perm> &S, std::map<int, Perm> &orbit, Full_stabs_chain &ans);

std::set<Perm> make_gen(std::set<Perm> &S, std::map<int, Perm> &orbit, Full_stabs_chain &ans);

void schreier_sims(std::set<Perm> S, Full_stabs_chain &ans, int n); /* Алгоритм Шраера-Симса */

bool in_group(Full_stabs_chain &chain, Perm g, std::vector<Perm> &ans); /* Проверка принадлежности g группе*/

void make_map(const std::set<Perm> &S, Full_stabs_chain &ans);

bool check_full_stab_chain(Full_stabs_chain& chain); /* Проверка полной цепочки стабилизаторов */


#endif //SS_SCHREIER_H
