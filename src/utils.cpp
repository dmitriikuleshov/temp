#include "utils.h"

void save(const set_t &array, const std::string &filename) {
    std::ofstream fs(filename);
    fs << array.size() << std::endl;
    for (auto &n : array)
        n->save(fs);
    fs.flush();
    fs.close();
}

set_t load(const std::string &filename) {
    set_t result;
    std::ifstream is(filename);
    if (is.good() && is.is_open()) {
        int count;
        is >> count;
        for (int i = 0; i < count; ++i)
            result.insert(Factory::CreateNPC(is));
        is.close();
    } else
        std::cerr << "Error: " << std::strerror(errno) << std::endl;
    return result;
}

set_t fight(const set_t &array, size_t distance) {
    set_t dead_list;
    for (const auto &attacker : array)
        for (const auto &defender : array)
            if ((attacker != defender) &&
                (attacker->is_close(defender, distance))) {
                bool success = defender->accept(attacker);
                if (success)
                    dead_list.insert(defender);
            }
    return dead_list;
}

std::ostream &operator<<(std::ostream &os, const set_t &array) {
    for (auto &n : array)
        n->print();
    return os;
}