#include "npc.h"

NPC::NPC(NpcType t, const std::string &_name, int _x, int _y)
    : type(t), name(_name), x(_x), y(_y) {}

NPC::NPC(NpcType t, std::istream &is) : type(t) {
    is >> name;
    is >> x;
    is >> y;
}

void NPC::subscribe(const std::shared_ptr<IFightObserver> &observer) {
    observers.push_back(observer);
}

void NPC::fight_notify(const std::shared_ptr<NPC> defender, bool win) const {
    for (auto &o : observers)
        o->on_fight(std::const_pointer_cast<NPC>(shared_from_this()), defender,
                    win);
}

bool NPC::is_close(const std::shared_ptr<NPC> &other, size_t distance) const {
    if ((std::pow(x - other->x, 2) + std::pow(y - other->y, 2)) <=
        std::pow(distance, 2))
        return true;
    else
        return false;
}

NpcType NPC::get_type() { return type; }

std::pair<int, int> NPC::position() const { return {x, y}; }

void NPC::save(std::ostream &os) {
    os << name << std::endl;
    os << x << std::endl;
    os << y << std::endl;
}

std::ostream &operator<<(std::ostream &os, NPC &npc) {
    os << npc.name << " "
       << "{ x:" << npc.x << ", y:" << npc.y << "} ";
    return os;
}
