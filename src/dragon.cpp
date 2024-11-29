#include "dragon.h"

Dragon::Dragon(const std::string &name, int x, int y)
    : NPC(DragonType, name, x, y) {}
Dragon::Dragon(std::istream &is) : NPC(DragonType, is) {}

bool Dragon::accept(const std::shared_ptr<NPC> &attacker) const {
    std::shared_ptr<Visitor> attacker_visitor =
        VisitorFactory::create_visitor(attacker->get_type());

    std::shared_ptr<NPC> defender =
        std::const_pointer_cast<NPC>(shared_from_this());

    bool result = attacker_visitor->visit(defender);
    attacker->fight_notify(defender, result);
    return result;
}

void Dragon::print() { std::cout << *this; }

void Dragon::save(std::ostream &os) {
    os << DragonType << std::endl;
    NPC::save(os);
}

std::ostream &operator<<(std::ostream &os, Dragon &dragon) {
    os << "Dragon: " << *static_cast<NPC *>(&dragon) << std::endl;
    return os;
}