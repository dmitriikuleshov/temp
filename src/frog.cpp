#include "frog.h"

Frog::Frog(const std::string &name, int x, int y) : NPC(FrogType, name, x, y) {}
Frog::Frog(std::istream &is) : NPC(FrogType, is) {}

bool Frog::accept(const std::shared_ptr<NPC> &attacker) const {
    std::shared_ptr<Visitor> attacker_visitor =
        VisitorFactory::create_visitor(attacker->get_type());

    std::shared_ptr<NPC> defender =
        std::const_pointer_cast<NPC>(shared_from_this());

    bool result = attacker_visitor->visit(defender);
    attacker->fight_notify(defender, result);
    return result;
}

void Frog::print() { std::cout << *this; }

void Frog::save(std::ostream &os) {
    os << FrogType << std::endl;
    NPC::save(os);
}

std::ostream &operator<<(std::ostream &os, Frog &Frog) {
    os << "Frog: " << *static_cast<NPC *>(&Frog) << std::endl;
    return os;
}