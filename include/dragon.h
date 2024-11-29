#pragma once
#include "npc.h"
#include "visitor.h"

struct Dragon : public NPC {
    Dragon(const std::string &name, int x, int y);
    Dragon(std::istream &is);

    void print() override;

    bool accept(const std::shared_ptr<NPC> &attacker) const override;

    void save(std::ostream &os) override;

    friend std::ostream &operator<<(std::ostream &os, Dragon &dragon);
};