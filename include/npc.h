#pragma once
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <string>

// type for npcs
class NPC;
class Dragon;
class Knight;
class Frog;
using set_t = std::set<std::shared_ptr<NPC>>;

enum NpcType { Unknown = 0, DragonType = 1, KnightType = 2, FrogType = 3 };

struct IFightObserver {
    virtual void on_fight(const std::shared_ptr<NPC> attacker,
                          const std::shared_ptr<NPC> defender, bool win) = 0;
};

class NPC : public std::enable_shared_from_this<NPC> {
  protected:
    NpcType type;
    std::string name{""};
    int x{0};
    int y{0};

    std::vector<std::shared_ptr<IFightObserver>> observers;

  public:
    std::shared_ptr<NPC> get_shared() {
        return shared_from_this(); // Non-const version
    }
    NPC(NpcType t, const std::string &name, int _x, int _y);
    NPC(NpcType t, std::istream &is);

    NpcType get_type();

    virtual bool accept(const std::shared_ptr<NPC> &attacker) const = 0;

    void subscribe(const std::shared_ptr<IFightObserver> &observer);

    void fight_notify(const std::shared_ptr<NPC> defender, bool win) const;

    virtual bool is_close(const std::shared_ptr<NPC> &other,
                          size_t distance) const;

    virtual void print() = 0;

    virtual void save(std::ostream &os);

    std::pair<int, int> position() const;

    friend std::ostream &operator<<(std::ostream &os, NPC &npc);
};