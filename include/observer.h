#pragma once
#include "npc.h"

class TextObserver final : public IFightObserver {
  private:
    TextObserver() {};

  public:
    static std::shared_ptr<IFightObserver> get() {
        static TextObserver instance;
        return std::shared_ptr<IFightObserver>(&instance,
                                               [](IFightObserver *) {});
    }

    void on_fight(const std::shared_ptr<NPC> attacker,
                  const std::shared_ptr<NPC> defender, bool win) override {
        if (win) {
            std::cout << std::endl << "Murder --------" << std::endl;
            std::cout << "killer: ";
            attacker->print();
            std::cout << "victim: ";
            defender->print();
        }
    }
};

class FileObserver final : public IFightObserver {
  private:
    FileObserver() {};

  public:
    static std::shared_ptr<IFightObserver> get() {
        static FileObserver instance;
        return std::shared_ptr<IFightObserver>(&instance,
                                               [](IFightObserver *) {});
    }

    void on_fight(const std::shared_ptr<NPC> attacker,
                  const std::shared_ptr<NPC> defender, bool win) override {
        if (win) {
            std::ofstream fs("log.txt", std::ios::app);
            fs << std::endl
               << "Murder --------" << std::endl
               << "killer: " << *attacker << std::endl
               << "victim: " << *defender;
            fs.close();
        }
    }
};