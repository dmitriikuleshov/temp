#pragma once
#include "dragon.h"
#include "frog.h"
#include "knight.h"
#include "npc.h"
#include "observer.h"

class Factory {
  public:
    static std::shared_ptr<NPC>
    CreateNPC(const NpcType &type, const std::string &name, int x, int y) {
        std::shared_ptr<NPC> result;
        switch (type) {
        case NpcType::DragonType:
            result = std::static_pointer_cast<NPC>(
                std::make_shared<Dragon>(name, x, y));
            break;

        case NpcType::FrogType:
            result = std::static_pointer_cast<NPC>(
                std::make_shared<Frog>(name, x, y));
            break;

        case NpcType::KnightType:
            result = std::static_pointer_cast<NPC>(
                std::make_shared<Knight>(name, x, y));
            break;
        default:
            break;
        }

        if (result) {
            result->subscribe(TextObserver::get());
            result->subscribe(FileObserver::get());
        }

        return result;
    }

    static std::shared_ptr<NPC> CreateNPC(std::istream &is) {
        std::shared_ptr<NPC> result;
        int type{0};
        if (is >> type) {
            switch (type) {
            case NpcType::DragonType:
                result =
                    std::static_pointer_cast<NPC>(std::make_shared<Dragon>(is));
                break;

            case NpcType::FrogType:
                result =
                    std::static_pointer_cast<NPC>(std::make_shared<Frog>(is));
                break;

            case NpcType::KnightType:
                result =
                    std::static_pointer_cast<NPC>(std::make_shared<Knight>(is));
                break;
            default:
                break;
            }
        } else {
            std::cerr << "Unexpected NPC type:" << type << std::endl;
        }

        if (result) {
            result->subscribe(TextObserver::get());
            result->subscribe(FileObserver::get());
        }

        return result;
    }
};
