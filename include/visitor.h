#pragma once
#include "npc.h"

class Visitor {
  public:
    virtual bool visit(const std::shared_ptr<NPC> &) const = 0;
};

class DragonVisitor final : public Visitor {
  public:
    bool visit(const std::shared_ptr<NPC> &npc) const override {
        bool result{};
        switch (npc->get_type()) {
        case NpcType::DragonType:
            result = false;
            break;

        case NpcType::FrogType:
            result = false;
            break;

        case NpcType::KnightType:
            result = true;
            break;
        }
        return result;
    }
};

class FrogVisitor final : public Visitor {
    bool visit(const std::shared_ptr<NPC> &npc) const override {
        bool result{};
        switch (npc->get_type()) {
        case NpcType::DragonType:
            result = true;
            break;

        case NpcType::FrogType:
            result = true;
            break;

        case NpcType::KnightType:
            result = true;
            break;
        }
        return result;
    }
};

class KnightVisitor final : public Visitor {
    bool visit(const std::shared_ptr<NPC> &npc) const override {
        bool result{};
        switch (npc->get_type()) {
        case NpcType::DragonType:
            result = true;
            break;

        case NpcType::FrogType:
            result = false;
            break;

        case NpcType::KnightType:
            result = false;
            break;
        }
        return result;
    }
};

class VisitorFactory {
  public:
    static std::shared_ptr<Visitor> create_visitor(const NpcType &type) {
        std::shared_ptr<Visitor> result;
        switch (type) {
        case NpcType::DragonType:
            result = std::static_pointer_cast<Visitor>(
                std::make_shared<DragonVisitor>());
            break;

        case NpcType::FrogType:
            result = std::static_pointer_cast<Visitor>(
                std::make_shared<FrogVisitor>());
            break;

        case NpcType::KnightType:
            result = std::static_pointer_cast<Visitor>(
                std::make_shared<KnightVisitor>());
            break;
        default:
            break;
        }
        return result;
    }
};
