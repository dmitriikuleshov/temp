#include <gtest/gtest.h>

#include "dragon.h"
#include "frog.h"
#include "knight.h"
#include "npc.h"
#include "visitor.h"

TEST(NPCTest, IsCloseTest) {
    auto knight = std::make_shared<Knight>("Knight1", 0, 14);
    auto frog = std::make_shared<Frog>("Frog1", 40, 11);
    auto dragon = std::make_shared<Dragon>("Dragon1", 88, 15);

    ASSERT_TRUE(knight->is_close(frog, 50));
    ASSERT_TRUE(frog->is_close(dragon, 100));
    ASSERT_FALSE(dragon->is_close(knight, 25));
}

TEST(NPCTest, AcceptTest) {
    auto knight = std::make_shared<Knight>("Knight1", 0, 14);
    auto frog = std::make_shared<Frog>("Frog1", 40, 11);
    auto dragon = std::make_shared<Dragon>("Dragon1", 88, 15);

    // FROG
    ASSERT_TRUE(knight->accept(frog));
    ASSERT_TRUE(dragon->accept(frog));
    ASSERT_TRUE(frog->accept(frog));

    // DRAGON
    ASSERT_TRUE(knight->accept(dragon));
    ASSERT_FALSE(dragon->accept(dragon));
    ASSERT_FALSE(frog->accept(dragon));

    // KNIGHT
    ASSERT_TRUE(dragon->accept(knight));
    ASSERT_FALSE(knight->accept(knight));
    ASSERT_FALSE(frog->accept(knight));
}

TEST(visitor_test, visit_test) {
    std::shared_ptr<NPC> knight, frog, dragon;
    knight = std::make_shared<Knight>("Knight1", 0, 14);
    frog = std::make_shared<Frog>("Frog1", 40, 11);
    dragon = std::make_shared<Dragon>("Dragon1", 88, 15);
    std::shared_ptr<Visitor> knight_visitor, frog_visitor, dragon_visitor;
    knight_visitor = std::make_shared<KnightVisitor>();
    frog_visitor = std::make_shared<FrogVisitor>();
    dragon_visitor = std::make_shared<DragonVisitor>();

    // FROG
    ASSERT_TRUE(frog_visitor->visit(knight));
    ASSERT_TRUE(frog_visitor->visit(dragon));
    ASSERT_TRUE(frog_visitor->visit(frog));

    // DRAGON
    ASSERT_TRUE(dragon_visitor->visit(knight));
    ASSERT_FALSE(dragon_visitor->visit(dragon));
    ASSERT_FALSE(dragon_visitor->visit(frog));

    // KNIGHT
    ASSERT_TRUE(knight_visitor->visit(dragon));
    ASSERT_FALSE(knight_visitor->visit(frog));
    ASSERT_FALSE(knight_visitor->visit(knight));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}