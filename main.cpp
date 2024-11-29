#include "utils.h"

int main() {
    set_t array;
    std::vector<std::string> names = {
        "Aeliana", "Thalion", "Gorath", "Kaelen", "Seraphina",
        "Darian",  "Eryndor", "Liriel", "Valtor", "Amara"};
    std::cout << "=== NPC Generator ===" << std::endl;
    std::cout << "Generating NPCs..." << std::endl;

    for (size_t i = 0; i < 10; ++i) {
        auto npc = Factory::CreateNPC(NpcType(std::rand() % 3 + 1), names[i],
                                      std::rand() % 100, std::rand() % 100);
        array.insert(npc);
        std::cout << "Created: " << *npc << std::endl;
    }
    std::cout << "\n\n";
    std::cout << "=== Saving NPCs to File ===" << std::endl;
    save(array, "npc.txt");
    std::cout << "=== Loading NPCs from File ===" << std::endl;
    array = load("npc.txt");
    std::cout << "=== Initial NPC List ===" << std::endl;
    std::cout << array;
    std::cout << "=== Starting Battles ===\n\n" << std::endl;

    for (size_t distance = 20; (distance <= 100) && !array.empty();
         distance += 10) {
        auto dead_list = fight(array, distance);

        for (auto &d : dead_list)
            array.erase(d);

        std::cout << "\n*** Fight Results ***" << std::endl;
        std::cout << "Distance: " << distance << " meters" << std::endl;
        std::cout << "Killed: " << dead_list.size() << " NPC(s)" << std::endl;
        std::cout << "Remaining: " << array.size() << " NPC(s)" << std::endl;
        std::cout << "**********************\n" << std::endl;
    }

    std::cout << "=== Final Survivors ===" << std::endl;
    if (!array.empty()) {
        std::cout << array;
    } else {
        std::cout << "No survivors remain." << std::endl;
    }

    return 0;
}
