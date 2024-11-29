#pragma once
#include "factory.h"
#include "npc.h"

void save(const set_t &array, const std::string &filename);

set_t load(const std::string &filename);

set_t fight(const set_t &array, size_t distance);

std::ostream &operator<<(std::ostream &os, const set_t &array);