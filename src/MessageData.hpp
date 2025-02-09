#pragma once

#include <string>
#include <vector>

class MessageData {
  public:
    int id;
    int len1;
    int len2;
    std::vector<std::string> data;
};
