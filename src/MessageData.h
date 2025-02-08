#pragma once

#include <vector>
#include <string>

class MessageData
{
public:
    int id;
    int len1;
    int len2;
    std::vector<std::string> data;
};
