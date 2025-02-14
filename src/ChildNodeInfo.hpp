#pragma once

#include "MessageTypes.hpp"
#include <thread>

class ChildNodeInfo {

  public:
    ChildNodeInfo(pid_t pid, int messageRecieverPort, int childRegisterPort)
        : Pid(pid), ReceiverPort(messageRecieverPort),
          RegisterPort(childRegisterPort) {}
    pid_t Pid;
    int RegisterPort;

    int ReceiverPort;
};
