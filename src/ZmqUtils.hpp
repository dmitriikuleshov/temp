#pragma once

#include "zmq.hpp"
#include <string>

class ZmqUtils {
  public:
    static const int PORT_TO_BIND_FROM = 30000;

    inline static const char *INPUT_URL_TEMPLATE = "tcp://*:";
    inline static const char *OUTPUT_URL_TEMPLATE = "tcp://localhost:";
    // in-process communication between threads in running manage node
    inline static const char *MESSAGE_PROCESSOR_URL = "inproc://processor";
    inline static const char *MESSAGE_SENDER_URL = "inproc://sender";

    static std::string getInputAddress(int port) {
        return INPUT_URL_TEMPLATE + std::to_string(port);
    }

    static std::string getOutputAddress(int port) {
        return OUTPUT_URL_TEMPLATE + std::to_string(port);
    }

    static int occupyPort(zmq::socket_t &socket) {
        int currentPort = PORT_TO_BIND_FROM;
        while (true) {
            try {
                socket.bind(getInputAddress(currentPort));
                break;
            } catch (const zmq::error_t &error) {
                currentPort++;
            }
        }
        return currentPort;
    }
};
