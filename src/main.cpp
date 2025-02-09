#include <iostream>

#include "ServerNode.hpp"
#include "SpringBootApplication.hpp"
#include <sys/wait.h>

void child(int sig) {
    pid_t pid;
    pid = wait(nullptr);
}

int main() {
    signal(SIGCHLD, child);
    SpringBootApplication app;
    app.run();
}
