cmake_minimum_required(VERSION 3.20)
project(lab6)

set(CMAKE_CXX_STANDARD 17)

find_package(cppzmq REQUIRED)

# Основные исполняемые файлы
add_executable(lab6 main.cpp AbstractNode.hpp ChildNodeInfo.hpp ZmqUtils.hpp Message.hpp MessageTypes.hpp MessageBuilder.hpp MessageData.hpp SpringBootApplication.hpp)
target_link_libraries(lab6 cppzmq pthread)

add_executable(server server.cpp AbstractNode.hpp ChildNodeInfo.hpp ZmqUtils.hpp Message.hpp MessageTypes.hpp MessageBuilder.hpp MessageData.hpp SpringBootApplication.hpp)
target_link_libraries(server cppzmq pthread)

# Добавляем тесты
add_executable(tests tests.cpp AbstractNode.hpp ChildNodeInfo.hpp ZmqUtils.hpp Message.hpp MessageTypes.hpp MessageBuilder.hpp MessageData.hpp SpringBootApplication.hpp)
target_link_libraries(tests cppzmq pthread gtest gtest_main)

# Добавляем тесты в CTest
add_test(NAME Lab6Tests COMMAND tests)