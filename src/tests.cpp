#include "Message.hpp"
#include "MessageBuilder.hpp"
#include <chrono>
#include <gtest/gtest.h>
#include <thread>
#include <zmq.hpp>

class MessageTest : public ::testing::Test {
  protected:
    void SetUp() override {
        context = new zmq::context_t(1);
        sender = new zmq::socket_t(*context, ZMQ_PUSH);
        receiver = new zmq::socket_t(*context, ZMQ_PULL);
        sender->bind("tcp://*:5555");
        receiver->connect("tcp://localhost:5555");
    }

    void TearDown() override {
        delete sender;
        delete receiver;
        delete context;
    }

    zmq::context_t *context;
    zmq::socket_t *sender;
    zmq::socket_t *receiver;
};

TEST_F(MessageTest, DefaultConstructor) {
    Message msg;
    EXPECT_EQ(msg.messageType, MessageTypes::EMPTY);
    EXPECT_EQ(msg.senderId, -1);
    EXPECT_EQ(msg.recieverId, -1);
    EXPECT_EQ(msg.sizeOfBody, 0);
    EXPECT_EQ(msg.body, nullptr);
}

TEST_F(MessageTest, ParameterizedConstructor) {
    int data = 123;
    Message msg(MessageTypes::CREATE_REQUEST, 1, 2, sizeof(data), &data);
    EXPECT_EQ(msg.messageType, MessageTypes::CREATE_REQUEST);
    EXPECT_EQ(msg.senderId, 1);
    EXPECT_EQ(msg.recieverId, 2);
    EXPECT_EQ(msg.sizeOfBody, sizeof(data));
    EXPECT_NE(msg.body, nullptr);
    EXPECT_EQ(memcmp(msg.body, &data, sizeof(data)), 0);
}

TEST_F(MessageTest, MoveConstructor) {
    int data = 123;
    Message original(MessageTypes::CREATE_REQUEST, 1, 2, sizeof(data), &data);
    Message moved(std::move(original));
    EXPECT_EQ(moved.messageType, MessageTypes::CREATE_REQUEST);
    EXPECT_EQ(moved.senderId, 1);
    EXPECT_EQ(moved.recieverId, 2);
    EXPECT_EQ(moved.sizeOfBody, sizeof(data));
    EXPECT_NE(moved.body, nullptr);
    EXPECT_EQ(original.body, nullptr);
}

TEST_F(MessageTest, MoveAssignmentOperator) {
    int data = 123;
    Message original(MessageTypes::CREATE_REQUEST, 1, 2, sizeof(data), &data);
    Message moved;
    moved = std::move(original);
    EXPECT_EQ(moved.messageType, MessageTypes::CREATE_REQUEST);
    EXPECT_EQ(moved.senderId, 1);
    EXPECT_EQ(moved.recieverId, 2);
    EXPECT_EQ(moved.sizeOfBody, sizeof(data));
    EXPECT_NE(moved.body, nullptr);
    EXPECT_EQ(original.body, nullptr);
}

TEST_F(MessageTest, SendAndReceiveMessage) {
    int data = 123;
    Message sendMsg(MessageTypes::CREATE_REQUEST, 1, 2, sizeof(data), &data);
    sendMsg.sendMessage(*sender);

    Message recvMsg;
    recvMsg.receiveMessage(*receiver);
    EXPECT_EQ(recvMsg.messageType, MessageTypes::CREATE_REQUEST);
    EXPECT_EQ(recvMsg.senderId, 1);
    EXPECT_EQ(recvMsg.recieverId, 2);
    EXPECT_EQ(recvMsg.sizeOfBody, sizeof(data));
    EXPECT_NE(recvMsg.body, nullptr);
    EXPECT_EQ(memcmp(recvMsg.body, &data, sizeof(data)), 0);
}

TEST_F(MessageTest, ReceiveMessageWithTimeoutSuccess) {
    int data = 123;
    Message sendMsg(MessageTypes::CREATE_REQUEST, 1, 2, sizeof(data), &data);
    sendMsg.sendMessage(*sender);

    Message recvMsg;
    bool received =
        recvMsg.receiveMessage(*receiver, std::chrono::milliseconds(1000));
    EXPECT_TRUE(received);
    EXPECT_EQ(recvMsg.messageType, MessageTypes::CREATE_REQUEST);
    EXPECT_EQ(recvMsg.senderId, 1);
    EXPECT_EQ(recvMsg.recieverId, 2);
    EXPECT_EQ(recvMsg.sizeOfBody, sizeof(data));
    EXPECT_NE(recvMsg.body, nullptr);
    EXPECT_EQ(memcmp(recvMsg.body, &data, sizeof(data)), 0);
}

TEST_F(MessageTest, ReceiveMessageWithTimeoutFail) {
    Message recvMsg;
    bool received =
        recvMsg.receiveMessage(*receiver, std::chrono::milliseconds(100));
    EXPECT_FALSE(received);
}

// Фикстура для тестов
class MessageBuilderTestFixture : public ::testing::Test {
  protected:
    void SetUp() override {
        // Здесь можно добавить общую настройку для всех тестов
    }

    void TearDown() override {
        // Здесь можно добавить общую очистку для всех тестов
    }
};

// Ваши предыдущие тесты
TEST_F(MessageBuilderTestFixture, PreviousTest1) {
    // Пример вашего предыдущего теста
    Message message = MessageBuilder::buildTestMessage();
    EXPECT_EQ(message.messageType, MessageTypes::TEST);
    EXPECT_EQ(message.senderId, -1);
    EXPECT_EQ(message.recieverId, -1);
    EXPECT_EQ(message.sizeOfBody, strlen("first") + 1);
    EXPECT_STREQ((char *)message.body, "first");
}

TEST_F(MessageBuilderTestFixture, PreviousTest2) {
    // Пример вашего предыдущего теста
    int id = 1;
    Message message = MessageBuilder::buildExitMessage(id);
    EXPECT_EQ(message.messageType, MessageTypes::QUIT);
    EXPECT_EQ(message.senderId, id);
    EXPECT_EQ(message.recieverId, id);
    EXPECT_EQ(message.sizeOfBody, 0);
    EXPECT_EQ(message.body, nullptr);
}

// Новые тесты для MessageBuilder
TEST_F(MessageBuilderTestFixture, BuildCreateMessage) {
    int id = 1;
    int data[2] = {10, 20};
    std::istringstream input("10 20");
    std::cin.rdbuf(input.rdbuf());

    Message message = MessageBuilder::buildCreateMessage(id);

    EXPECT_EQ(message.messageType, MessageTypes::CREATE_REQUEST);
    EXPECT_EQ(message.senderId, id);
    EXPECT_EQ(message.recieverId, id);
    EXPECT_EQ(message.sizeOfBody, sizeof(data));
    EXPECT_EQ(memcmp(message.body, data, sizeof(data)), 0);
}

TEST_F(MessageBuilderTestFixture, BuildPingRequest) {
    int time = 100;
    int id = 1;

    Message message = MessageBuilder::buildPingRequest(time, id);

    EXPECT_EQ(message.messageType, MessageTypes::HEARTBIT_REQUEST);
    EXPECT_EQ(message.senderId, id);
    EXPECT_EQ(message.recieverId, -1);
    EXPECT_EQ(message.sizeOfBody, sizeof(int));
    EXPECT_EQ(*(int *)message.body, time);
}

TEST_F(MessageBuilderTestFixture, BuildPingMessage) {
    unsigned long long time = 123456789;
    int id = 1;

    Message message = MessageBuilder::buildPingMessage(time, id);

    EXPECT_EQ(message.messageType, MessageTypes::HEARTBIT_RESULT);
    EXPECT_EQ(message.senderId, id);
    EXPECT_EQ(message.recieverId, -1);
    EXPECT_EQ(message.sizeOfBody, sizeof(unsigned long long));
    EXPECT_EQ(*(unsigned long long *)message.body, time);
}

TEST_F(MessageBuilderTestFixture, BuildExecMessage) {
    int id = 1;
    std::vector<std::string> data = {"command1", "command2"};
    std::istringstream input("1 command1 command2");
    std::cin.rdbuf(input.rdbuf());

    Message message = MessageBuilder::buildExecMessage();

    EXPECT_EQ(message.messageType, MessageTypes::EXEC_REQUEST);
    EXPECT_EQ(message.senderId, -1);
    EXPECT_EQ(message.recieverId, id);
    EXPECT_GT(message.sizeOfBody, 0);

    MessageData deserializedData = MessageBuilder::deserialize(message.body);
    EXPECT_EQ(deserializedData.id, id);
    EXPECT_EQ(deserializedData.data[0], data[0]);
    EXPECT_EQ(deserializedData.data[1], data[1]);
}

TEST_F(MessageBuilderTestFixture, GetSize) {
    std::vector<std::string> data = {"test1", "test2"};
    int expectedSize = 3 * sizeof(int) + data[0].size() + data[1].size();

    int size = MessageBuilder::getSize(data);

    EXPECT_EQ(size, expectedSize);
}
