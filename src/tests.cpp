#include "Message.hpp"
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

TEST_F(MessageTest, UpdateMessage) {
    Message msg(MessageTypes::CREATE_REQUEST, 1, 2);
    msg.update(3, 4);
    EXPECT_EQ(msg.senderId, 3);
    EXPECT_EQ(msg.recieverId, 4);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}