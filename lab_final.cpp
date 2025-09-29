#include <iostream>

#include <algorithm>
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

struct SuffixTreeNode;
struct SuffixTreeEdge;

/**
 * Represents a node in the suffix tree.
 *
 * Each node contains:
 * - Edges to child nodes (keyed by starting characters)
 * - A suffix link for efficient traversal during construction
 * - An identifier for leaf nodes (represents the starting index of the suffix)
 */
struct SuffixTreeNode {
    // Static counter for assigning unique IDs to leaf nodes
    // Leaf node IDs represent the starting index of the suffix they contain
    static size_t currentMaxSuffixIndex;

    // Map of starting characters to outgoing edges
    // Each edge represents a substring starting with the key character
    std::map<char, SuffixTreeEdge *> next;

    // Pointer to another node that represents the same path minus the first
    // character Used for efficient suffix transitions during tree construction
    // (Rule 3)
    SuffixTreeNode *suffixLink = nullptr;

    // For leaf nodes: represents the starting index of the suffix in the input
    // string For internal nodes: typically 0 or unused
    size_t id = 0;

    // Default constructor - creates an empty node
    SuffixTreeNode();

    // Constructor that creates a node with a single edge
    // Used when creating new leaf edges during suffix insertion
    SuffixTreeNode(char value, size_t start, size_t *end);

    // Add a new leaf edge from this node
    // Creates a new edge starting with 'value' that leads to a new leaf node
    void addEdge(char value, size_t start, size_t *end);

    // Add an edge that connects to an existing node
    // Used during edge splitting to reconnect the remaining part of the edge
    void addEdge(SuffixTreeNode *node, char value, size_t start, size_t *end);

    // Set the suffix link for this node and return the target node
    // Chainable method used during suffix link maintenance (Rule 2)
    SuffixTreeNode *setSuffixLink(SuffixTreeNode *node);

    // Check if this node has any outgoing edges
    // Leaf nodes will be empty, internal nodes have at least 2 edges
    bool empty();

    // Clean up all edges and nodes recursively
    ~SuffixTreeNode();

    // Debug output operator for visualizing the node structure
    friend std::ostream &operator<<(std::ostream &os, const SuffixTreeNode &node);
};
size_t SuffixTreeNode::currentMaxSuffixIndex = 0;

/**
 * Represents an edge in the suffix tree.
 *
 * Each edge contains:
 * - Start and end indices defining the substring it represents
 * - A pointer to the node at the end of the edge
 * - For leaf edges: end pointer is shared and updated globally
 * - For internal edges: end pointer is fixed
 */
struct SuffixTreeEdge {
    // Starting index of the substring in the input string
    size_t start = 0;

    // Ending index (exclusive) of the substring
    // For leaf edges: points to a shared global counter that increases as tree
    // grows For internal edges: fixed value representing the actual end
    // position
    size_t *end = nullptr;

    // The node that this edge leads to
    // For leaf edges: points to an empty leaf node
    // For internal edges: points to a branching node
    SuffixTreeNode *node = nullptr;

    // Constructor creates an edge with the given start and end positions
    // Automatically creates a new node at the end of the edge
    SuffixTreeEdge(size_t start, size_t *end);

    SuffixTreeEdge(size_t start, size_t *end, SuffixTreeNode *node);

    // Split the edge at the specified position, creating a new internal node
    // length: how far along the edge to split (0 < length < edge length)
    // newChar: the character that caused the split (will form new leaf edge)
    // differentChar: the character that already exists on the edge
    // Returns the new internal node created by the split
    SuffixTreeNode *split(size_t *currentEnd, size_t length, char newChar, char differentChar);

    // Calculate the length of the substring represented by this edge
    // For leaf edges: length increases as the global end counter increases
    // For internal edges: fixed length
    size_t getLength();

    // Clean up the node if it exists
    ~SuffixTreeEdge();

    // Debug output operator for visualizing the edge
    friend std::ostream &operator<<(std::ostream &os, const SuffixTreeEdge &edge);
};

SuffixTreeNode::SuffixTreeNode() {}

SuffixTreeNode::SuffixTreeNode(char value, size_t start, size_t *end) {
    this->next[value] = new SuffixTreeEdge(start, end);
}

void SuffixTreeNode::addEdge(char value, size_t start, size_t *end) {
    this->next[value] = new SuffixTreeEdge(start, end);
    this->next[value]->node->id = SuffixTreeNode::currentMaxSuffixIndex++;
}

void SuffixTreeNode::addEdge(SuffixTreeNode *node, char value, size_t start, size_t *end) {
    SuffixTreeEdge *newEdge = new SuffixTreeEdge(start, end, node);
    this->next[value] = newEdge;
}

SuffixTreeNode *SuffixTreeNode::setSuffixLink(SuffixTreeNode *node) {
    this->suffixLink = node;
    return node;
}

bool SuffixTreeNode::empty() { return this->next.empty(); }

SuffixTreeNode::~SuffixTreeNode() {
    for (auto &[key, edge] : next) {
        delete edge;
    }
}

std::ostream &operator<<(std::ostream &os, const SuffixTreeNode &node) {
    os << "next:\n";
    for (const auto &[key, value] : node.next) {
        os << "\t" << key << " --- " << *value << '\n';
    }
    os << "suffixLink: " << node.suffixLink;
    return os;
}

SuffixTreeEdge::SuffixTreeEdge(size_t start, size_t *end)
    : start(start), end(end), node(new SuffixTreeNode()) {}

SuffixTreeEdge::SuffixTreeEdge(size_t start, size_t *end, SuffixTreeNode *node)
    : start(start), end(end), node(node) {}

SuffixTreeNode *SuffixTreeEdge::split(size_t *currentEnd, size_t length, char newChar,
                                      char differentChar) {
    // Length should not be greater than the substring's length
    assert(!(length == 0 || length >= this->getLength()));

    SuffixTreeNode *internalNode = new SuffixTreeNode(newChar, *currentEnd - 1, currentEnd);
    internalNode->next[newChar]->node->id = SuffixTreeNode::currentMaxSuffixIndex++;
    internalNode->addEdge(this->node, differentChar, this->start + length, this->end);

    this->node = internalNode;
    this->end = new size_t(this->start + length);

    return internalNode;
}

size_t SuffixTreeEdge::getLength() { return *this->end - this->start; }

SuffixTreeEdge::~SuffixTreeEdge() {
    if (!node->empty()) {
        delete end;
    }
    delete node;
}

std::ostream &operator<<(std::ostream &os, const SuffixTreeEdge &edge) {
    os << ' ' << edge.start << ":" << *edge.end << " (" << edge.end << ") node: " << edge.node;
    return os;
}

// Implementation of suffix tree
class SuffixTree {

  public:
    SuffixTree(const std::string &text) {
        this->text = text + '\0';
        this->textLength = this->text.length();
        this->buildTree();
    }

    ~SuffixTree() {
        delete root;
        delete end;
    }

  protected:
    std::string text;

    size_t textLength = 0;

    // Root of the suffix tree
    SuffixTreeNode *root = new SuffixTreeNode();

    // An active point (activeNode, activeEdge, activeLength)

    // The current node we're positioned at in the suffix tree.
    SuffixTreeNode *activeNode = this->root;

    // The current edge we're traversing from activeNode
    SuffixTreeEdge *activeEdge = nullptr;

    // How far we've progressed along the current active edge
    size_t activeLength = 0;

    // The length of the current suffix we're actively working with
    // in the current phase
    size_t depth = 0;

    // An integer indicating how many new suffixes we need to insert
    size_t remainder = 0;

    // Represents the current phase or how many characters
    // we've processed so far
    size_t *end = new size_t(0);

    // Default constructor, is used in descendant classes
    SuffixTree() {}

    /**
     * Rule1
     *
     * Applies whenever the activeNode is root
     *
     * - activeNode remains root;
     * - activeEdge is set to the first character of the
     *   new suffix we need to insert;
     * - activeLength is reduced by 1;
     *
     * Rule 1 is essentially saying: "When at root, to insert the next shorter
     * suffix, just drop the first character and continue from there." The depth
     * reduction represents this shortening of the suffix.
     */
    void ruleOne(size_t currentCharIndex) {

        assert(!(this->activeNode != this->root || this->depth == 0));

        // Calculate the first character of the next suffix we need to insert
        // remainder = number of suffixes still waiting to be inserted
        // currentCharIndex - remainder + 1 = start position of next suffix
        char nextRemainingChar = this->text[currentCharIndex - this->remainder + 1];

        // Apply the rule: reduce depth by 1 (we're moving to next suffix)

        // If the new activeLength is larger than the current edge's length,
        // we need canonizalization.
        this->depth--;

        // Try to find an edge from activeNode starting with nextRemainingChar
        const auto entry = this->activeNode->next.find(nextRemainingChar);
        if (entry != this->activeNode->next.end()) {
            // Edge exists - set active point to this edge
            this->activeEdge = entry->second;
            this->activeLength = this->depth; // New active length

            // Canonicalize to ensure active point is valid
            // (activeLength might be longer than current edge)
            this->canonicize(currentCharIndex);
        } else {
            // No edge exists - reset active point to node
            this->activeLength = 0;
            this->activeEdge = nullptr;
        }
    }

    /**
     * Rule3
     *
     * When we need to insert the next shorter suffix and we're not at the root,
     * Rule3 uses suffix links to efficiently jump to the next position instead
     * of restarting from the root.
     */
    void ruleThree(size_t currentCharIndex) {

        // Rule3 should only be applied when not a root
        assert(activeNode != this->root);

        char nextRemainingChar = '\0';

        if (this->remainder - 1 > this->activeLength) {
            // We have more suffixes remaining than our current edge position
            nextRemainingChar = this->text[currentCharIndex - this->activeLength];
        } else {
            // We're near the end of the suffix list
            nextRemainingChar = this->text[currentCharIndex - this->remainder + 1];
        }

        // Apply the rule

        if (this->activeNode->suffixLink != nullptr) {
            // A suffix link from node X (representing string "s") points to
            // node Y (representing string "s" without its first character)
            // This allows us to efficiently jump to the next shorter suffix
            this->activeNode = this->activeNode->suffixLink;
            this->depth--;
        } else {
            // No suffix link, trying to find next suffix starting from the root
            this->activeNode = this->root;
            this->ruleOne(currentCharIndex);
            return;
        }

        // Find the edge starting with our calculated character and set up the
        // new active point.
        const auto entry = this->activeNode->next.find(nextRemainingChar);
        if (entry != this->activeNode->next.end()) {
            this->activeEdge = entry->second;
            this->canonicize(currentCharIndex);
        } else {
            this->activeEdge = nullptr;
        }
    }

    /**
     * Use Rule1 or Rule3 depending on whether the activeNode is the root
     */
    void ruleOneOrThree(size_t currentCharIndex) {
        if (this->activeNode == this->root && this->activeLength != 0) { //! rule 1
            this->ruleOne(currentCharIndex);
        } else if (this->activeNode != this->root) { //! rule 3
            this->ruleThree(currentCharIndex);
        }
    }

    /**
     * Ensures that the active point (activeNode, activeEdge, activeLength) is
     * in its "canonical" form - meaning the active point should never point
     * beyond the current edge's length.
     *
     * It fixes situations where activeLength
     * has become larger than the current edge's length
     */
    void canonicize(size_t currentCharIndex) {
        size_t activeEdgeLength = this->activeEdge->getLength();

        while (this->activeLength >= activeEdgeLength) { // substring is too small, go further
                                                         // through the edges
            this->activeNode = this->activeEdge->node;
            this->activeLength = this->activeLength - activeEdgeLength;

            if (this->activeLength == 0) {
                this->activeEdge = nullptr;
                return;
            } else {
                this->activeEdge =
                    this->activeNode->next[text[currentCharIndex - this->activeLength]];
                activeEdgeLength = this->activeEdge->getLength();
            }
        }
    }

    /**
     * Processes all remaining suffixes (remainder > 0) for the
     * current phase, handling both edge creation and splitting while
     * maintaining suffix links
     *
     * previousInternalNode - Tracks the last created internal node for suffix
     * linking
     *
     * currentCharIndex - The current character position being processed
     */
    void buildSuffixLinks(SuffixTreeNode *&previousInternalNode, size_t currentCharIndex) {
        char currentChar = this->text[currentCharIndex];

        while (this->remainder > 0) {
            if (this->activeLength == 0) {
                // Inserting new edge to the active node

                const auto entry = this->activeNode->next.find(currentChar);
                if (entry != this->activeNode->next.end()) {
                    // Character already exists - extend matching
                    this->activeEdge = entry->second;
                    this->activeLength++;
                    this->depth++;
                    this->canonicize(currentCharIndex);
                    break;
                }
                // We found an edge starting with the current character
                this->activeNode->addEdge(currentChar, currentCharIndex, this->end);

                // Rule 2: When we create a new edge from a non-root node, we
                // Set up suffix links between internal nodes.
                if (this->activeNode != this->root) { // Rule2
                    previousInternalNode = previousInternalNode->setSuffixLink(this->activeNode);
                }

                this->remainder--;

            } else { // Splitting active edge

                char checkedChar = this->text[this->activeEdge->start + this->activeLength];
                if (checkedChar == currentChar) {
                    break; // Character matches - no split needed
                }

                // Create a new internal node and set up suffix links
                SuffixTreeNode *currentInternalNode = this->activeEdge->split(
                    this->end, this->activeLength, currentChar, checkedChar);
                previousInternalNode =
                    previousInternalNode->setSuffixLink(currentInternalNode); // Rule2
                this->remainder--;
            }

            // Rule1 or Rule3 (or none)
            ruleOneOrThree(currentCharIndex);
        }
    }

    /**
     * Main method that builds the suffix tree using Ukkonen's algorithm.
     *
     * Processes the input string character by character (each character is a
     * phase). In each phase, it adds the current character to all relevant
     * suffixes while maintaining the active point and suffix links for
     * efficiency.
     */
    void buildTree() {
        for (size_t phase = 0; phase < this->textLength; ++phase) {

            // PHASE INITIALIZATION:
            // Update the global end counter - this automatically extends all
            // leaf edges
            (*this->end)++;
            // We have one new suffix to process (the entire string up to
            // current phase)
            this->remainder++;

            // Current character being processed in this phase
            size_t currentCharIndex = phase;
            char currentChar = this->text[currentCharIndex];

            // CASE 1: Active point is exactly on a node (not in the middle of
            // an edge)
            if (this->activeLength == 0) {
                // Check if there's an existing edge from the current node
                // starting with currentChar
                const auto entry = this->activeNode->next.find(currentChar);
                if (entry != this->activeNode->next.end()) {
                    // SUB-CASE 1A: Edge exists - we can extend the match along
                    // this edge
                    this->activeEdge = entry->second;
                    this->activeLength++;
                    this->depth++;
                    this->canonicize(currentCharIndex);

                    // Note: We break here because we found a matching edge to
                    // extend.
                    // The remainder will be processed in subsequent phases or
                    // via suffix link
                } else {

                    // SUB-CASE 1B: No edge exists - create a new leaf edge
                    this->activeNode->addEdge(currentChar, phase, this->end);
                    this->remainder--; // One suffix has been processed

                    // If we're not at the root, we need to handle suffix links
                    if (this->activeNode != this->root) {
                        SuffixTreeNode *activeNodePreRule = this->activeNode;
                        // Move to the next shorter suffix using Rule 3 (since
                        // we're not at root)
                        this->ruleThree(currentCharIndex);
                        this->buildSuffixLinks(activeNodePreRule, currentCharIndex);
                    }

                    // If we're at the root, the remainder decrease above is
                    // sufficient
                    // and we'll continue with the next phase
                }

            } else {

                // CASE 2: Active point is in the middle of an edg (activeLength > 0)

                // Check the character at our current position on the active edge
                char checkedChar = this->text[this->activeEdge->start + this->activeLength];
                if (currentChar == checkedChar) {
                    // SUB-CASE 2A: Characters match - simply extend the match

                    this->activeLength++;
                    this->depth++;
                    this->canonicize(currentCharIndex);

                    // Continue to next phase - the match extension is sufficient

                } else {
                    // SUB-CASE 2B: Characters don't match - we need to split the edge

                    // Split the edge at the current position, creating a new internal node
                    SuffixTreeNode *previousInternalNode = this->activeEdge->split(
                        this->end, this->activeLength, currentChar, checkedChar);

                    this->remainder--; // One suffix has been processed (via splitting)

                    // Move to the next shorter suffix using appropriate rule
                    this->ruleOneOrThree(currentCharIndex);

                    // Process any remaining suffixes, maintaining suffix links
                    this->buildSuffixLinks(previousInternalNode, currentCharIndex);
                }
            }
        }
    }
};

class SuffixTreeWithLcs : public SuffixTree {
  private:
    std::string secondText;

  public:
    SuffixTreeWithLcs(const std::string &text1, const std::string &text2) : secondText(text2) {

        this->text = text1;

        // Build suffix tree using larger text
        if (text.length() > secondText.length()) {
            std::string temp = std::move(secondText);
            secondText = std::move(text);
            text = std::move(temp);
        }

        this->textLength = this->text.length();
        this->buildTree();
    }

    std::pair<std::vector<size_t>, size_t> findAllLcs() {
        const size_t patternLength = secondText.length();

        std::vector<size_t> indecies;
        size_t maxSubstringLength = 0;

        for (size_t suffixIndex = 0; suffixIndex < patternLength;
             ++suffixIndex) { // search for ith pattern's suffix
                              // via dfs from root
            this->activeNode = this->root;
            this->activeEdge = nullptr;
            this->activeLength = 0;

            size_t index;
            for (index = suffixIndex; index < patternLength; ++index) {
                const size_t substringLength = index - suffixIndex;
                const char currentChar = secondText[index];

                if (this->activeLength == 0) {
                    const auto entry = this->activeNode->next.find(currentChar);
                    if (entry != this->activeNode->next.end()) {
                        // There is an edge starting with this char
                        this->activeEdge = entry->second; // enter the edge
                        this->activeLength++;
                        if (this->activeLength == this->activeEdge->getLength()) {
                            this->activeNode = this->activeEdge->node;
                            this->activeEdge = nullptr;
                            this->activeLength = 0;
                        }
                    } else {
                        if (maxSubstringLength == substringLength) {
                            indecies.push_back(suffixIndex);
                        } else if (maxSubstringLength < substringLength) {
                            maxSubstringLength = substringLength;
                            indecies.clear();
                            indecies.push_back(suffixIndex);
                        }
                        break;
                    }
                } else {
                    const char checkedChar =
                        this->text[this->activeEdge->start + this->activeLength];
                    if (currentChar == checkedChar) {
                        // Char is in the edge's substring
                        this->activeLength++;
                        if (this->activeLength == this->activeEdge->getLength()) {
                            this->activeNode = this->activeEdge->node;
                            this->activeEdge = nullptr;
                            this->activeLength = 0;
                        }
                    } else {
                        // Char is not in the edge's substring, should
                        // Split and create a new internal node
                        if (maxSubstringLength == substringLength) {
                            indecies.push_back(suffixIndex);
                        } else if (maxSubstringLength < substringLength) {
                            maxSubstringLength = substringLength;
                            indecies.clear();
                            indecies.push_back(suffixIndex);
                        }
                        break;
                    }
                }
            }

            if (index == patternLength) {
                size_t substringLength = index - suffixIndex;
                if (maxSubstringLength == substringLength) {
                    indecies.push_back(suffixIndex);
                } else if (maxSubstringLength < substringLength) {
                    maxSubstringLength = substringLength;
                    indecies.clear();
                    indecies.push_back(suffixIndex);
                }
                break;
            }

            if (maxSubstringLength >= patternLength - suffixIndex) {
                // There are no strings, with length bigger
                // Than the current maxSubstringLength
                break;
            }
        }

        return std::make_pair(indecies, maxSubstringLength);
    }
};

// For leetcode task
class SuffixTreeWithScoreChecking : public SuffixTree {
  public:
    SuffixTreeWithScoreChecking(std::string text) : SuffixTree(text) {}

    size_t dfsLeavesCount(SuffixTreeNode *currentNode) {
        auto edgeMap = currentNode->next;
        if (edgeMap.empty())
            return 1;

        size_t leavesCount = 0;
        for (auto &[c, edge] : edgeMap) {
            leavesCount += dfsLeavesCount(edge->node);
        }
        return leavesCount;
    }

    size_t dfsLeavesCount(SuffixTreeNode *currentNode, size_t currnentDepth, size_t wordIndex) {
        auto edgeMap = currentNode->next;
        size_t leavesCount = 0;
        for (auto &[c, edge] : edgeMap) {
            if (edge->start == wordIndex)
                continue;

            leavesCount += dfsLeavesCount(edge->node);
        }
        return leavesCount;
    }

    size_t iterativeLeavesCount(SuffixTreeNode *currentNode, size_t currentDepth) {

        size_t leavesCount = 0;
        std::stack<SuffixTreeNode *> stack;

        auto &edgeMap = currentNode->next;

        for (auto &[c, edge] : edgeMap) {
            if (edge->start == currentDepth)
                continue;
            stack.push(edge->node);
        }

        while (!stack.empty()) {
            auto currentNode = stack.top();
            stack.pop();

            auto &edgeMap = currentNode->next;
            if (edgeMap.empty()) {
                leavesCount++;
                continue;
            }

            for (auto &[c, edge] : edgeMap) {
                stack.push(edge->node);
            }
        }
        return leavesCount;
    }

    size_t findScore() {
        size_t accumulatedScore = 0;

        // Start from the root node
        SuffixTreeNode *currentNode = this->root;
        size_t currentDepth = 0;

        while (currentDepth < textLength) {
            // Get the current character from word A
            char currentChar = text[currentDepth];

            // Check if current node has an edge starting with current character
            auto edgeIt = currentNode->next.find(currentChar);

            // Found matching edge - traverse it
            SuffixTreeEdge *currentEdge = edgeIt->second;

            // Compare the edge's substring with the remaining part of word A
            size_t edgeLength = currentEdge->getLength();
            size_t compareLength = std::min(edgeLength, textLength - currentDepth);

            // Update position in word A
            currentDepth += compareLength;

            // Move to the next node at the end of this edge
            currentNode = currentEdge->node;

            // If this is an internal node (has children), add current depth multiplied by edge
            // count

            // Use IterativeLeavesCount or dfsLeavesCount
            accumulatedScore += currentDepth * iterativeLeavesCount(currentNode, currentDepth);
        }
        // If we've processed all characters of word A, check if we're at a valid ending
        return accumulatedScore + textLength - 1;
    }
};

class Solution {
  public:
    long long sumScores(std::string s) {

        SuffixTreeWithScoreChecking tree(s);

        long long totalScore = tree.findScore();

        return totalScore;
    }
};

// int main() {
//     std::string s1, s2;
//     std::cin >> s1 >> s2;
//     auto tree = SuffixTreeWithLcs(s1, s2);
//     auto [startIndexes, length] = tree.findAllLcs();

//     if (length == 0) {
//         return 0;
//     }

//     std::set<std::string> substrings;

//     for (const size_t &start : startIndexes) {
//         substrings.insert(s2.substr(start, length));
//     }

//     std::cout << length << '\n';
//     for (const std::string &substring : substrings) {
//         std::cout << substring << '\n';
//     }
//     return 0;
// }
