#include <iostream>

#include <cassert>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <vector>

#include <algorithm>

#include <memory>
#include <stdexcept>

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
    SuffixTreeEdge *newEdge = new SuffixTreeEdge(start, end);
    newEdge->node = node;
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

SuffixTreeNode *SuffixTreeEdge::split(size_t *currentEnd, size_t length, char newChar,
                                      char differentChar) {
    if (length == 0 || length >= this->getLength()) {
        throw std::invalid_argument("length is greater than the substring's length");
    }

    SuffixTreeNode *internalNode = new SuffixTreeNode(newChar, *currentEnd - 1, currentEnd);
    internalNode->next[newChar]->node->id = SuffixTreeNode::currentMaxSuffixIndex++;
    internalNode->addEdge(this->node, differentChar, this->start + length, this->end);

    this->node = internalNode;
    this->end = new size_t(this->start + length);

    return internalNode;
}

size_t SuffixTreeEdge::getLength() { return *this->end - this->start; }

SuffixTreeEdge::~SuffixTreeEdge() { delete node; }

std::ostream &operator<<(std::ostream &os, const SuffixTreeEdge &edge) {
    os << ' ' << edge.start << ":" << *edge.end << " (" << edge.end << ") node: " << edge.node;
    return os;
}

// Implementation of suffix tree
class SuffixTree {
  public:
    std::string inputString;
    size_t inputStringLength = 0;

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

  public:
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
        char nextRemainingChar = this->inputString[currentCharIndex - this->remainder + 1];

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
            nextRemainingChar = this->inputString[currentCharIndex - this->activeLength];
        } else {
            // We're near the end of the suffix list
            nextRemainingChar = this->inputString[currentCharIndex - this->remainder + 1];
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
                    this->activeNode->next[inputString[currentCharIndex - this->activeLength]];
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
        char currentChar = this->inputString[currentCharIndex];

        while (this->remainder > 0) {
            if (this->activeLength == 0) {
                // inserting new edge to the active node

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
                // set up suffix links between internal nodes.
                if (this->activeNode != this->root) { // Rule2
                    previousInternalNode = previousInternalNode->setSuffixLink(this->activeNode);
                }

                this->remainder--;

            } else { // splitting active edge

                char checkedChar = this->inputString[this->activeEdge->start + this->activeLength];
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

            ruleOneOrThree(currentCharIndex); // Rule1 or Rule3 (or none)
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
        for (size_t phase = 0; phase < this->inputStringLength; ++phase) {

            // PHASE INITIALIZATION:
            // Update the global end counter - this automatically extends all
            // leaf edges
            (*this->end)++;
            // We have one new suffix to process (the entire string up to
            // current phase)
            this->remainder++;

            // Current character being processed in this phase
            size_t currentCharIndex = phase;
            char currentChar = this->inputString[currentCharIndex];

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
                char checkedChar = this->inputString[this->activeEdge->start + this->activeLength];
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

  protected:
    SuffixTree() {}

  public:
    SuffixTree(const std::string &inputString) {
        this->inputString = inputString + '\0';
        this->inputStringLength = this->inputString.length();
        this->buildTree();
    }

    ~SuffixTree() {
        delete root;
        delete end;
    }
};

#include <functional>
#include <set>
#include <unordered_map>

class SuffixTreeWithLCS : public SuffixTree {
  public:
    SuffixTreeWithLCS(const std::string &s1, const std::string &s2) {
        this->inputString = s1 + "#" + s2 + "$";
        this->inputStringLength = this->inputString.length();
        this->buildTree();
    }

    void findLongestCommonSubstrings() {
        size_t len1 = this->inputString.find('#');
        size_t len2 = this->inputString.length() - len1 - 2; // -2 for '#' and '$'

        std::unordered_map<SuffixTreeNode *, std::pair<bool, bool>> flagsMap;
        size_t maxLength = 0;

        std::function<std::pair<bool, bool>(SuffixTreeNode *, size_t)> dfs1;
        dfs1 = [&](SuffixTreeNode *node, size_t currentLength) -> std::pair<bool, bool> {
            bool isS1 = false;
            bool isS2 = false;

            // Если это листовой узел
            if (node->empty()) {
                if (node->id < len1) {
                    isS1 = true;
                } else if (node->id > len1 && node->id <= len1 + len2 + 1) {
                    isS2 = true;
                }
                flagsMap[node] = {isS1, isS2};
                return {isS1, isS2};
            }

            for (auto &[ch, edge] : node->next) {
                size_t edgeLength = edge->getLength();
                auto [childS1, childS2] = dfs1(edge->node, currentLength + edgeLength);
                isS1 = isS1 || childS1;
                isS2 = isS2 || childS2;
            }

            flagsMap[node] = {isS1, isS2};
            if (isS1 && isS2 && currentLength > maxLength) {
                maxLength = currentLength;
            }
            return {isS1, isS2};
        };

        dfs1(this->root, 0);

        std::set<std::string> results;
        if (maxLength > 0) {
            std::function<void(SuffixTreeNode *, size_t, std::string)> dfs2;
            dfs2 = [&](SuffixTreeNode *node, size_t currentLength, std::string currentString) {
                if (currentLength == maxLength) {
                    if (flagsMap[node].first && flagsMap[node].second) {
                        results.insert(currentString);
                    }
                    return;
                }

                for (auto &[ch, edge] : node->next) {
                    size_t edgeLength = edge->getLength();
                    if (currentLength + edgeLength > maxLength) {
                        size_t take = maxLength - currentLength;
                        std::string newString =
                            currentString + this->inputString.substr(edge->start, take);
                        if (flagsMap[edge->node].first && flagsMap[edge->node].second) {
                            results.insert(newString);
                        }
                    } else {
                        std::string newString =
                            currentString + this->inputString.substr(edge->start, edgeLength);
                        dfs2(edge->node, currentLength + edgeLength, newString);
                    }
                }
            };
            dfs2(this->root, 0, "");
        }

        std::cout << maxLength << std::endl;
        for (const auto &str : results) {
            std::cout << str << std::endl;
        }
    }
};

int main() {
    std::string s1, s2;
    std::cin >> s1 >> s2;

    SuffixTreeWithLCS tree(s1, s2);
    tree.findLongestCommonSubstrings();

    return 0;
}
