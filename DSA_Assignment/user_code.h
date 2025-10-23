/*
Name: Suhas Kamath
Email: suhaskamath@iisc.ac.in
Student ID: 25945
Course: Introduction to Scalable Systems (DS221)
Assignment 1
*/

#ifndef USER_CODE_H
#define USER_CODE_H

// Feel free to include more library functions
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
#include <set>
#include <limits>
#include <algorithm>

using namespace std;

/*
Details about the project:
This project contains solutions to three distinct problems related to parcel management and routing in a delivery network.
Each problem is encapsulated in its own function, with clear input and output specifications.
The solutions utilize various data structures and algorithms, including hash maps, sets, tree construction, and graph traversal techniques.
The code is written in C++ and adheres to the C++17 standard.

How Github Copilot was used:

Github Copilot was used to assist in writing some parts of the code.
However, I reviewed and modified the code to ensure correctness and efficiency.
It was not blindly accepted.
The parts where Copilot was used are mentioned in the comments above the respective functions or code blocks.
The prompts used for Copilot are also included in the comments.
All block comments added above the function and class definitions were written by Github Copilot.
*/

/**
 * @brief Finds all duplicate parcel IDs and their minimum weights.
 * 
 * For each parcel ID that appears more than once, returns a vector {id, min_weight}.
 * 
 * @param parcels Vector of {id, weight} pairs.
 * @return Vector of {id, min_weight} for duplicate IDs.
 */
vector<vector<int> > question_one(const vector<vector<int> >& parcels)
{
    map<int,int> parcels_min_weight;
    set<int> dup_package_ids;
    vector<vector<int> > result;

    for (int i = 0; i < parcels.size(); i++)
    {
        //Check if the package id already exists
        //Prompt: Write code to check if the package id already exists in the map
        if (parcels_min_weight.find(parcels[i][0]) != parcels_min_weight.end())
        {
            //Add duplicate package id to the list 
            dup_package_ids.insert(parcels[i][0]);
            //Update the min weight
            parcels_min_weight[parcels[i][0]] = min(parcels_min_weight[parcels[i][0]], parcels[i][1]);
        }
        //If the package id does not exist, add it to the map
        else
            parcels_min_weight[parcels[i][0]] = parcels[i][1];
    }

    //Prompt: Iterate through the list of duplicate package ids and add the min weight to the result vector
    for (set<int>::iterator it = dup_package_ids.begin(); it != dup_package_ids.end(); ++it) {
        vector<int> row;
        row.push_back(*it);
        row.push_back(parcels_min_weight[*it]);
        result.push_back(row);
    }

    //Prompt: Sort the result vector based on the package id
    struct CmpId {
        bool operator()(const vector<int>& a, const vector<int>& b) const { return a[0] < b[0]; }
    };
    sort(result.begin(), result.end(), CmpId());
    
    return result;
}

/**
 * @brief Binary tree node used for Question 2.
 * 
 * Prompt used: Define a binary tree node structure.
 */
struct TreeNode {
    int val;
    TreeNode* left;
    TreeNode* right;
    TreeNode(int v) : val(v), left(NULL), right(NULL) {}
};

/**
 * @brief Builds a binary tree from preorder and inorder traversals.
 *
 * @param preorder Preorder traversal array.
 * @param preStart Start index in preorder (inclusive).
 * @param preEnd End index in preorder (inclusive).
 * @param inorder Inorder traversal array.
 * @param inStart Start index in inorder (inclusive).
 * @param inEnd End index in inorder (inclusive).
 * @param inorderIndexMap Map from node value to its index in inorder.
 * @return TreeNode* Pointer to the root of the constructed subtree.
 */
TreeNode* buildTree(
    const vector<int>& preorder, int preStart, int preEnd,
    const vector<int>& inorder, int inStart, int inEnd,
    map<int, int>& inorderIndexMap
) {
    if (preStart > preEnd || inStart > inEnd) return NULL;

    int rootVal = preorder[preStart];
    TreeNode* root = new TreeNode(rootVal);
    int inRootIndex = inorderIndexMap[rootVal];
    int numsLeft = inRootIndex - inStart;

    root->left = buildTree(preorder, preStart + 1, preStart + numsLeft,
                                    inorder, inStart, inRootIndex - 1, inorderIndexMap);
    root->right = buildTree(preorder, preStart + numsLeft + 1, preEnd,
                                     inorder, inRootIndex + 1, inEnd, inorderIndexMap);
    return root;
}

//  Global Variables for LCA 
const int MAXN = 100005;   // adjust depending on constraints
const int LOG = 20;        // enough for n <= 1e6
vector<int> depth;
vector<vector<int> > up;    // up[v][i] = 2^i-th ancestor of v
map<int, TreeNode*> nodeMap; // map value → node pointer
map<int, int> nodeId;        // map value → integer id
vector<int> idToVal;                   // reverse mapping
int curId = 0;

/**
 * @brief DFS to initialize depth and binary-lifting ancestor table for LCA.
 *
 * @param root Current tree node.
 * @param parent Parent id of the current node (or -1 for root).
 * @param d Current depth.
 */
void dfs(TreeNode* root, int parent, int d)
{
    if (!root) return;
    int v = nodeId[root->val];
    depth[v] = d;

    // immediate parent
    up[v][0] = parent;

    // binary lifting ancestors
    for (int i = 1; i < LOG; i++) {
        if (up[v][i-1] != -1)
            up[v][i] = up[up[v][i-1]][i-1];
        else
            up[v][i] = -1;
    }

    if (root->left) dfs(root->left, v, d + 1);
    if (root->right) dfs(root->right, v, d + 1);
}

/**
 * @brief Computes the Lowest Common Ancestor (LCA) of two nodes.
 *
 * @param a Node id (compact index) of the first node.
 * @param b Node id (compact index) of the second node.
 * @return int Compact id of the LCA node.
 */
int lca(int a, int b)
{
    if (depth[a] < depth[b]) swap(a, b);

    // Lift a up to same depth as b
    int diff = depth[a] - depth[b];
    for (int i = LOG-1; i >= 0; i--) {
        if ((diff >> i) & 1)
            a = up[a][i];
    }

    if (a == b) return a;

    // Lift both up until just before LCA
    for (int i = LOG-1; i >= 0; i--) {
        if (up[a][i] != -1 && up[a][i] != up[b][i]) {
            a = up[a][i];
            b = up[b][i];
        }
    }
    return up[a][0]; // parent is LCA
}

/**
 * @brief Answers Question 2: for each query, find the junction containing all parcels.
 *
 * @param preorder Preorder traversal of the tree nodes (by value).
 * @param inorder Inorder traversal of the tree nodes (by value).
 * @param leafParcels Parcels present at each leaf (in left-to-right leaf order).
 * @param query_list List of queries; each query is a list of parcel ids.
 * @return vector<int> A value (node id) answer for each query.
 */
vector<int> question_two(
    const vector<int>& preorder,
    const vector<int>& inorder,
    const vector<vector<int> >& leafParcels,
    const vector<vector<int> >& query_list
) {
    int n = preorder.size();
    if (n == 0) return vector<int>();

    // Prompt used: Map node value to inorder index
    map<int, int> inorderIndexMap;
    for (int i = 0; i < n; i++)
        inorderIndexMap[inorder[i]] = i;

    // Build tree
    TreeNode* root = buildTree(preorder, 0, n - 1, inorder, 0, n - 1, inorderIndexMap);

    // Assign integer ids to nodes for easier indexing
    nodeId.clear(); idToVal.clear(); curId = 0;
    // recursive helper to assign IDs
    struct AssignIdsHelper {
        static void run(TreeNode* node, map<int,int>& nodeId, vector<int>& idToVal, int& curId) {
            if (!node) return;
            nodeId[node->val] = curId++;
            idToVal.push_back(node->val);
            run(node->left, nodeId, idToVal, curId);
            run(node->right, nodeId, idToVal, curId);
        }
    };
    AssignIdsHelper::run(root, nodeId, idToVal, curId);

    // Initialise LCA structures
    depth.assign(n, 0);
    up.assign(n, vector<int>(LOG, -1));
    dfs(root, -1, 0);

    // Map each parcel → leaf node id
    map<int, int> parcelToLeaf;
    int leafIndex = 0;
    // recursive helper to map parcels
    struct MapParcelsHelper {
        static void run(TreeNode* node, const vector<vector<int> >& leafParcels, map<int,int>& parcelToLeaf, map<int,int>& nodeId, int& leafIndex) {
            if (!node) return;
            if (!node->left && !node->right) {
                // leaf node
                if (leafIndex < (int)leafParcels.size()) {
                    const vector<int>& parcels = leafParcels[leafIndex];
                    for (size_t i = 0; i < parcels.size(); ++i) {
                        parcelToLeaf[parcels[i]] = nodeId[node->val];
                    }
                }
                leafIndex++;
            }
            run(node->left, leafParcels, parcelToLeaf, nodeId, leafIndex);
            run(node->right, leafParcels, parcelToLeaf, nodeId, leafIndex);
        }
    };
    MapParcelsHelper::run(root, leafParcels, parcelToLeaf, nodeId, leafIndex);

    // Answer queries
    vector<int> result;
    for (size_t qi = 0; qi < query_list.size(); ++qi) {
        const vector<int>& query = query_list[qi];
        if (query.empty()) {
            result.push_back(-1);
            continue;
        }
        int current = parcelToLeaf[query[0]];
        for (size_t i = 1; i < query.size(); i++) {
            current = lca(current, parcelToLeaf.find(query[i])->second);
        }
        result.push_back(idToVal[current]);
    }
    return result;
}

/**
 * @brief Builds an adjacency list representation of the graph.
 * 
 * @param edges List of edges {u, v, w}.
 * @param n Number of nodes.
 * @return Adjacency list graph.
 * 
 * Prompt used: Write a function to build an adjacency list from edge list.
 */
vector<vector<pair<int, int> > > build_graph(const vector<vector<int> >& edges, int n)
{
    vector<vector<pair<int, int> > > graph(n + 1);
    for (size_t idx = 0; idx < edges.size(); ++idx) {
        const vector<int>& e = edges[idx];
        int u = e[0], v = e[1], w = e[2];
        graph[u].push_back(std::make_pair(v, w));
        graph[v].push_back(std::make_pair(u, w));
    }
    return graph;
}

/**
 * @brief Modified Dijkstra's algorithm to handle booster logic.
 * 
 * Computes shortest times from start to all nodes, considering whether the booster has been used.
 * 
 * @param start Starting node.
 * @param graph Adjacency list.
 * @param metro Set of metro cities (booster stations).
 * @param n Number of nodes.
 * @return dist[node][0]: min time to node without booster, dist[node][1]: min time after booster.
 * 
 * Prompt used: Implement Dijkstra's algorithm.
 * I modified the code generated by Copilot to include the booster logic.
 */
vector<vector<long long> > dijkstra
(
    int start,
    const vector<vector<pair<int, int> > >& graph,
    const set<int>& metro,
    int n
)
{
    vector<vector<long long> > dist(n + 1, vector<long long>(2, numeric_limits<long long>::max()));
    struct State { long long t; int u; int booster; };
    struct StateCmp {
        bool operator()(const State& a, const State& b) const { return a.t > b.t; }
    };
    priority_queue<State, vector<State>, StateCmp> pq;
    dist[start][0] = 0;
    pq.push((State){0, start, 0}); // (time, node, booster_used)

    while (!pq.empty()) {
        State top = pq.top(); pq.pop();
        long long cur_time = top.t; int u = top.u; int booster = top.booster;
        if (cur_time > dist[u][booster]) continue;

        // If at a metro city and booster not used, can use booster here
        if (booster == 0 && metro.find(u) != metro.end()) {
            if (dist[u][1] > cur_time) {
                dist[u][1] = cur_time;
                pq.push((State){cur_time, u, 1});
            }
        }

        for (size_t ei = 0; ei < graph[u].size(); ++ei) {
            int v = graph[u][ei].first;
            int w = graph[u][ei].second;
            long long next_time = cur_time + (booster ? w / 2 : w);
            if (dist[v][booster] > next_time) {
                dist[v][booster] = next_time;
                pq.push((State){next_time, v, booster});
            }
        }
    }
    return dist;
}

/**
 * @brief Finds the earliest possible meeting time for two trucks starting from node 1 and node n.
 * 
 * Each truck can use a booster at a metro city to halve all subsequent travel times (once only).
 * Returns -1 if no meeting is possible.
 * 
 * @param edges List of edges {u, v, w}.
 * @param metro_cities List of metro city node numbers.
 * @return Minimum meeting time, or -1 if not possible.
 * 
 * No Github Copilot was used for this function.
 */
long long question_three
(
    const vector<vector<int> >& edges,
    const vector<int>& metro_cities
)
{
    if (edges.empty()) return -1;
    int n = 0;
    for (size_t i = 0; i < edges.size(); ++i) n = max(n, max(edges[i][0], edges[i][1]));
    vector<vector<pair<int,int> > > graph = build_graph(edges, n);
    set<int> metro(metro_cities.begin(), metro_cities.end());

    // Truck 1: from node 1
    vector<vector<long long> > dist1 = dijkstra(1, graph, metro, n);
    // Truck 2: from node n
    vector<vector<long long> > dist2 = dijkstra(n, graph, metro, n);

    long long res = numeric_limits<long long>::max();
    for (int i = 1; i <= n; ++i) {
        long long t1 = std::min(dist1[i][0], dist1[i][1]);
        long long t2 = std::min(dist2[i][0], dist2[i][1]);
        if (t1 < numeric_limits<long long>::max() && t2 < numeric_limits<long long>::max()) {
            res = std::min(res, std::max(t1, t2));
        }
    }
    return (res == numeric_limits<long long>::max() ? -1 : res);
}

#endif // USER_CODE_H