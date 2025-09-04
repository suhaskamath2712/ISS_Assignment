## **Function:** `question_one`  

---

Given a list of parcels, each represented as $[id, weight]$ , this function  
identifies parcels with duplicate IDs and weights and selects the one with the  
minimum weight for each such ID.  

$\textbf{Input:}$
- $\texttt{parcels}$ : vector of vectors where  
  - $\texttt{parcels[i][0]}=$ ID of the $i$ -th parcel  
  - $\texttt{parcels[i][1]}=$ weight of the $i$ -th parcel  


$\textbf{Input constraints:}$ 
- $1 \leq n \leq 10^{8}$, where $n$ is the number of parcels
- $0 \leq w \leq 10^{8}$, where $w$ is the weight of a parcel
- $0 \leq i \leq 10^{8}$, where $i$ is the id of a parcel

$\textbf{Output:}$
- $\texttt{output}$ : vector of vectors where each entry corresponds to a parcel  
  - $\texttt{output[i][0]}=$ parcel ID that appears more than once in the input  
  - $\texttt{output[i][1]}=$ minimum weight among all parcels with ID $\texttt{output[i][0]}$ 

**Function Signature:**  

```cpp
vector<vector<int>> question_one(const vector<vector<int>>& parcels);
```

<br><br>



## **Function:** `question_two`  

---

The warehouse conveyor belt network is structured as a binary tree:  
- Each node represents a junction where parcels may be routed left or right.  
...

$\textbf{Input:}$
1. **Two vectors representing tree traversals**  
   - $\texttt{preorder}$ : preorder traversal of the binary tree formed by the conveyor belt network  
   - $\texttt{inorder}$ : inorder traversal of the binary tree formed by the conveyor belt network  
   <!-- - These uniquely determine the binary tree structure.   -->

2. **A mapping of parcels at leaf nodes**  
   - $\texttt{leafParcels[i]}$ : vector of parcel IDs at the $i$ -th loading junction (leaf node)  
   - Parcels on internal nodes can be inferred from their subtrees.  

3. **A list of query lists of parcel IDs**  
   - $\texttt{query}$ : vector of vectors of parcel IDs to check  
   - $\texttt{query[i]}$ : $i$ -th vector of parcel IDs to check  


$\textbf{Input constraints:}$ 
- $1 \leq n \leq 10^{6}$, where $n$ is the number of loading junctions (leaf nodes)  
- $1 \leq m \leq 100$, where $m$ is the number of parcels on each loading junction  
- $1 \leq k \leq 10^{5}$, where $k$ is the number of queries  
- Junctions/nodes are numbered from 1 to $n$

$\textbf{Output:}$
- A vector of non-negative integers  

**Function Signature:**  

```cpp
vector<int> question_two(
    const vector<int>& preorder,
    const vector<int>& inorder,
    const vector<vector<int>>& leafParcels,
    const vector<vector<int>>& query
);
```



<br><br>


## **Function:** `question_three`  

---

The road map of the country is represented as a weighted graph:  
- Each city is a node.  
- Each road between two cities is an edge with a weight representing the time taken to travel the road with normal fuel.  
...

$\textbf{Input:}$
1. **Graph of the country**  
   - $\texttt{edges}$ : a vector of vectors where  
     - $\texttt{edges[i][0]}=$ starting city of the $i$ -th road (positive integer)  
     - $\texttt{edges[i][1]}=$ destination city of the $i$ -th road (positive integer)  
     - $\texttt{edges[i][2]}=$ time taken to travel from $\texttt{edges[i][0]}$ to $\texttt{edges[i][1]}$ with normal fuel  
     - Roads are bidirectional (undirected edges) and have even positive weights  

2. **Metro cities**  
   - $\texttt{metro-cities}$ : a vector containing metro cities (node numbers of metro cities)  

$\textbf{Input constraints:}$ 
- $2 \leq n \leq 2 \cdot 10^{5}$ , where $n$ is the number of cities (nodes in the graph)  
- $1 \leq m \leq 2 \cdot 10^{5}$ , where $m$ is the number of roads (edges in the graph)  
- $1 \leq k \leq n$, where $k$ is the number of metro cities 
 
$\textbf{Output:}$
- A non-negative integer or -1 if it is not possible for two friends to meet  

**Function Signature:**  

```cpp
long long question_three(
    const vector<vector<int>>& edges,
    const vector<int>& metro_cities
);

```




