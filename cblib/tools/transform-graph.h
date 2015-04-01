// A C++ program to find strongly connected components in a given
// directed graph using Tarjan's algorithm (single DFS)
//
// Original author: Anurag Singh
// Original source: http://www.geeksforgeeks.org
//
#include <list>
#include <stack>
#include <deque>
 
// A class that represents an directed graph
class Graph
{
    long long int V;    // No. of vertices
    std::list<long long int> *adj;    // A dynamic array of adjacency lists
 
    // A Recursive DFS based function used by SCC()
    void SCCUtil(long long int u, long long int disc[], long long int low[],
                 std::stack<long long int> *st, bool stackMember[], std::deque<long long int> *scc);
public:
    Graph(long long int V);   // Constructor
    ~Graph();
    void addArc(long long int v, long long int w);   // function to add an edge to graph
    void SCC(long long int from, std::deque<long long int> *scc);    // fetches strongly connected components
};