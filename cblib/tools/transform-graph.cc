#include "transform-graph.h"
#include <algorithm>
#include <stdio.h>

using namespace std;

#define NIL -1


Graph::Graph(long long int V)
{
  this->V = V;
  adj = new list<long long int>[V];
}

Graph::~Graph()
{
  delete[] adj;
}

void Graph::addArc(long long int v, long long int w)
{
  // Skip if arc exists already
  if (std::find(adj[v].begin(), adj[v].end(), w) == adj[v].end()) {
    adj[v].push_back(w);
  }
}

// A recursive function that finds and prints strongly connected
// components using DFS traversal
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
//             discovery time) that can be reached from subtree
//             rooted with current vertex
// *st -- >> To store all the connected ancestors (could be part
//           of SCC)
// stackMember[] --> bit/index array for faster check whether
//                  a node is in stack
void Graph::SCCUtil(long long int u, long long int disc[], long long int low[], stack<long long int> *st,
                    bool stackMember[], deque<long long int> *scc)
{
  // A static variable is used for simplicity, we can avoid use
  // of static variable by passing a pointer.
  static long long int time = 0;

  // Initialize discovery time and low value
  disc[u] = low[u] = ++time;
  st->push(u);
  stackMember[u] = true;

  // Go through all vertices adjacent to this
  list<long long int>::iterator i;
  for (i = adj[u].begin(); i != adj[u].end(); ++i) {
    long long int v = *i;  // v is current adjacent of 'u'

    // If v is not visited yet, then recur for it
    if (disc[v] == -1) {
      SCCUtil(v, disc, low, st, stackMember, scc);

      // Check if the subtree rooted with 'v' has a
      // connection to one of the ancestors of 'u'
      // Case 1 (per above discussion on Disc and Low value)
      low[u]  = min(low[u], low[v]);
    }

    // Update low value of 'u' only of 'v' is still in stack
    // (i.e. it's a back edge, not cross edge).
    // Case 2 (per above discussion on Disc and Low value)
    else if (stackMember[v] == true)
      low[u]  = min(low[u], disc[v]);
  }

  // head node found, pop the stack and print an SCC
  long long int w = 0, wcnt = 0;  // To store stack extracted vertices and their count 
  if (low[u] == disc[u]) {
    if (st->top() == u) {
      
      // Do not store single elements
      w = (long long int) st->top();
      stackMember[w] = false;
      st->pop();

    } else {
      
      while (st->top() != u) {
        w = (long long int) st->top();
        scc->push_back(w);
        stackMember[w] = false;
        st->pop();
        --wcnt;
      }
      
      w = (long long int) st->top();
      scc->push_back(w);
      stackMember[w] = false;
      st->pop();
      --wcnt;

      // Strong component separator and size
      scc->push_back(wcnt);
    }
  }
}

// The function to do DFS traversal. It uses SCCUtil()
void Graph::SCC(long long int from, deque<long long int> *ssc)
{
  long long int *disc = new long long int[V];
  long long int *low = new long long int[V];
  bool *stackMember = new bool[V];
  stack<long long int> *st = new stack<long long int>();

  // Initialize disc and low, and stackMember arrays
  for (long long int i = 0; i < V; i++) {
    disc[i] = NIL;
    low[i] = NIL;
    stackMember[i] = false;
  }

  // Call the recursive helper function to find strongly
  // connected components in DFS tree with vertex 'i'
  //
  // TODO: Only consider integer variables (large gain possilbe) and 
  // diverging bounds (indicator of infeasibility)
  //
  for (long long int i = from; i < V; i++)
    if (disc[i] == NIL)
      SCCUtil(i, disc, low, st, stackMember, ssc);
      
  delete[] disc;
  delete[] low;
  delete[] stackMember;
  delete st;
}
