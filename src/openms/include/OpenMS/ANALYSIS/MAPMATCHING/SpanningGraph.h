#pragma once
#include <stdio.h>
#include <stdlib.h>


using namespace OpenMS;
using namespace std;

class SpanningGraph
{  
    
public:
    int V;   
    vector<vector<int>> adj; 
    SpanningGraph(int V);   
    void addEdge(int v, int w);   
    bool containsCycle();  
    bool DFS(int v, vector<bool>& visited, int source); 
};
 
SpanningGraph::SpanningGraph(int V)
{
    this->V = V;
    vector<vector<int>> adj;
    for (Size i = 0; i < V; i++)
    {
      vector<int> a;
      adj.push_back(a);
    }
    this->adj = adj;
}
 
void SpanningGraph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

bool SpanningGraph::containsCycle()
{
  vector<bool> visited(V);
  for (int i = 0; i < V; i++)
  {
    visited[i] = false;
  }
    
  for (int u = 0; u < V; u++)
  {
    if (!visited[u])
    {
      if (DFS(u, visited, -1))
      {
        return true;
      }
    }
  } 
  return false;

}

//DFS helper
bool SpanningGraph::DFS(int v, vector<bool>& visited, int source)
{
    visited[v] = true;
    for (vector<int>::iterator it = adj[v].begin(); it != adj[v].end(); ++it)
    {
      if (!visited[*it])
      {
        if (DFS(*it, visited, v)) {return true;}
      }
      else if (*it != source) {return true;}
    }
    return false;
}

