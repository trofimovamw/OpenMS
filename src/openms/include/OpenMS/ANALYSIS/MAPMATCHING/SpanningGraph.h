#pragma once
#include <stdio.h>
#include <stdlib.h>


using namespace OpenMS;
using namespace std;

class SpanningGraph
{
    int V;   
    vector<int> *adj;
    bool isCyclic(int v, bool visited[], int source);    
    
public:
    SpanningGraph(int V);   
    void addEdge(int v, int w);   
    bool containsCycle();  
    
};
 
SpanningGraph::SpanningGraph(int V)
{
    this->V = V;
    adj = new vector<int>[V];
}
 
void SpanningGraph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}

bool SpanningGraph::containsCycle()
{
  bool *visited = new bool[V];
  for (int i = 0; i < V; i++)
  {
    visited[i] = false;
  }
    
  for (int u = 0; u < V; u++)
  {
    if (!visited[u])
    {
      if (isCyclic(u, visited, -1))
      {
        return true;
      }
    }
  } 
  return false;

}

bool SpanningGraph::isCyclic(int v, bool visited[], int source)
{
    visited[v] = true;
    vector<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
    {
      if (!visited[*i])
      {
        if (isCyclic(*i, visited, v)) {return true;}
      }
      else if (*i != source) {return true;}
    }
    return false;
}

