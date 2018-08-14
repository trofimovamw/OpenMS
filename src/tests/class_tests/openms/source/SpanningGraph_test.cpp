
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/SpanningGraph.h>


using namespace OpenMS;
using namespace std;

START_TEST(SpanningGraph, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


int V = 4;
SpanningGraph graph = new SpanningGraph(V);
graph.addEdge(1,2);
graph.addEdge(1,3);
graph.addEdge(0,1);
graph.addEdge(2,3);


START_SECTION(void addEdge(int v, int w))
  SpanningGraph new_graph = graph;
  new_graph.addEdge(2,0);
  TEST_EQUAL(graph.adj[2],3);
END_SECTION


START_SECTION(void containsCycle())
  TEST_EQUAL(graph.containsCycle(),true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST