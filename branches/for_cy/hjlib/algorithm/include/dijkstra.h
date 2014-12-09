#ifndef HJ_ALG_DIJKSTRA_H_
#define HJ_ALG_DIJKSTRA_H_

#include "heap.h"

namespace hj { namespace algorithm {

//! @brief: start from seed, propagate the distance until reach the
//! Voronoi border.
// class graph
// {
// public:
//   size_t node_num(void) const;
//   typedef vector<pair<size_t, real_t> >::const_iterator const_iterator;
//   const_iterator begin(size_t h) const;
//   const_iterator end(size_t h) const;
// };

template <typename GRAPH, typename RanCon>
void distance(const GRAPH &g, size_t seed, RanCon &dist)
{
  typedef typename RanCon::value_type value_t;
  heap_index<typename RanCon::iterator, value_t, std::greater<value_t> > heap(dist.begin(), dist.end());
  dist[seed] = 0;
  heap.make();
  std::vector<bool> visited(g.node_num(), false);
  while(!heap.empty()) {
    assert(heap.is_valid());
    const size_t cur = heap.top();
    visited[cur] = true;
    heap.pop();
    for(typename GRAPH::const_iterator i = g.begin(cur); i != g.end(cur); ++i) { // for each neighbor
      if(visited[i->first]) continue;
      const value_t alt = dist[cur] + i->second;
      if(alt < dist[i->first]) {
        dist[i->first] = alt;
        heap.update(i->first);
      }
    }
  }
}

}}

#endif
