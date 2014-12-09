#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <set>
#include <vector>

#define BOOST_TEST_MODULE container

#include "../../include/container.h"
#include "../../include/operation.h"

using namespace std;
using namespace hj::half_edge;

template <template <typename> class Con>
int add_del()
{
  Con<int> c;
  set<int> s;

  vector<typename Con<int>::const_iterator> itrs;
  const size_t N = 10;
  for(size_t i = 0; i < N; ++i) {
    itrs.push_back(c.add(i));
    s.insert(i);
    BOOST_REQUIRE(itrs[i]);
  }

  BOOST_CHECK_EQUAL(c.size(), N);
  for(typename Con<int>::const_iterator i(c.begin());
      i != c.end(); ++i) {
    BOOST_REQUIRE(s.find(*i) != s.end());
  }

  const size_t idx[] = {0, N-1, N/2};
  vector<int> is_deleted(N, false);
  for(size_t j = 0; j < sizeof(idx)/sizeof(size_t); ++j) {
    s.erase(idx[j]);
    BOOST_REQUIRE(is_valid(c, itrs[idx[j]]));
    c.del(itrs[idx[j]]);
    BOOST_REQUIRE(!is_valid(c, itrs[idx[j]]));
    is_deleted[idx[j]] = true;

    BOOST_CHECK_EQUAL(c.size(), N-j-1);

    for(typename Con<int>::const_iterator i = c.begin();
        i != c.end(); ++i) {
      BOOST_REQUIRE(s.find(*i) != s.end());
    }

    for(size_t i = 0; i < itrs.size(); ++i) {
      if(is_deleted[i] == true) {
        for(typename Con<int>::const_iterator i1 = c.begin();
            i1 != c.end(); ++i1)
          BOOST_REQUIRE(i1 != itrs[i]);
      }
      else {
        BOOST_REQUIRE(*itrs[i] == i);
      }
    }
  }

  const typename Con<int>::const_iterator null;
  BOOST_REQUIRE(!null);

  return 0;
}

BOOST_AUTO_TEST_CASE(test_add_del)
{
  BOOST_CHECK_EQUAL(add_del<std_list>(), 0);
  BOOST_CHECK_EQUAL(add_del<std_vector>(), 0);
}
