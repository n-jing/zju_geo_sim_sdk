#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>

#include "heap.h"

using namespace std;

template <typename T>
ostream &operator << (ostream &os, const vector<T> &v)
{
	copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	os << endl;
	return os;
}


template <typename T>
inline bool is_heap(const vector<T> &v)
{
	return hj::algorithm::is_heap(v.begin(), v.end());
}

void simple_test(vector<int> v)
{
	cout << "# simple_test: " << endl;
	hj::algorithm::make_heap(v.begin(), v.end());
	cout << v << is_heap(v) << "==1" << endl;

	hj::algorithm::pop_heap(v.begin(), v.end());
	cout << v << hj::algorithm::is_heap(v.begin(), v.end()-1) << "==1" << endl;
	hj::algorithm::push_heap(v.begin(), v.end());
	cout << v << is_heap(v) << "==1" << endl;
}

template <typename Container>
class less_value_address
{
public:
	less_value_address(const Container &dat)
		:_dat(dat) {
	}
	template <typename T>
	bool operator()(T a, T b) const {
		return _dat[a] < _dat[b];
	}
private:
	const Container &_dat;
};

template <typename OS, typename Value, typename Value_Address>
void print_heap_value_address(OS &os, const Value &value, const Value_Address &value_address)
{
	for(size_t i = 0; i < value_address.size(); ++i) {
		os << value[value_address[i]] << " ";
	}
	os << endl;
}

void test_heap_value_address(const vector<int> &v)
{
	cout << "# heap_value_address: " << endl;
	vector<int> value_address(v.size());
	for(size_t i = 0; i < v.size(); ++i) {
		value_address[i] = i;
	}

	print_heap_value_address(cout, v, value_address);
	cout << hj::algorithm::is_heap(value_address.begin(), value_address.end(), value_address.begin(),
                                 less_value_address<vector<int> >(v)) << "==0" << endl;

	hj::algorithm::make_heap(value_address.begin(), value_address.end(),
							 less_value_address<vector<int> >(v), hj::algorithm::default_swap_functor());
	print_heap_value_address(cout, v, value_address);
	cout << hj::algorithm::is_heap(value_address.begin(), value_address.end(), value_address.begin(), less_value_address<vector<int> >(v))
       << "==1" << endl;

	hj::algorithm::pop_heap(value_address.begin(), value_address.end(),
							less_value_address<vector<int> >(v), hj::algorithm::default_swap_functor());;
	print_heap_value_address(cout, v, value_address);
	cout << hj::algorithm::is_heap(value_address.begin(), value_address.end()-1, value_address.begin(),
                                 less_value_address<vector<int> >(v))
       << "==1" << endl;

	hj::algorithm::push_heap(value_address.begin(), value_address.end(),
							 less_value_address<vector<int> >(v), hj::algorithm::default_swap_functor());;
	print_heap_value_address(cout, v, value_address);
	cout << hj::algorithm::is_heap(value_address.begin(), value_address.end(), value_address.begin(),
                                 less_value_address<vector<int> >(v))
       << "==1" << endl;
}

template <typename Container>
class swap_dual_adress_0
{
public:
	swap_dual_adress_0(Container &heap_address)
		:_heap_address(heap_address) {
	}
	template <typename T>
	inline void operator()(T &a, T &b) const {
		std::swap(_heap_address[a], _heap_address[b]);
		std::swap(a, b);
	}
private:
	Container &_heap_address;
};

bool verify_heap_address(const vector<int> &v, const vector<int> &value_address, const vector<int> &heap_address)
{
	for(size_t i = 0; i < v.size(); ++i) {
		if(i != value_address[heap_address[i]])
			return false;
	}
	return true;
}

void test_heap_dual_address_0(const vector<int> &v)
{
	cout << "# heap_dual_address_0: " << endl;
	vector<int> value_address(v.size()), heap_address(v.size());
	for(size_t i = 0; i < v.size(); ++i) {
		value_address[i] = i;
		heap_address[i] = i;
	}

	hj::algorithm::make_heap(value_address.begin(), value_address.end(),
							 less_value_address<vector<int> >(v),
							 swap_dual_adress_0<vector<int> >(heap_address));
	print_heap_value_address(cout, v, value_address);
	cout << hj::algorithm::is_heap(value_address.begin(), value_address.end(), value_address.begin(), less_value_address<vector<int> >(v)) << endl;
	cout << verify_heap_address(v, value_address, heap_address) << endl;

	hj::algorithm::pop_heap(value_address.begin(), value_address.end(),
							less_value_address<vector<int> >(v),
							swap_dual_adress_0<vector<int> >(heap_address));
	print_heap_value_address(cout, v, value_address);
	cout << hj::algorithm::is_heap(value_address.begin(), value_address.end()-1, value_address.begin(),
								   less_value_address<vector<int> >(v)) << endl;
	cout << verify_heap_address(v, value_address, heap_address) << endl;

	hj::algorithm::push_heap(value_address.begin(), value_address.end(),
							 less_value_address<vector<int> >(v),
							 swap_dual_adress_0<vector<int> >(heap_address));
	print_heap_value_address(cout, v, value_address);
	cout << hj::algorithm::is_heap(value_address.begin(), value_address.end(), value_address.begin(),
								   less_value_address<vector<int> >(v)) << endl;
	cout << verify_heap_address(v, value_address, heap_address) << endl;
}

typedef vector<pair<const int *, size_t> > dual_address; // value and heap address

template <typename OS, typename Value>
void print_heap(OS &os, const Value &value, const dual_address &address)
{
	for(size_t i = 0; i < address.size(); ++i)
		os << *address[i].first << " ";
	os << endl;
}

bool verify_heap_address(const vector<int> &v, const dual_address &address)
{
	for(size_t i = 0; i < v.size(); ++i) {
		if(&v[i] != address[address[i].second].first);
			return false;
	}
	return true;
}

class heap_dual_address_op
{
public:
	class less
	{
	public:
		less(){}
		bool operator()(const pair<const int *, size_t> &a,
						const pair<const int *, size_t> &b) const {
			return *a.first < *b.first;
		}
	};
	class swap
	{
	public:
		swap(dual_address &da, const int *first)
			:_da(da), _first(first) {
		}
		void operator()(pair<const int *, size_t> &a,
						pair<const int *, size_t> &b) const {
			
			std::swap(_da[a.first-_first].second, _da[b.first-_first].second);
			std::swap(a.first, b.first);
		}
	private:
		dual_address &_da;
		const int *_first;
	};
};

void test_heap_dual_address(const vector<int> &v)
{
	cout << "# heap_dual_address" << endl;
	dual_address address; // value and heap address;
	for(size_t i = 0; i < v.size(); ++i) {
		address[i] = make_pair(&v[i], i);
	}

	hj::algorithm::make_heap(address.begin(), address.end(),
							 heap_dual_address_op::less(),
							 heap_dual_address_op::swap(address, &v[0]));
	print_heap(cout, v, address);
	cout << hj::algorithm::is_heap(address.begin(), address.end(), address.begin(),
								   heap_dual_address_op::less()) << endl;
	cout << verify_heap_address(v, address) << endl;

	hj::algorithm::pop_heap(address.begin(), address.end(),
							 heap_dual_address_op::less(),
							heap_dual_address_op::swap(address, &v[0]));
	print_heap(cout, v, address);
	cout << hj::algorithm::is_heap(address.begin(), address.end()-1, address.begin(),
								   heap_dual_address_op::less()) << endl;
	cout << verify_heap_address(v, address) << endl;

	hj::algorithm::push_heap(address.begin(), address.end(),
							 heap_dual_address_op::less(),
							 heap_dual_address_op::swap(address, &v[0]));
	print_heap(cout, v, address);
	cout << hj::algorithm::is_heap(address.begin(), address.end(), address.begin(),
								   heap_dual_address_op::less()) << endl;
	cout << verify_heap_address(v, address) << endl;
}

template <typename RandomIterator>
class heap_dual_address
{
public:
	typedef pair<size_t, size_t> heap_entry;
	typedef vector<heap_entry> dual_address;

	heap_dual_address(RandomIterator first, RandomIterator last, size_t size)
		:_address(size), _less(first), _swap(_address) {
		for(size = 0; size < _address.size(); ++size, ++first)
			_address[size] = make_pair(size, size);
	}
	class less_functor
	{
	public:
		less_functor(RandomIterator first)
			:_first(first) {
		}
		bool operator()(const heap_entry &a, const heap_entry &b) const {
			return *(_first+a.first) < *(_first+b.first);
		}
		const RandomIterator _first;
	};
	less_functor _less;

	class swap_functor
	{
	public:
		swap_functor(dual_address &da)
			:_da(da) {
		}
		void operator()(heap_entry &a, heap_entry &b) const {
			std::swap(_da[a.first].second, _da[b.first].second);
			std::swap(a.first, b.first);
		}
	private:
		dual_address &_da;
	};
	swap_functor _swap;

	dual_address _address;

	template <typename OS>
	void print(OS &os) {
		for(size_t i = 0; i < _address.size(); ++i)
			os << *(_less._first+_address[i].first) << " ";
		os << endl;
	}
	bool verify_heap_address(void) {
		for(size_t i = 0; i < _address.size(); ++i)
			if(i != _address[_address[i].first].second)
				return false;
		return true;
	}
};

void test_heap_dual_address_oo(const vector<int> &v)
{
	cout << "# test_heap_dual_address_oo" << endl;
	{
		vector<int> v1 = v;
		hj::algorithm::heap<vector<int>::iterator> h(v1.begin(), v1.end());
		h.make();
		cout << v1 << h.is_valid() << endl;
		h.pop();
		cout << v1 << h.is_valid() << endl;
		h.push(80);
		cout << v1 << h.is_valid() << endl;
		v1[4] = 20;
		cout << v1;
		h.update(4);
		cout << v1 << h.is_valid() << endl;
    cout << "begin to pop." << endl;
    while(!h.empty()) {
      cout << h.top() << endl;
      h.pop();
    }
		cout << "# test_heap_dual_address_oo." << endl;
	}

	heap_dual_address<const int *> address_(&v[0], &v[v.size()-1], v.size());

	hj::algorithm::heap<heap_dual_address<const int *>::dual_address::iterator, int,
		heap_dual_address<const int *>::less_functor,
		heap_dual_address<const int *>::swap_functor> h(address_._address.begin(),
														address_._address.end(),
														address_._less,
														address_._swap);
	h.make();
	address_.print(cout);
	cout << hj::algorithm::is_heap(address_._address.begin(), address_._address.end(), address_._address.begin(),
								   address_._less) << endl;
	cout << address_.verify_heap_address() << endl;
}

template <typename OS, typename Values, typename RandomIterator, typename ValueType>
void print_heap_index(OS &os, const Values &v, const hj::algorithm::heap_index<RandomIterator, ValueType> &heap)
{
  os << "val : ";
	for(size_t i = 0; i < heap.size(); ++i)
		os << v[i] << " ";
	os << endl;
  os << "h[i]: ";
	for(size_t i = 0; i < heap.size(); ++i)
		os << heap[i] << " ";
	os << endl;
}

void test_heap_dual_address_oo_2(const vector<int> &v)
{
	vector<int> v1 = v;
	cout << "# test_heap_dual_address_oo_2" << endl;
	hj::algorithm::heap_index<vector<int>::const_iterator> h(v1.begin(), v1.end());
	h.make();
	print_heap_index(cout, v1, h);
	cout << h.is_valid() << "==1" << endl;
	cout << h.top() << " " << h[0] << "==" << v1[h.top()] << endl;

	h.pop();
	print_heap_index(cout, v1, h);
	cout << h.is_valid() << "==1" << endl;
	cout << h.top() << " " << h[0] << "==" << v1[h.top()] << endl;

	v1[6] = 20;
	h.update(6);
	print_heap_index(cout, v1, h);
	cout << h.is_valid() << endl;
  cout << "begin to pop" << endl;
  while(!h.empty()) {
    cout << h.top() << "\t" << v1[h.top()] << endl;
    h.pop();
  }
	cout << "# end test_heap_dual_address_oo_2" << endl;
}

int main(int argc, char *argv[])
{
	vector<int> v(15);
	for(size_t i = 0; i < v.size(); ++i)
		v[i] = rand() % 100;
	cout << v << is_heap(v) << endl;
	
	simple_test(v);
	test_heap_value_address(v);
	test_heap_dual_address_0(v);
	test_heap_dual_address_oo(v);
	test_heap_dual_address_oo_2(v);

	cout << "# success." << endl;
	return 0;
}
