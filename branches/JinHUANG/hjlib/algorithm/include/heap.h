#ifndef HJ_HEAP_H_
#define HJ_HEAP_H_

#include <cassert>
#include <vector>
#include <utility>

namespace hj { namespace algorithm {

/// generic programming style

class default_swap_functor
{
public:
	template <typename T>
	inline void operator()(T &a, T &b) const {
		std::swap(a, b);
	}
};

/// is_heap

template <typename RandomIterator, typename LessPred>
inline bool is_heap(RandomIterator first, RandomIterator last, RandomIterator root, const LessPred &less_pred)
{
	const RandomIterator left = first + (root-first)*2+1, right = left+1;
	if(left < last) { // has left child
		if(less_pred(*root, *left) || !is_heap(first, last, left, less_pred))
			return false;
		if(right < last) { // has right child
			if(less_pred(*root, *right) || !is_heap(first, last, right, less_pred))
				return false;
		}
	}
	// no child at all
	return true;
}

#ifndef __GXX_EXPERIMENTAL_CXX0X__
template <typename RandomIterator, typename LessPred>
inline bool is_heap(RandomIterator first, RandomIterator last, const LessPred &less_pred)
{
	return is_heap(first, last, first, less_pred);
}
#endif

template <typename RandomIterator>
inline bool is_heap(RandomIterator first, RandomIterator last)
{
	return is_heap(first, last, std::less<typename RandomIterator::value_type>());
}

/// update the heap between first and last with the value at root
template <typename RandomIterator, typename LessPred, typename SwapOp>
inline void update_heap(RandomIterator first, RandomIterator last, RandomIterator root,
						const LessPred &less_pred, const SwapOp &swap)
{
  assert(first <= last);
  //	assert(root >= first && root < last); // can update an empty heap
	{ // try up heap
		RandomIterator parent;
		while(1) {
			parent = first + (root-first-1)/2;
			if(parent >= first && less_pred(*parent, *root)) {
				swap(*parent, *root);
				assert(less_pred(*root, *parent));
				root = parent;
			}
			else
				break;
		}
	}
	{ // try down heap
		RandomIterator left, right, max_child;
		while(1) {
			left = first + (root-first)*2+1;
			if(left >= last) // no child
				break;
			right = left+1;
			if(right < last) // has two children
				max_child = less_pred(*left, *right)?right:left;
			else
				max_child = left;
			if(less_pred(*root, *max_child)) {
				swap(*root, *max_child);
				root = max_child;
				continue;
			}
			break;
		}
	}
}

template <typename RandomIterator, typename LessPred>
inline void update_heap(RandomIterator first, RandomIterator last, RandomIterator root,
						const LessPred &less_pred)
{
	update_heap(first, last, root, less_pred, default_swap_functor());
}

template <typename RandomIterator>
inline void update_heap(RandomIterator first, RandomIterator last, RandomIterator root)
{
	update_heap(first, last, root,
				std::less<typename RandomIterator::value_type>(), default_swap_functor());
}

/// pop_heap
template <typename RandomIterator, typename LessPred, typename SwapOp>
inline void pop_heap(RandomIterator first, RandomIterator last,
					 const LessPred &less_pred, const SwapOp &swap)
{
	swap(*first, *(last-1));
	update_heap(first, last-1, first, less_pred, swap);
}

template <typename RandomIterator, typename LessPred>
inline void pop_heap(RandomIterator first, RandomIterator last,
					 const LessPred &less_pred)
{
	pop_heap(first, last, less_pred, default_swap_functor());
}

template <typename RandomIterator>
inline void pop_heap(RandomIterator first, RandomIterator last)
{
	pop_heap(first, last, 
			 std::less<typename RandomIterator::value_type>(),
			 default_swap_functor());
}

/// push the element at last-1 into the heap
template <typename RandomIterator, typename LessPred, typename SwapOp>
inline void push_heap(RandomIterator first, RandomIterator last,
					  const LessPred &less_pred, const SwapOp &swap)
{
	if(first == last) return;
	update_heap(first, last, last-1, less_pred, swap);
}

template <typename RandomIterator, typename LessPred>
inline void push_heap(RandomIterator first, RandomIterator last,
					  const LessPred &less_pred)
{
	push_heap(first, last, less_pred, default_swap_functor());
}

template <typename RandomIterator>
inline void push_heap(RandomIterator first, RandomIterator last)
{
	push_heap(first, last,
			 std::less<typename RandomIterator::value_type>(),
			 default_swap_functor());
}

/// make_heap
template <typename RandomIterator, typename LessPred, typename SwapOp>
inline void make_heap(RandomIterator first, RandomIterator last,
					  const LessPred &less_pred, const SwapOp &swap)
{
	for(RandomIterator current = first; current <= last; ++current)
		push_heap(first, current, less_pred, swap);
}

template <typename RandomIterator, typename LessPred>
inline void make_heap(RandomIterator first, RandomIterator last,
					  const LessPred &less_pred)
{
	make_heap(first, last, less_pred, default_swap_functor());
}

template <typename RandomIterator>
inline void make_heap(RandomIterator first, RandomIterator last)
{
	make_heap(first, last,
			  std::less<typename RandomIterator::value_type>(),
			  default_swap_functor());
}

//! @brief a simple heap class provide the feature of update
template <typename RandomIterator,
		  typename ValueType = typename RandomIterator::value_type,
		  typename LessPred = std::less<ValueType>,
		  typename SwapOp = default_swap_functor>
class heap
{
public:
	heap(RandomIterator first, RandomIterator last,
		 const LessPred &less_pred, const SwapOp &swap)
		:_first(first), _last(last), _less_pred(less_pred), _swap(swap) {
	}
	heap(RandomIterator first, RandomIterator last,
		 const LessPred &less_pred)
		:_first(first), _last(last), _less_pred(less_pred), _swap(default_swap_functor()) {
	}
	heap(RandomIterator first, RandomIterator last)
		:_first(first), _last(last),
		 _less_pred(LessPred()), _swap(default_swap_functor()) {
	}
	bool is_valid(void) const {
		return is_heap(_first, _last, _less_pred);
	}
	void update(size_t pos) {
		update_heap(_first, _last, _first+pos, _less_pred, _swap);
	}
	void make(void) {
		make_heap(_first, _last, _less_pred, _swap);
	}
	const ValueType &top(void) const { // max according to LessPred
		return *_first;
	}
	void pop(void) {
		if(empty()) return;
		pop_heap(_first, _last, _less_pred, _swap);
		--_last;
	}
	// be carefull about exceed the range
	void push(const ValueType &v) {
		*_last++ = v;
		push_heap(_first, _last, _less_pred, _swap);
	}
	bool empty(void) const {
		return _first == _last;
	}
	size_t size(void) const {
		return _last - _first;
	}
private:
	RandomIterator _first, _last;
	const LessPred _less_pred;
	const SwapOp _swap;
};

//! @brief a index based heap without reordering the value.
template <typename RandomIterator,
          typename ValueType = typename RandomIterator::value_type,
          typename LessPred = std::less<ValueType> >
class heap_index
{
public:
	typedef heap_index<RandomIterator, ValueType> heap_type;

	heap_index(const RandomIterator first, const RandomIterator last)
		:_first(first), _last(last), _value_index(last-first), _heap_index(last-first),
		 _less(_first), _swap(_heap_index),
		 _heap(_value_index.begin(), _value_index.end(), _less, _swap) {
		for(size_t i = 0; i < _value_index.size(); ++i) {
			_value_index[i] = i;
			_heap_index[i] = i;
		}
	}
	const std::vector<size_t> &get_value_index(void) const {
		return _value_index;
	}
	const std::vector<size_t> &get_heap_index(void) const {
		return _heap_index;
	}

	const ValueType &operator[](size_t i) const {
		return *(_first+_value_index[i]);
	}

	bool is_valid(void) const {
		if(!_heap.is_valid())
			return false;
		for(size_t i = 0; i < _heap.size(); ++i) {
			if(i != _heap_index[_value_index[i]])
				return false;
		}
		return true;
	}
	void update(size_t pos) {
		_heap.update(_heap_index[pos]); // TODO
	}
	void make(void) {
		_heap.make();
	}
	const size_t &top(void) const {
		return _heap.top();
	}
	void pop(void) {
		_heap.pop();
	}
	// push is non-trival
	bool empty(void) const {
		return _heap.empty();
	}
	size_t size(void) const {
		return _heap.size();
	}
private:
	const RandomIterator _first, _last; // value iterator
	std::vector<size_t> _value_index; // point to range of RandomIterator
	std::vector<size_t> _heap_index; // point to _value_index

	class less_functor
	{
	public:
		less_functor(const RandomIterator &first)
			:_first(first), _less_pred(LessPred()) {
		}
		inline bool operator()(size_t a, size_t b) const {
      return _less_pred(*(_first+a), *(_first+b));
		}
	private:
		const RandomIterator &_first;
    const LessPred _less_pred;
	};

	class swap_functor
	{
	public:
		swap_functor(std::vector<size_t> &heap_index)
			:_heap_index(heap_index) {
		}
		inline void operator()(size_t &a, size_t &b) const {
			std::swap(_heap_index[a], _heap_index[b]);
			std::swap(a, b);
		}
	private:
		std::vector<size_t> &_heap_index;
	};

	const less_functor _less;
	const swap_functor _swap;

	heap<std::vector<size_t>::iterator, size_t, less_functor, swap_functor> _heap;
};

}} // namespace hj::algorithm

#endif
