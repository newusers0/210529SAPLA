
#pragma once
#ifndef DOUBLY_LINKEDLIST_H
#define DOUBLY_LINKEDLIST_H
#include <iostream>
#include <initializer_list>
#include <cassert>

using namespace std;

template <typename T>
class doubly_linkedlist_iterator;

template <typename T>
class DoublyLinkedList;

/*=====================================================================================*/
template <typename T>
class DoublyListNode {

public:
	T _value;
	DoublyListNode<T> *_next = nullptr, *_prev = nullptr;
	/*..............................*/
	friend DoublyLinkedList<T>;
	friend doubly_linkedlist_iterator<T>;

public:
	DoublyListNode(const T& const) noexcept ;
	DoublyListNode() noexcept : _value(), _next(nullptr), _prev(nullptr){};//190828
	DoublyListNode<T>& operator =(const T&);
	DoublyListNode<T>& operator =(const DoublyListNode<T>&);
	DoublyListNode<T>& operator =(const DoublyListNode<T>*);

	template <typename T>
	friend std::ostream& operator <<(std::ostream&, const DoublyListNode<T>&);

	template <typename T>
	friend std::ostream& operator <<(std::ostream&, const DoublyLinkedList<T>&);

	//200301 return next node
	inline DoublyListNode<T>*& const next_node() { return this->_next; }

	//200301 return previous node
	inline DoublyListNode<T>*& const previous_node() { return this->_prev; }

	//200301 return next value
	inline T& const next_value() { return this->_next->_value; }

	//200301 return previous value
	inline T& const previous_value() { return this->_prev->_value; }
};
/*.......................................................................................*/

/*=========================================================================================*/
//190828 Helps "get" method, by saving last position
template<typename T>
struct LastNode
{
	// isCached should be set to FALSE everytime the list suffer changes
	bool isCached;

	// Helps "get" method, by saving last position
	int last_index;
	DoublyListNode<T> *last_node = nullptr;

	//LastNode(const bool& const isCached = false) : isCached(false){};//190828
	LastNode() noexcept : isCached(false), last_index(NULL), last_node(nullptr) {};//190828
	//190828
	~LastNode(){ 
		isCached = false;
		last_index = NULL;
		last_node = nullptr;
	};
	//LastNode(const bool& const isCached = false,const int& const last_index = NULL,const DoublyListNode<T>* last_node = nullptr) : isCached(isCached),last_index(last_index), last_node(last_node) {};//190828

};
/*................................................................................................*/

template <typename T>
class DoublyLinkedList {
private:
	DoublyListNode<T> *_head = nullptr, *_tail = nullptr;
	size_t _size;
	LastNode<T> _last_node;
public:
	friend class doubly_linkedlist_iterator<T>;

	//class doubly_linkedlist_iterator;//191113

	DoublyLinkedList() noexcept;
	~DoublyLinkedList();
	DoublyLinkedList(const std::initializer_list<T>&);
	size_t size() const noexcept;
	bool empty() const;
	bool contains(const T&) const;

	void resize(size_t container_size); //191121

	void push_back(const T&);// 191104 Appends a new element to the end of the container.  the same as emplace_back(T), add(T)
	void emplace_back(const T&);// 191104 Appends a new element to the end of the container.  the same as push_back(T), add(T)
	void emplace_front(const T&);//200301 insert value at head
	void add(const T& const);//191104 Appends a new element to the end of the container.  the same as push_back() & emplace_back()
	void add(size_t, const T&);
	void addAll(const std::initializer_list<T>&);
	//190903
	inline T& insertValueBeforeNode(const T& const, DoublyListNode<T>& const);
	//190912
	inline T& insertNodeBeforeNode(DoublyListNode<T>& const newElement, DoublyListNode<T>& const node);

	//200301 insert new node after current node
	inline T& insert_new_node_after_current_node(DoublyListNode<T>& const newElement, DoublyListNode<T>& const node);

	//200302 insert new value after current node
	inline T& insert_new_value_after_current_node(T& const value, DoublyListNode<T>& const node);

	bool remove(const T&);
	T remove(size_t);
	bool removeAll(const std::initializer_list<T>&);

	//190910
	inline T& removeNode(DoublyListNode<T>& const node);

	//200713
	inline T& pop_back();

	void clear();

	//190828 Get node of list
	DoublyListNode<T>& getNode(const size_t& const index) const;
	inline T& get(size_t) const;
	T set(size_t, const T&);
	size_t indexOf(const T&) const;
	size_t lastIndexOf(const T&) const;
	DoublyLinkedList<T>& subList(size_t, size_t) const;
	bool swap(size_t, size_t);

	//191029  
	DoublyLinkedList<T>& swap(DoublyLinkedList<T>& const list);

	//191029  
	DoublyLinkedList<T>& copy(const DoublyLinkedList<T>& const list);

	void sort();

	//DoublyListNode<T>& operator[](size_t) const;
	// add support to array brakets [] operator
	inline T& operator[](const int& const index);
	inline T& operator[](const size_t& const i) { return this->get(i); }
	inline const T& operator[]( const size_t& const i) const { return get(i); }

	template <typename T>
	friend std::ostream& operator << (std::ostream&, const DoublyLinkedList<T>&);
	/*191113
	// Return by reference so that it can be used in left hand side of the assignment expression
	*/
	inline DoublyListNode<T>*& const get_head_node() const noexcept { return this->_head; }
	/*191113
		// Return by reference so that it can be used in left hand side of the assignment expression
	*/
	inline DoublyListNode<T>*& get_tail_node() noexcept { return this->_tail; }
	/*190828
		get front/begin of linked list
	*/
	inline T& front() noexcept;
	/*191211
		get const front/begin of linked list
	*/
	inline const T& const front() const noexcept;
	/*190828
		get tail/end of linked list
	*/
	inline T& back() const noexcept;
	/*191113
		// Root of LinkedList wrapped in Iterator type
	*/
	inline doubly_linkedlist_iterator<T> begin() const noexcept { return doubly_linkedlist_iterator<T>(this->_head);}
	/*190828
		get tail/end of linked list, should return nullptr
	*/
	inline const doubly_linkedlist_iterator<T> end() const noexcept {return doubly_linkedlist_iterator<T>();}
	/*191113
		// Root of LinkedList wrapped in Iterator type
	*/
	inline doubly_linkedlist_iterator<T> cbegin() const noexcept { return doubly_linkedlist_iterator<T>(this->_head);}
	/*190828
		get tail/end of linked list
	*/
	inline const doubly_linkedlist_iterator<T> cend() const noexcept {return doubly_linkedlist_iterator<T>();}
	/*191113
		// Root of LinkedList wrapped in Iterator type
	*/
	inline doubly_linkedlist_iterator<T> rbegin() const noexcept { return doubly_linkedlist_iterator<T>(_tail);}
	/*190828
		get tail/end of linked list
	*/
	inline const doubly_linkedlist_iterator<T> rend() const noexcept {return doubly_linkedlist_iterator<T>();}
	/*191113
	traverse doubly linked list
	*/
	void traverse_linkedlist();
};

template <typename T>
class doubly_linkedlist_iterator {
private:
	const DoublyListNode<T>* current_node_pointer = nullptr;
public:
	
	//Iterator() noexcept : current_node_pointer(_head) { }
	doubly_linkedlist_iterator() noexcept : current_node_pointer(nullptr) { }

	doubly_linkedlist_iterator(const DoublyListNode<T>* node) noexcept : current_node_pointer(node)  { }

	doubly_linkedlist_iterator(const DoublyLinkedList<T>& const linked_list) noexcept : current_node_pointer(linked_list._head) { }

	doubly_linkedlist_iterator(const doubly_linkedlist_iterator<T>& const list_iterator) noexcept : current_node_pointer(list_iterator.current_node_pointer) { }

	//bool hasNext()
	inline doubly_linkedlist_iterator<T>& const operator = (const DoublyListNode<T>* const node) const;
	doubly_linkedlist_iterator<T>& operator++();
	doubly_linkedlist_iterator<T> operator++(int);
	doubly_linkedlist_iterator<T>& operator--();
	doubly_linkedlist_iterator<T> operator--(int);
	inline bool operator != (const doubly_linkedlist_iterator& const iterator) const;
	inline bool operator == (const doubly_linkedlist_iterator& const iterator) const;
	inline const T& operator * () noexcept;
	inline const T* operator -> () noexcept;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
================================= Realization of class DoublyListNode =================================
*/

template <typename T>
DoublyListNode<T>::DoublyListNode(const T& const value) noexcept : _value(value), _next(nullptr), _prev(nullptr){}

template <typename T>
std::ostream& operator << (std::ostream &out, const DoublyListNode<T> &element) {
	out << element._value;
	return out;
}

template <typename T>
DoublyListNode<T>& DoublyListNode<T>::operator = (const T &value) { _value = value; return *this; }

template <typename T>
DoublyListNode<T>& DoublyListNode<T>::operator = (const DoublyListNode<T> &element) { _value = element._value; return *this; }

template <typename T>
DoublyListNode<T>& DoublyListNode<T>::operator = (const DoublyListNode<T> *element) { _value = element->_value; return *this; }
/*..................................................................................................................*/



/*
================================= Realization of class DoublyLinkedList =================================
*/

template <typename T>
DoublyLinkedList<T>::DoublyLinkedList() noexcept : _head(nullptr), _tail(nullptr), _size(0) {}

template <typename T>
DoublyLinkedList<T>::~DoublyLinkedList() {
	clear();
	delete _head;
	delete _tail;
	_head = _tail = nullptr;

	_last_node.~LastNode();
}

template <typename T>
DoublyLinkedList<T>::DoublyLinkedList(const std::initializer_list<T> &elements) {
	for (T element : elements)
		add(element);
}

template <typename T>
size_t DoublyLinkedList<T>::size() const noexcept { return _size; }

template <typename T>
bool DoublyLinkedList<T>::empty() const { return _size == 0; }

template <typename T>
bool DoublyLinkedList<T>::contains(const T &value) const {
	DoublyListNode<T>* node = _head;
	while (node != nullptr) {
		if (node->_value != value)
			node = node->_next;
		else
			return true;
	}
	return false;
}

//************************************
// Method:resize
// Qualifier: Resizes the container so that it contains n elements.
// Input: 
// Output: 
// date:191121
//************************************
template <typename T>
void DoublyLinkedList<T>::resize(size_t container_size) { //191121
	T temp_element;
	for (int i = 0; i < container_size; i++) {
		add(temp_element);
	}
}


//************************************
// Method:push_back
// Qualifier: Appends a new element to the end of the container. the same as emplace_back(T), add(T)
// Input: value, node
// Output: refreshed doubly linked list
// date:191104
//************************************
template <typename T>
void DoublyLinkedList<T>::push_back(const T& value) {// 191104 Appends a new element to the end of the container.  the same as emplace_back(T), add(T)
	DoublyListNode<T>* newElement = new DoublyListNode<T>(value);
	if (_size != 0) {
		newElement->_prev = _tail;
		_tail->_next = newElement;
		_tail = newElement;
	}
	else {
		_head = _tail = newElement;
	}
	++_size;

	_last_node.isCached = false;
}

//************************************
// Method:emplace_back
// Qualifier: Appends a new element to the end of the container.the same as push_back(T), add(T)
// Input: value, node
// Output: refreshed doubly linked list
// date:191104
//************************************
template <typename T>
void DoublyLinkedList<T>::emplace_back(const T& value) {// 191104 push back value at the tail of linke list
	DoublyListNode<T>* newElement = new DoublyListNode<T>(value);
	if (_size != 0) {
		newElement->_prev = _tail;
		_tail->_next = newElement;
		_tail = newElement;
	}
	else {
		_head = _tail = newElement;
	}
	++_size;

	_last_node.isCached = false;
}

//************************************
// Method:emplace_front
// Qualifier: Appends a new element to the head of the container.the same as push_front(T), add(T)
// Input: value, node
// Output: refreshed doubly linked list
// date:200301
//************************************
template <typename T>
void DoublyLinkedList<T>::emplace_front(const T& value) {//200301 insert value at head
	DoublyListNode<T>* newElement = new DoublyListNode<T>(value);
	if (_size != 0) {
		newElement->_next = _head;
		_head->_prev = newElement;
		_head = newElement;
	}
	else {
		_head = _tail = newElement;
	}
	++_size;

	_last_node.isCached = false;
}

template <typename T>
void DoublyLinkedList<T>::add(const T& const value) {
	DoublyListNode<T>* newElement = new DoublyListNode<T>(value);
	if (_size != 0) {
		newElement->_prev = _tail;
		_tail->_next = newElement;
		_tail = newElement;
	}
	else {
		_head = _tail = newElement;
	}
	++_size;

	_last_node.isCached = false;
}

template <typename T>
void DoublyLinkedList<T>::add(size_t index, const T &value) {
	if (_size == 0) {
		add(value);
		return;
	}
	if (index >= 0 && index < _size) {
		DoublyListNode<T>* node = _head;
		for (size_t i = 0; i < index; ++i)
			node = node->_next;
		DoublyListNode<T>* newElement = new DoublyListNode<T>(value);
		newElement->_next = node;
		newElement->_prev = node->_prev;
		if (node->_prev != nullptr) node->_prev->_next = newElement;
		node->_prev = newElement;
		if (index == 0) _head = newElement;
		if (index == _size - 1) _tail = newElement->_next;
		++_size;
		_last_node.isCached = false;
	}
	else
		throw std::out_of_range("DoublyLinkedList :: add(index, value)");
}

template <typename T>
void DoublyLinkedList<T>::addAll(const std::initializer_list<T> &elements) {
	for (T element : elements)
		add(element);
	_last_node.isCached = false;
}

//************************************
// Method:insertPreviousNode
// Qualifier: Insert new node before specific node.
// Input: value, node
// Output: refreshed doubly linked list
// date:190903
//************************************
template <typename T>
inline T& DoublyLinkedList<T>::insertValueBeforeNode(const T& const value, DoublyListNode<T>& const node) {
	DoublyListNode<T>* newElement = new DoublyListNode<T>(value);
	return insertNodeBeforeNode(*newElement, node);
}

//************************************
// Method:insertPreviousNode
// Qualifier: Insert new node before specific node.
// Input: value, node
// Output: refreshed doubly linked list
// date:190912
//************************************
template <typename T>
inline T& DoublyLinkedList<T>::insertNodeBeforeNode(DoublyListNode<T>& const newElement, DoublyListNode<T>& const node) {
	if (_size == 0) {

		assert(0);
	}
	else if (_size == 1) {

	}

	newElement._next = &node;
	newElement._prev = node._prev;
	if (node._prev != nullptr) node._prev->_next = &newElement;
	node._prev = &newElement;
	if (_head == &node) _head = &newElement;
	if (_tail == &node) _tail = newElement._next;
	++_size;
	_last_node.isCached = false;

	return newElement._value;
}


//200301 insert new node after current node
//************************************
// Method:insert_new_node_after_current_node
// Qualifier: Insert new node after specific node.
// Input: value, node
// Output: refreshed doubly linked list
// date:200301
//************************************
template <typename T>
inline T& DoublyLinkedList<T>::insert_new_node_after_current_node(DoublyListNode<T>& const newElement, DoublyListNode<T>& const node) {
	if (_size == 0) {

		assert(0);
	}
	else if (_size == 1) {

	}

	newElement._next = node._next;
	newElement._prev = &node;
	if (node._next != nullptr) node._next->_prev = &newElement;
	node._next = &newElement;
	//if (_head == &node) _head = &newElement;
	if (_tail == &node) _tail = &newElement;
	++_size;
	_last_node.isCached = false;

	return newElement._value;
}

//200302 insert new value after current node
//200301 insert new node after current node
//************************************
// Method:insert_new_node_after_current_node
// Qualifier: Insert new node after specific node.
// Input: value, node
// Output: refreshed doubly linked list
// date:200301
// author:
//************************************
template <typename T>
inline T& DoublyLinkedList<T>::insert_new_value_after_current_node(T& const value, DoublyListNode<T>& const node) {
	DoublyListNode<T>* newElement = new DoublyListNode<T>(value);
	return insert_new_node_after_current_node(*newElement, node);
}

template <typename T>
bool DoublyLinkedList<T>::remove(const T &value) {
	DoublyListNode<T> *node = _head;
	bool isDeleted = false;
	while (node != nullptr) {
		if (node->_value != value)
			node = node->_next;
		else {
			if (node->_prev != nullptr) node->_prev->_next = node->_next;
			if (node->_next != nullptr) node->_next->_prev = node->_prev;
			if (_tail == node && node->_prev != nullptr) _tail = node->_prev;
			if (_head == node && node->_next != nullptr) _head = node->_next;
			if (_head == _tail && _tail == node) _head = _tail = nullptr;
			DoublyListNode<T> *tmp = node->_next;
			node->_next = nullptr;
			node->_prev = nullptr;
			node = tmp;
			tmp = nullptr;
			--_size;
			isDeleted = true;
			_last_node.isCached = false;
		}
	}
	return isDeleted;
}

template <typename T>
T DoublyLinkedList<T>::remove(size_t index) {
	if (index >= 0 && index < _size) {
		DoublyListNode<T> *node = _head;
		for (size_t i = 0; i < index; ++i)
			node = node->_next;
		if (node->_prev != nullptr) node->_prev->_next = node->_next;
		if (node->_next != nullptr) node->_next->_prev = node->_prev;
		if (_tail == node && node->_prev != nullptr) _tail = node->_prev;
		if (_head == node && node->_next != nullptr) _head = node->_next;
		if (_head == _tail && _tail == node) _head = _tail = nullptr;
		node->_next = nullptr;
		node->_prev = nullptr;
		--_size;
		_last_node.isCached = false;
		return node->_value;
	}
	else
		throw std::out_of_range("DoublyLinkedList :: remove(index)");
}

template <typename T>
bool DoublyLinkedList<T>::removeAll(const std::initializer_list<T> & elements) {
	// TODO : removeAll()
}

//190910
template <typename T>
inline T& DoublyLinkedList<T>::removeNode(DoublyListNode<T>& const node) {
	if (node._prev != nullptr) node._prev->_next = node._next;
	if (node._next != nullptr) node._next->_prev = node._prev;

	if (_tail == &node && node._prev != nullptr) _tail = node._prev;
	if (_head == &node && node._next != nullptr) _head = node._next;
	if (_head == _tail && _tail == &node) _head = _tail = nullptr;

	node._next = nullptr;
	node._prev = nullptr;
	--_size;
	_last_node.isCached = false;

	return node._value;
}

//200713
template <typename T>
inline T& DoublyLinkedList<T>::pop_back() {
	//removeNode(get_tail_node());
}

template <typename T>
void DoublyLinkedList<T>::clear() {
	if (_size != 0) {
		DoublyListNode<T> *node = _tail;
		while (node != nullptr) {
			if (node->_next != nullptr) {
				delete node->_next;
				node->_next = nullptr;
			}
			if (node->_prev != nullptr)
				node = node->_prev;
			else {
				delete node;
				node = nullptr;
			}
		}
		_head = _tail = nullptr;
		_size = 0;
	}
}

//190828
template <typename T>
DoublyListNode<T>& DoublyLinkedList<T>::getNode(const size_t& const index) const {

	/*=======================================================================*/
	int _segment_index = 0;
	DoublyListNode<T>* current_node = _head;
	// Check if the node trying to get is immediatly AFTER the previous got one
	if (_last_node.isCached) {
		_segment_index = _last_node.last_index;
		current_node = _last_node.last_node;

		while (_segment_index > index && current_node) {
			current_node = current_node->_prev;
			--_segment_index;
		}
	}

	while (_segment_index < index && current_node) {
		current_node = current_node->_next;
		++_segment_index;
	}

	if (_segment_index == index) {
		const_cast<DoublyLinkedList*>(this)->_last_node.isCached = true;
		const_cast<DoublyLinkedList*>(this)->_last_node.last_index = index;
		const_cast<DoublyLinkedList*>(this)->_last_node.last_node = current_node;
	}
	else {
		assert(0);
	}
	/*........................................................................*/

	/*------------------------------------------------------------------------*/
	if (index >= 0 && index < _size) {
		DoublyListNode<T> *node = _head;
		for (size_t i = 0; i < index; ++i)
			node = node->_next;
		return *node;
	}
	else
		throw std::out_of_range("DoublyLinkedList :: getNode (index)");
	/*......................................................................*/
}

template <typename T>
inline T& DoublyLinkedList<T>::get(size_t index) const{
	
	return getNode(index)._value;
}

template <typename T>
T DoublyLinkedList<T>::set(size_t index, const T &value) {
	if (index >= 0 && index < _size) {
		DoublyListNode<T> *node = _head;
		for (size_t i = 0; i < index; ++i)
			node = node->_next;
		T tmp = node->_value;
		node->_value = value;
		return tmp;
	}
	else
		throw std::out_of_range("DoublyLinkedList :: set(index, value)");
}

template <typename T>
size_t DoublyLinkedList<T>::indexOf(const T &element) const {
	// TODO : indexOf()
	return 0;
}

template <typename T>
size_t DoublyLinkedList<T>::lastIndexOf(const T &element) const {
	// TODO : lastIndexOf()
	return 0;
}

template <typename T>
DoublyLinkedList<T>& DoublyLinkedList<T>::subList(size_t fromIndex, size_t toIndex) const {
	// TODO : subList()
	return *(new DoublyLinkedList);
}

template <typename T>
bool DoublyLinkedList<T>::swap(size_t index1, size_t index2) {
	if (index1 >= 0 && index1 < _size && index2 >= 0 && index2 < _size) {
		DoublyListNode<T> tmp = (*this)[index1];
		(*this)[index1] = (*this)[index2];
		(*this)[index2] = tmp;
		_last_node.isCached = false;
		return true;
	}
	else
		return false;
}

//************************************
// Method:swap
// Qualifier: //191029   swap two Linked List
// Input: list list.
// Output: list
// date:191029
// author:
//************************************
template <typename T>
DoublyLinkedList<T>& DoublyLinkedList<T>::swap(DoublyLinkedList<T>& const list) {

	DoublyLinkedList<T> temp_List;
	temp_List._head = _head;
	temp_List._tail = _tail;
	temp_List._size = _size;
	temp_List._last_node = _last_node;

	_head = list._head;
	_tail = list._tail;
	_size = list._size;
	_last_node = list._last_node;

	list._head = temp_List._head;
	list._tail = temp_List._tail;
	list._size = temp_List._size;
	list._last_node = temp_List._last_node;

	temp_List._head = nullptr;
	temp_List._tail = nullptr;

	return *this;
}

//************************************
// Method:copy
// Qualifier: //191029   swap two Linked List
// Input: list 
// Output: copyed list
// date:191030
// author:
//************************************
template <typename T>
DoublyLinkedList<T>& DoublyLinkedList<T>::copy(const DoublyLinkedList<T>& const list) {

	if (_size == 0) {
		for (int node_id = 0; node_id < list.size(); node_id++) {
			this->add(list[node_id]);
		}
	}
	else if (_size == list.size()) {
		DoublyListNode<T>* current_node = _head;
		for (auto&& au : list) {
			current_node->_value = au;
			current_node = current_node->_next;
		}
	}
	else {
		assert(0);
		this->clear();
	}

	return *this;
}

template <typename T>
void DoublyLinkedList<T>::sort() {
	// TODO : sort()
}

//************************************
// Method:front
// Qualifier: first elemet of linked list
// date:190821 12:54
// author:
//************************************
template <typename T>
inline T& DoublyLinkedList<T>::operator[](const int& const index) {
	return getNode(index)._value;
}

//************************************
// Method:operator <<
// Qualifier: print result
// date:190828 12:54
// author:
//************************************
template <typename T>
std::ostream& operator << (std::ostream &out, const DoublyLinkedList<T>& list) {
	DoublyListNode<T>* node = list._head;
	out << '[';
	while (node != nullptr) {
		if (node->_next != nullptr)
			out << node->_value << ", ";
		else
			out << node->_value;
		node = node->_next;
	}
	out << ']';
	node = nullptr;
	return out;
}

//************************************
// Method:front
// Qualifier: first elemet of linked list
// date:190821 12:54
// author:
//************************************
template<typename T>
inline T& DoublyLinkedList<T>::front() noexcept {
	return _head->_value;
}

//************************************
// Method:front
// Qualifier:const first elemet of linked list
// date:191211 12:54
// author:
//************************************
template<typename T>
inline const T& const DoublyLinkedList<T>::front() const noexcept {
	return _head->_value;
}


//************************************
// Method:back
// Qualifier: last element of linked list
// date:190821 12:54
// author:
//************************************
template<typename T>
inline T& DoublyLinkedList<T>::back() const noexcept {
	return _tail->_value;
}

//************************************
// Method:traverse_linkedlist
// Qualifier: traverse linked list
// date:191113 10:43
// author:
//************************************
template<typename T>
void DoublyLinkedList<T>::traverse_linkedlist() {
	DoublyListNode<T>* node = this->get_head_node();
	while (node) {
		node = node->_next;
	}
}


////////////////////////////////////////////////////////////////////////ITERATOR/////////////////////////////////////////////////////////////////////////////////////////////
//************************************
// Method:=
// Qualifier: 
// date:191112 12:54
// author:
//************************************
template<typename T>
inline doubly_linkedlist_iterator<T>& const doubly_linkedlist_iterator<T>::operator = (const DoublyListNode<T>* const node) const {
	this->current_node_pointer = node;
	return *this;
}

//************************************
// Method:++
// Qualifier: Prefix ++ overload 
// date:191113 12:54
// author:
//************************************
template<typename T>
doubly_linkedlist_iterator<T>& doubly_linkedlist_iterator<T>::operator++() {
	//if (this->current_node_pointer)
	this->current_node_pointer = this->current_node_pointer->_next;
	return *this;
}

//************************************
// Method:++
// Qualifier: Prefix ++ overload 
// date:191113 09:21
// author:
//************************************
template<typename T>
doubly_linkedlist_iterator<T> doubly_linkedlist_iterator<T>::operator++(int) {
	doubly_linkedlist_iterator<T> iterator(*this);
	operator++();
	//++*this;
	return iterator;
}

//************************************
// Method:--
// Qualifier: Prefix -- overload 
// date:191113 09:21
// author:
//************************************
template<typename T>
doubly_linkedlist_iterator<T>& doubly_linkedlist_iterator<T>::operator--() {
	//if (this->current_node_pointer)
	this->current_node_pointer = this->current_node_pointer->_prev;
	return *this;
}

//************************************
// Method:--
// Qualifier: Prefix -- overload 
// date:191113 09:21
// author:
//************************************
template<typename T>
doubly_linkedlist_iterator<T> doubly_linkedlist_iterator<T>::operator--(int) {
	doubly_linkedlist_iterator<T> iterator = *this;
	--* this;
	return iterator;
}

//************************************
// Method: !=
// Qualifier: Prefix != overload 
// date:191113 09:21
// author:
//************************************
template<typename T>
inline bool doubly_linkedlist_iterator<T>::operator != (const doubly_linkedlist_iterator<T>& const iterator)  const {
	return this->current_node_pointer != iterator.current_node_pointer;
}

//************************************
// Method: ==
// Qualifier: Prefix != overload 
// date:191114 09:21
// author:
//************************************
template<typename T>
inline bool doubly_linkedlist_iterator<T>::operator == (const doubly_linkedlist_iterator<T>& const iterator)  const {
	return this->current_node_pointer == iterator.current_node_pointer;
}

//************************************
// Method: & reference
// Qualifier: Prefix != overload 
// date:191113 09:21
// author:
//************************************
template<typename T>
inline const T& doubly_linkedlist_iterator<T>::operator * () noexcept {
	return this->current_node_pointer->_value;
}

//************************************
// Method: ->
// Qualifier: pointer* -> 
// date:191113 09:21
// author:
//************************************
template<typename T>
inline const T* doubly_linkedlist_iterator<T>::operator -> () noexcept {
	return &this->current_node_pointer->_value;
}

#endif
