#pragma once
#ifndef RTREE_PARTITION_H
#define RTREE_PARTITION_H

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <assert.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include<vector>
#include <limits>

#include <iostream>

//#include "APCA.h"
/*-----210618-------*/
#include "CAPLA.h"
/*------------------*/

using namespace std;
using std::vector;

#define ASSERT assert // RTree_partition uses ASSERT( condition )
#define INF std::numeric_limits<double>::infinity()

#ifndef Min
#define Min std::min
#endif //Min
#ifndef Max
#define Max std::max
#endif //Max

//int VOLUME_CHOOSE = DBL_MAX;// 181206 for representation option: PAA, APCA,PLA,CHEBY

//unsigned int memory_account[20] = { 0 };

//
// RTree_partition.h
//

//#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
//#define RTREE_PARTITION_QUAL RTree_partition<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, class ELEMTYPEREAL>
#define RTREE_PARTITION_QUAL RTree_partition<DATATYPE, ELEMTYPE, ELEMTYPEREAL>

#define RTREE_DONT_USE_MEMPOOLS // This version does not contain a fixed memory allocator, fill in lines with EXAMPLE to implement one.

//#define RTREE_USE_APCA_SPHERICAL_VOLUME // Better split classification, may be slower on some systems

// Fwd decl
class RTFileStreamPartition;  // File I/O helper class, look below for implementation and notes.
//******class RTree_partition;

/// class RTree_partition
/// Implementation of RTree_partition, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use RTree_partition<Object*, float, 3> myTree;
///
/// DATATYPE Referenced data, should be int, void*, obj* etc. no larger than sizeof<void*> and simple type
/// ELEMTYPE Type of element such as int or float
/// NUMDIMS Number of dimensions such as 2 or 3
/// ELEMTYPEREAL Type of element that allows fractional and large values such as float or double, for use in volume calcs
///
/// NOTES: Inserting and removing data requires the knowledge of its constant Minimal Bounding Rectangle.
///        This version uses new/delete for nodes, I recommend using a fixed size allocator for efficiency.
///        Instead of using a callback function for returned results, I recommend an efficient pre-sized, grow-only memory
///        array similar to MFC CArray or STL Vector for returning search query result.
///
template<class DATATYPE, class ELEMTYPE, class ELEMTYPEREAL = ELEMTYPE>
class RTree_partition : virtual public APLA {
public:
	struct Node;  // Fwd decl.  Used by other internal structs and iterator
	int NUMDIMS = NULL;
	int TMAXNODES = NULL;
	int TMINNODES = NULL;
	int representation_option = NULL;// 180916 for representation option: PAA, APCA,PLA,CHEBY
	RTree_partition(int NUMDIMS, int& TMAXNODESs);
	RTree_partition(int NUMDIMS, int& TMAXNODESs, const int& representation_option);// 180916 for representation option: PAA, APCA,PLA,CHEBY

	// These constant must be declared after Branch and before Node struct
	// Stuck up here for MSVC 6 compiler.  NSVC .NET 2003 is much happier.
	//enum
	//{
	//	TMAXNODES = TMAXNODES,                         ///< Max elements in node
	//	TMINNODES = TMINNODES,                         ///< Min elements in node
	//};

	//210623 partition For PLA Coefficient
	//struct APLA_COEFFICIENT {
	//public:
	//	long double a = INF;//ax+b
	//	long double b = INF;//ax+b
	//	//double right_endpoint = INF;
	//	//int segmentNum = INF; //the dimension of the index(eg. MBR, imput parameter)

	//	long double a_minuend = INF;//(l-1)/2
	//	long double a_divisor = INF;//l(l-1)(l+1)
	//	long double b_minuend = INF;//2l-1
	//	long double b_divisor = INF;//l(l+1)

	//	APLA_COEFFICIENT() {
	//		a = INF;//ax+b
	//		b = INF;//ax+b
	//		//right_endpoint = INF;
	//		//segmentNum = INF; //the dimension of the index(eg. MBR, imput parameter)
	//		//because initial segment length is 2, l=2
	//		a_minuend = 0.5;//(l-1)/2
	//		a_divisor = 6;//l(l-1)(l+1)
	//		b_minuend = 3;//2l-1
	//		b_divisor = 6;//l(l+1)
	//	}
	//};

	//210618 Rtree cover convex hull
	struct AREA_COEFFICIENT_RTREE_NODE {
	public:
		long double right_endpoint = INF; //The right end point of segment;
		long double rectangle_width = INF; //width of rectangle

		APLA::APLA_COEFFICIENT apla;//a & b

		AREA_COEFFICIENT_RTREE_NODE() {
			right_endpoint = INF;
			rectangle_width = INF;
		}
		~AREA_COEFFICIENT_RTREE_NODE() {
			//right_endpoint = INF; //The right end point of segment;
			//apla.~APLA_COEFFICIENT();
		}
	};

	//211008 rectangle in leaf node
	template<typename T>
	struct RECTANGLE_ENTRY_STRUCT;

	//211008 min max point 
	template<typename T>
	struct MINMAX_POINT_STRUCT;

	//211008 min max global 
	template<typename T, typename Y, typename U>
	struct MINMAX_GLOBAL_STRUCT;

	//211008 min max local 
	template<typename T, typename Y, typename U>
	struct MINMAX_LOCAL_STRUCT;

	//211011 min max MBR 
	template<typename T, typename Y>
	struct MINMAX_MBR_STRUCT;

	//211008 vector of vector 
	template<typename T, typename Y>
	struct VECTOR_VECTOR_STRUCT;

	//211008 parameter in internal node
	struct RECTANGLE_INTER_NODE_STRUCT;

	// Minimal bounding rectangle (n-dimensional)
	struct Rect;

	/// May be data or may be another subtree
	/// The parents level determines this.
	/// If the parents level is 0, then this is data
	// 191111 store the MBR of sub nodes
	struct Branch
	{
		Rect m_rect;///< Bounds covering set
		Node* m_child = nullptr;///< onlye has one Child node
		DATATYPE m_data = INF;///< Data Id : Only for the leaf node-> one approximation point id

		/*------------- 210618 --------------------*/
		Node* branch_to_parent_node_pointer = nullptr; // poitner to parent node
		/*-----------------------------------------*/

		shared_ptr<Rect> node_rect_pointer;// 211226 for internal node
		/*~Branch() {
			delete[] m_rect.m_min;
			m_rect.m_min = nullptr;
		}*/
	};

	/// Node for each branch level : 191111 node is root, or leaf node, internal node.
	struct Node
	{
		inline bool IsInternalNode() { return (m_level > 0); } // Not a leaf, but a internal node
		inline bool IsLeaf() { return (m_level == 0); } // A leaf, contains data, sub are entires
		int m_count = NULL;                                  ///< Count  : 191111 the number of sub nodes/ branches. <= MAX NODES Limits: TMAXNODESs
		int m_level = NULL;                                  ///< Leaf is zero, others positive
		Branch* m_branch = nullptr;                    ///< Branch, multi branches. size == m_count

		/*------------- 210618 --------------------*/
		Rect rectangle_in_node;//210618 For all branches.
		Branch* node_to_parent_branch_pointer = nullptr;//210618 node parent branch
		vector<Branch*> branch_ponter_vector;//210813 vector to store address of branch
		/*-----------------------------------------*/
	};

public:

	Node* m_root = nullptr;                                    ///< Root of tree
	Branch branch_root;// 210618 parent branch of root node
	//Node* pointer[10];//210618 pointer is node
	RTree_partition();
	RTree_partition(const RTree_partition& other);
	virtual ~RTree_partition();

	//T type_tree = -1;// 0 R-tree, 1 distance tree, 2 MBR_SAPLA tree
	//Y type_representation = -1; // 0 SAPLA/APLA/PLA/PAA/PAALM/SAX/CHEBY, 3 APCA(average)
	//U type_distance = -1;// 0 SAPLA Eucliden, 1 triangle area
	//T type_volume = -1;//Type of volume. 0:global distance
	//T number_point = INF;// point number to scan in original time series.
	TOOL::OPTION_TREE<int, int, int> option_tree_struct;//211008 option in tree, representation option, tree type option
	void initialRTree(int NUMDIMS, int& TMAXNODESs, const int& representation_option);// 180917 for other class to use
	
	template<typename T, typename Y, typename U>
	void initialRTree(const T& const NUMDIMS, const Y& const TMAXNODESs, const U& const option_tree_struct);// 211008 add option tree.

	/// Insert entry
	/// \param a_min Min of bounding rect
	/// \param a_max Max of bounding rect
	/// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
	void Insert(const ELEMTYPE a_min[], const ELEMTYPE a_max[], const DATATYPE& a_dataId);
	//void Insert(ELEMTYPE a_min[], ELEMTYPE a_max[], const DATATYPE& a_dataId);
	// 191129
	void Insert(const vector<ELEMTYPE>& const a_min, const vector<ELEMTYPE>& const a_max, const DATATYPE& a_dataId);

	//210906 Add approximation point
	//template<typename T, typename Y, typename U, typename T1>
	//void Insert(const T a_min[], const T a_max[], const Y& a_dataId, const DoublyLinkedList<U>& const doubly_linked_list, const vector<T1>& const normalized_time_series_vector);

	//210618 Add approximation point
	template<typename T, typename Y, typename U, typename T1>
	void Insert(const vector<T>& const a_min, const vector<T>& const a_max, const Y& a_dataId, const DoublyLinkedList<U>& const doubly_linked_list, const vector<T1>& const normalized_time_series_vector);

	//211010 Add approximation point, no minmax point
	template<typename T, typename Y, typename U>
	void Insert(const T& a_dataId, const DoublyLinkedList<Y>& const doubly_linked_list, const vector<U>& const normalized_time_series_vector);

	//211010 Add approximation point, no minmax point
	template<typename T, typename Y, typename U, typename T1>
	void Insert(const T& a_dataId, const DoublyLinkedList<Y>& const doubly_linked_list, const vector<U>& const normalized_time_series_vector, T1& const output_argument);

	template<typename T, typename Y>//210815 assign, not copy
	inline void assign_node_rectangle_to_branch_rectangle(T& const node, Y& const rect_in_branch);

	//210618 rectangle copy origianl approxiamtion point.
	template<typename T, typename Y>
	void copy_approximation_point_original(const DoublyLinkedList<T>& const doubly_linked_list, Y& const m_rect);

	//210711 rectangle 1 vector emplace original time series and approximation. 2 shared_ptr point to original time series and approximation.
	template<typename T, typename Y>
	inline void copy_rect_original_series_approximation_vector_pointer(const T& const original_rectangle, Y& const m_rect);

	template<typename T, typename Y>
	void erase_rectangle_from_sub_rectangle(T& const rectangle_sub, Y& const rectangle);

	//210814 update rectangle, delete rectangle_sub_old and insert rectangle_sub_new
	template<typename T, typename Y, typename U>
	void update_rectangle_from_splited_rectangles(T& const rectangle_sub_old, Y& const rectangle_sub_new, U& const rectangle);

	/// Remove entry
	/// \param a_min Min of bounding rect
	/// \param a_max Max of bounding rect
	/// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
	void Remove(const ELEMTYPE a_min[], const ELEMTYPE a_max[], const DATATYPE& a_dataId);
	//void Remove( ELEMTYPE a_min[], ELEMTYPE a_max[], const DATATYPE& a_dataId);

	/// Find all within search rectangle
	/// \param a_min Min of search bounding rect
	/// \param a_max Max of search bounding rect
	/// \param a_searchResult Search result array.  Caller should set grow size. Function will reset, not append to array.
	/// \param a_resultCallback Callback function to return result.  Callback should return 'true' to continue searching
	/// \param a_context User context to pass as parameter to a_resultCallback
	/// \return Returns the number of entries found
	int Search(const ELEMTYPE a_min[], const ELEMTYPE a_max[], std::function<bool(const DATATYPE&)> callback) const;
	//int Search(const ELEMTYPE a_min[], const ELEMTYPE a_max[], std::function<bool(const DATATYPE&)> callback) const;
	/// Remove all entries from tree
	void RemoveAll();
	/// Count the data elements in this container.  This is slow as no internal counter is maintained.
	int Count();

	//211215 Rtree size: internal node size. leaf node number.
	template<typename T, typename Y>
	void count_node_size(Node* a_node, T& const count_node_internal, Y& const count_node_leaf);

	/// Load tree contents from file
	bool Load(const char* a_fileName);
	/// Load tree contents from stream
	bool Load(RTFileStreamPartition& a_stream);
	/// Save tree contents to file
	bool Save(const char* a_fileName);
	/// Save tree contents to stream
	bool Save(RTFileStreamPartition& a_stream);

	void PrintBranchRect(const Branch& m_branch);

	void PrintRect(const Rect& m_rect);

	void printRTree();// 180921

	template<typename T>
	void print_tree(T& const node);//210910

	/// Iterator is not remove safe.
	class Iterator
	{
	private:

		enum { MAX_STACK = 1024 }; //201031: 64  Max stack size. Default 32. Allows almost n^32 where n is number of branches in node

		struct StackElement //: what is m_branchIndex mean???
		{
			Node* m_node;
			int m_branchIndex; //: I think it is the branch id of parent's node.
		};

	public:

		int dimension;

		Iterator() { Init(); }
		~Iterator() { }
		/// Is iterator invalid
		bool IsNull() { return (m_tos <= 0); }
		/// Is iterator pointing to valid data
		bool IsNotNull() { return (m_tos > 0); }

		/// Access the current data element. Caller must be sure iterator is not NULL first.
		DATATYPE& operator*() // * signal
		{
			//cout << "m_tos: "<< m_tos << endl;
			ASSERT(IsNotNull());
			StackElement& curTos = m_stack[m_tos - 1];
			return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
		}

		/// Access the current data element. Caller must be sure iterator is not NULL first.
		const DATATYPE& operator*() const // * signal
		{
			//cout << "m_tos: "<< m_tos << endl;
			ASSERT(IsNotNull());
			StackElement& curTos = m_stack[m_tos - 1];
			return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
		}

		/// Find the next data element   //  ++
		bool operator++() {
			//cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl; 
			return FindNextData();
		}

		/// Delete Memory of Leaf Node
		void deleteLeafNodeMemory()
		{
			ASSERT(IsNotNull());

			StackElement& curTos = m_stack[m_tos - 1];
			Branch& curBranch = curTos.m_node->m_branch[curTos.m_branchIndex];

			/*delete[] curBranch.m_rect.m_min;
			curBranch.m_rect.m_min = nullptr;
			delete[] curBranch.m_rect.m_max;
			curBranch.m_rect.m_max = nullptr;*/

			free_minmax_vector_in_rectangle(curBranch.m_rect);
		}

		/// Get the bounds for this node
		void GetBounds(ELEMTYPE a_min[], ELEMTYPE a_max[], const int& NUMDIMS)
		{
			ASSERT(IsNotNull());

			StackElement& curTos = m_stack[m_tos - 1];
			Branch& curBranch = curTos.m_node->m_branch[curTos.m_branchIndex];

			for (int index = 0; index < NUMDIMS; ++index)
			{
				a_min[index] = curBranch.m_rect.m_min[index];
				a_max[index] = curBranch.m_rect.m_max[index];
			}
		}

	private:

		StackElement m_stack[MAX_STACK];              ///< Stack as we are doing iteration instead of recursion
		int m_tos;                                    ///< Top Of Stack index

		/// Reset iterator
		void Init() { m_tos = 0; }

		//***void Push(Node* a_node, int a_branchIndex);
		//***StackElement& Pop();

		/// Find the next data element in the tree (For internal use only)
		bool FindNextData()
		{
			for (;;)
			{
				if (m_tos <= 0)
				{
					return false;
				}
				StackElement curTos = Pop(); // Copy stack top cause it may change as we use it

				if (curTos.m_node->IsLeaf())
				{
					//cout << "Next Leaf Node:" << endl;
					// Keep walking through data while we can
					if (curTos.m_branchIndex + 1 < curTos.m_node->m_count)
					{
						// There is more data, just point to the next one
						Push(curTos.m_node, curTos.m_branchIndex + 1);
						return true;
					}
					// No more data, so it will fall back to previous level
				}
				else
				{
					//cout << "Next internal Node: " << endl;
					if (curTos.m_branchIndex + 1 < curTos.m_node->m_count)
					{
						// Push sibling on for future tree walk
						// This is the 'fall back' node when we finish with the current level
						Push(curTos.m_node, curTos.m_branchIndex + 1);
					}
					// Since cur node is not a leaf, push first of next level to get deeper into the tree
					Node* nextLevelnode = curTos.m_node->m_branch[curTos.m_branchIndex].m_child;
					Push(nextLevelnode, 0);

					// If we pushed on a new leaf, exit as the data is ready at TOS
					if (nextLevelnode->IsLeaf())
					{
						return true;
					}
				}
			}
		}

		/// Push node and branch onto iteration stack (For internal use only)
		void Push(Node* a_node, int a_branchIndex)
		{
			//cout << "Push: id: "<< a_branchIndex <<endl;
			m_stack[m_tos].m_node = a_node;
			m_stack[m_tos].m_branchIndex = a_branchIndex;
			++m_tos;
			ASSERT(m_tos <= MAX_STACK);
		}

		/// Pop element off iteration stack (For internal use only)
		StackElement& Pop()
		{
			ASSERT(m_tos > 0);
			--m_tos;
			return m_stack[m_tos];
		}

		friend class RTree_partition; // Allow hiding of non-public functions while allowing manipulation by logical owner
	};

	/// Get 'first' for iteration
	void GetFirst(Iterator& a_it)
	{
		a_it.Init();
		Node* first = m_root;
		while (first)
		{
			if (first->IsInternalNode() && first->m_count > 1)
			{
				//cout << "internal Node" << endl;
				a_it.Push(first, 1); // Descend sibling branch later  : what 1 mean?
			}
			else if (first->IsLeaf())
			{
				//cout << "Leaf Node" << endl;
				if (first->m_count)
				{
					a_it.Push(first, 0);//what 0 mean?
				}
				break;
			}
			first = first->m_branch[0].m_child; //: why just index 0???
		}
	}

	/// Get Next for iteration
	void GetNext(Iterator& a_it) { ++a_it; }

	/// Is iterator NULL, or at end?
	bool IsNull(Iterator& a_it) { return a_it.IsNull(); }

	/// Get object at iterator position
	DATATYPE& GetAt(Iterator& a_it) { return *a_it; }


	//************************************
	// Method:deleteAllLeafNodeMemory()
	// Qualifier:Delete Memory of all Leaf Node
	// date: 180921
	// author:
	//************************************
	void deleteAllLeafNodeMemory()
	{
		Iterator it;
		for (GetFirst(it); !IsNull(it); GetNext(it)) {
			it.deleteLeafNodeMemory();
		}
	}

public:

	/// Minimal bounding rectangle (n-dimensional)
	/*struct Rect
	{
		ELEMTYPE m_min[NUMDIMS];                      ///< Min dimensions of bounding box
		ELEMTYPE m_max[NUMDIMS];                      ///< Max dimensions of bounding box
	};*/

	/// May be data or may be another subtree
	/// The parents level determines this.
	/// If the parents level is 0, then this is data
	/*struct Branch
	{
		Rect m_rect;                                  ///< Bounds
		Node* m_child;                                ///< Child node
		DATATYPE m_data;                              ///< Data Id     :Only for the leaf node??
	};*/

	/// Node for each branch level
	//struct Node
	//{
	//	bool IsInternalNode() { return (m_level > 0); } // Not a leaf, but a internal node
	//	bool IsLeaf() { return (m_level == 0); } // A leaf, contains data

	//	int m_count;                                  ///< Count  : number of branches??/Nodes？？？
	//	int m_level;                                  ///< Leaf is zero, others positive
	//	Branch m_branch[MAXNODES];                    ///< Branch
	//};

	/// A link list of nodes for reinsertion after a delete operation
	struct ListNode
	{
		ListNode* m_next;                             ///< Next in list
		Node* m_node;                                 ///< Node
	};

	/// Variables for finding a split partition
	struct PartitionVars
	{
		enum { NOT_TAKEN = -1 }; // indicates that position
		/*~~~~~~~~~~~~~~~210706 vector instead pointer ~~~~~~~~~~~~~~~~~~~~~*/
		//int* m_partition;// = new int;     //  to store the node belong to which group? 0/1
		vector<int> m_partition;// to store the node belong to which group? 0/1. assign the group id of each branch
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		int m_total;                       // the number of the rectangle/branches. the number of total branches in m_ branchBuf = max node number + 1.
		int m_minFill;                     // is minimum nuber of branch in one node
		int m_count[2];                    // 0/1 for 2 groups, count how many branches in one group. the number of notes in every group
		Rect m_cover[2];                   // 0/1 for 2 groups
		ELEMTYPEREAL m_area[2];            //  the are of each two group

		Branch* m_branchBuf;// all branch in original node + new branch. size = m_ branchBuf = max node number + 1.

		int m_branchCount; //the number of total branches in m_branchBuf

		Rect m_coverSplit; // Covered volumes of current node and new branches.  Combined MBR (a_parVars->m_coverSplit) for combined Volume(a_parVars->m_coverSplitArea) from (original node + new branch)
		ELEMTYPEREAL m_coverSplitArea;// Volume of m_coverSplit, from m_ branchBuf (original node + new branch)

		/*~~~~~~~~  210811 For partition  ~~~~~~~~~*/
		int seed_0;
		int seed_1;

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		~PartitionVars() {
			/*~~~~~~~~~~~~~~~210706 vector instead pointer ~~~~~~~~~~~~~~~~~~~~~*/
			//delete[] m_partition;
			//m_partition = nullptr;
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/*======================= min max============================================*/
			/*if (m_cover[0].m_min != nullptr) {
				delete[] m_cover[0].m_min;
				m_cover[0].m_min = nullptr;
				delete[] m_cover[0].m_max;
				m_cover[0].m_max = nullptr;
				delete[] m_cover[1].m_min;
				m_cover[1].m_min = nullptr;
				delete[] m_cover[1].m_max;
				m_cover[1].m_max = nullptr;
			}

			if (m_coverSplit.m_min != nullptr) {
				delete[] m_coverSplit.m_min;
				m_coverSplit.m_min = nullptr;
			}
			if (m_coverSplit.m_max != nullptr) {
				delete[] m_coverSplit.m_max;
				m_coverSplit.m_max = nullptr;
			}*/


			/*===============================================================================*/
			/*for (int i = 0; i < m_total; i++) {
				if (m_ branchBuf[i].m_rect.m_min != nullptr) {
					delete[] m_ branchBuf[i].m_rect.m_min;
					m_ branchBuf[i].m_rect.m_min = nullptr;
				}
				if (m_ branchBuf[i].m_rect.m_max !=nullptr) {
					delete[] m_ branchBuf[i].m_rect.m_max;
					m_ branchBuf[i].m_rect.m_max = nullptr;
				}

			}*/

			delete[] m_branchBuf;
			m_branchBuf = nullptr;
		}
	};

	Node* AllocNode();

	template<typename T>
	inline void allocate_minmax_vector_in_rectangle(T& const a_rect);//210816 allocatge new memory to vector for rectangle

	template<typename T>
	inline void free_minmax_vector_in_rectangle(T& const a_rect);//210816 free new memory to vector for rectangle

	template<typename T>
	inline T& allocate_rectangle_entry(T& const a_rect);//210707 allocatge new memory to pointer or vector for entry

	template<typename T>
	inline T& allocate_rectangle_node(T& const a_rect);//210707 allocatge new memory to pointer or vector for node:

	template<typename T, typename Y>
	inline Y& allocate_rectangle(const T& const option_rectangle, Y& const a_rect);//210707 allocatge new memory to pointer or vector

	template<typename T>
	inline T& clear_vector_in_rectangle(T& const a_rect);//21802 No free memory, only clear vector.

	template<typename T>
	inline T& free_rectangle(T& const a_rect);//210730 free memory to pointer or vectotr

	template<typename T>//210804  For rectangle in struct PartitionVars, reset poitner for convex hull
	inline T& clear_rectangle_PartitionVars(T& const a_rect);

	//************************************
	// Method:AllocBranch
	// Qualifier: Because rect was changed from array to pointer. pointer needs new memory. So whenever apply new Branch, needs this frunction for new memory
	// date: personal function
	// author:
	//************************************
	template<typename T>
	inline void AllocBranch(T& const a_branch);

	template<typename T>
	inline void free_shared_ptr_linked_list(shared_ptr<T>& const smart_pointer);

	template<typename T>
	inline void reset_shared_ptr_linked_list(shared_ptr<T>& const smart_pointer);

	//210730
	template<typename T>
	inline void free_pointer_vector(vector<T>& const ponter_vector);

	void FreeNode(Node* a_node);
	void InitNode(Node* a_node);
	void InitRect(Rect* a_rect);//make every value is 0
	bool InsertRectRec(Branch& a_branch, Node* a_node, Node** a_newNode, int a_level);
	bool InsertRect(Branch& a_branch, Node** a_root, int a_level);
	//Rect& NodeCover(Node* a_node);
	Rect& NodeCover(Node* a_node, Rect& a_rect);
	bool AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode);
	void DisconnectBranch(Node* a_node, int a_index);
	int PickBranch(Rect* a_rect, Node* a_node);

	template<typename T, typename Y>//210810 pick branch by distance
	int pick_branch(T& const a_rect, Y& const a_node);

	Rect& CombineRect(Rect* a_rectA, Rect* a_rectB);
	Rect CombineRect(Rect* a_rectA, Rect* a_rectB, Rect* a_rect);

	//The purpuse of function is get volume of combined rectangle for pick branch
	template<typename T, typename Y, typename U>//210728 return difference of volume
	long double combine_compute_partition_rectangle_by_entry(T& const rectangle_entry, Y& const rectangle_node, U& const rectangle_combine);

	//The purpuse of function is get volume of combined rectangle for pick branch
	template<typename T, typename Y>//210905 return difference of volume. rectangle_entry is entry. No need temp combine rectangle, speed up
	long double combine_compute_partition_rectangle_by_entry_210905(T& const rectangle_entry, Y& const rectangle_node);

	//211225 The purpuse of function is get volume of combined rectangle for pick branch
	template<typename T, typename Y>// return difference of volume. rectangle_entry is entry. No need temp combine rectangle, speed up
	long double combine_compute_partition_rectangle_by_convex(T& const rectangle_entry, Y& const rectangle_node);

	//211227 The purpuse of function is get volume of internal nodt for pick branch
	template<typename T, typename Y>// return difference of volume. rectangle_entry is entry. No need temp combine rectangle, speed up
	long double combine_compute_partition_node_by_convex(T& const rectangle_node_0, Y& const rectangle_node_1);

	////211227 The purpuse of function is get volume of internal nodt for pick branch
	//template<typename T, typename Y,>// return difference of volume. rectangle_entry is entry. No need temp combine rectangle, speed up
	//long double get_combine_partition_node_by_convex(T& const rectangle_node_0, Y& const rectangle_node_1);

	//Rect& CombineRect( Rect* a_rectA,  Rect* a_rectB, Rect& a_rect);
	void assertRect(Rect& a_rect);
	void SplitNode(Node* a_node, const Branch* a_branch, Node** a_newNode);
	/*=======================================================       Volume Calculation              =======================================================================*/
	ELEMTYPEREAL RectSphericalVolume(Rect* a_rect);
	ELEMTYPEREAL APCARectSphericalVolume(Rect* a_rect);
	ELEMTYPEREAL RectVolume(Rect* a_rect);
	ELEMTYPEREAL APCARectVolume(Rect* a_rect);
	//191128 accumulation, min - min - 2
	ELEMTYPEREAL apca_rect_accumulation_volume(Rect* a_rect);
	ELEMTYPEREAL PLARectVolume(Rect* a_rect);
	ELEMTYPEREAL ChebyshevRectVolume(Rect* a_rect);
	//191117  apla segment volume calculation
	ELEMTYPEREAL apla_rect_volume(Rect* a_rect);

	//***************************************************************
	// Method: partition_rect_volume
	// Qualifier:   
	// Input: 1 rectangle_new: sub rectangle. 2 node: leaf node.
	// Output: compute volume of APLA apprximation.
	// date:210719
	// author:
	//***************************************************************
	template<typename T, typename Y>
	long double partition_rect_volume(const T& const rectangle_new, Y& const node);

	template<typename T, typename Y>//210806 Scan each sub node to get variant.
	Y& scan_each_node_for_variant(T& const node, Y& const rectangle);

	template<typename T, typename Y, typename U>// entry in branch.
	long double& get_convex_hull_volume_by_minmax_reconstruct_point(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined);

	template<typename T>//210806 vector in node
	long double& get_convex_hull_by_minmax_reconstruct_point(T& const rectangle_in_node);

	template<typename T, typename Y, typename U>// entry in branch. In Leaf Node. 
	long double& get_convex_hull_volume_by_distance_approximation_global(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined);

	//211227
	template<typename T, typename Y, typename U>// entry in branch. In intranl Node. 
	long double& get_convex_hull_volume_by_distance_internal_node_global(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined);

	template<typename T>//210806 vector in node
	long double& get_convex_hull_by_distance_approximation_global(T& const rectangle_in_node);

	

	template<typename T, typename Y, typename U> // entry in branch
	long double& get_convex_hull_volume_by_distance_approximation_local(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined);

	template<typename T>//210806 vector in node
	long double& get_convex_hull_by_distance_approximation_local(T& const rectangle_in_node);

	template<typename T, typename Y>//210811 vector in node
	long double compute_distance_max_between_rectangle_node(T& const rectangle_in_node_0, Y& const rectangle_in_node_1);

	template<typename T, typename Y, typename U>//211011 entry in branch. In Leaf Node. 
	long double& get_convex_hull_volume_by_MBR_local(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined);

	template<typename T, typename Y, typename U>
	long double get_convex_hull_by_entry(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined);

	template<typename T>
	long double get_convex_hull_by_node(T& const node);

	template<typename T>//210804
	long double get_convex_hull_by_node_rectangle(T& const node);

	template<typename T>//210804
	long double get_convex_hull_by_rectangle(T& const rectangle);

	ELEMTYPEREAL CalcRectVolume(Rect* a_rect);
	/*=================================================================================================================================================================*/

	void GetBranches(Node* a_node, const Branch* a_branch, PartitionVars* a_parVars);
	void ChoosePartition(PartitionVars* a_parVars, int a_minFill);

	//210729 choose_partition_by_ distance, already get two seeds.
	template<typename T, typename Y>
	void choose_partition_by_distance(const T& const a_minFill, Y& const a_parVars);

	void LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars);
	inline void InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill);
	//[Calculate inefficiency of grouping entries together.] For each pair of entries and E2, compose a rectan? gle J including EVI and E2.l. Calcu? late d= area(</) - area (E^I) - area {E2.I).
	//[Choose the most wasteful pair.] Choose the pair with the largest
	void PickSeeds(PartitionVars* a_parVars);
	void Classify(int a_index, int a_group, PartitionVars* a_parVars);

	template<typename T>
	long double pick_seeds(T& const a_parVars);

	bool RemoveRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root);
	bool RemoveRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, ListNode** a_listNode);
	ListNode* AllocListNode();
	void FreeListNode(ListNode* a_listNode);
	bool Overlap(Rect* a_rectA, Rect* a_rectB) const;
	void ReInsert(Node* a_node, ListNode** a_listNode);
	bool Search(Node* a_node, Rect* a_rect, int& a_foundCount, std::function<bool(const DATATYPE&)> callback) const;
	void RemoveAllRec(Node* a_node);
	void reset_root(Node* a_node);
	void Reset();
	void CountRec(Node* a_node, int& a_count);
	bool SaveRec(Node* a_node, RTFileStreamPartition& a_stream);
	bool LoadRec(Node* a_node, RTFileStreamPartition& a_stream);
	void CopyRec(Node* current, Node* other);

	template<typename T>//210813
	void update_branch_pointer_vector(T& const node);

	/*===================================  210812 Is function===============================================*/
	template<typename T>
	inline bool is_rectangle_entry(const T& const rectangle);

	template<typename T>
	inline bool is_rectangle_node(const T& const rectangle);

	template<typename T>
	inline bool is_branch_entry(const T& const branch);

	template<typename T>
	inline bool is_branch_node(const T& const branch);

	template<typename T, typename Y>
	bool is_branch_pointer_all_entry(const T& const size_array, const Y& const branch_pointer);

	template<typename T, typename Y>
	bool is_branch_pointer_all_node(const T& const size_array, const Y& const branch_pointer);
	/*==============================================================================================*/

	/*====================================== 210618 Assert ================================================*/
	//211008 assert option tree
	template<typename T>
	inline bool assert_option_tree(const T& const option_tree_struct);

	//210725 two vectors have same value and size
	template<typename T, typename Y>
	bool assert_same_vector(const T& const vector_0, const Y& const vector_1);

	// 210729 two SAPLAs have same a&b, right endpoints and rectangle width.
	template<typename T, typename Y>
	bool assert_same_SAPLA(const DoublyLinkedList<T>& const doubly_linked_list_0, const DoublyLinkedList<Y>& const doubly_linked_list_1);

	//211227
	template<typename T>
	inline bool assert_rectangle_global_minmax(const T& const rectangle_in_node);

	//210728
	template<typename T>
	bool assert_rectangle(const T& const rectangle);

	//210804 Two rectangle has same values. branch.m_rect == node.rectangle_in_node
	template<typename T, typename Y>
	bool assert_rectangle_same(const T& const rectangle_0, const Y& const rectangle_1);

	//210728
	template<typename T>
	bool assert_branch(T& const branch);

	//210729
	template<typename T, typename Y>
	bool assert_branch_array(const T& const size_array, Y* branch_pointer);

	//210801
	template<typename T, typename Y>
	bool assert_node_branch_pointer(const T& const node, Y* branch_pointer);

	//210730
	template<typename T, typename Y>
	bool assert_branch_rect_in_node(const T& const rectangle, Y* branch_pointer);

	//210725
	template<typename T>
	bool assert_node_leaf(T& const a_node);

	//210803 from each internal node to each leaf node
	template<typename T>
	bool assert_node(T& const a_node);
	/*=====================================================================================================*/

	//Node* m_root;                                    ///< Root of tree
	ELEMTYPEREAL m_unitSphereVolume;                 ///< Unit sphere constant for required number of dimensions
};

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      RTFileStreamPartition Classs    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

// Because there is not stream support, this is a quick and dirty file I/O helper.
// Users will likely replace its usage with a Stream implementation from their favorite API.
class RTFileStreamPartition
{
	FILE* m_file;

public:


	RTFileStreamPartition()
	{
		m_file = NULL;
	}

	~RTFileStreamPartition()
	{
		Close();
	}

	bool OpenRead(const char* a_fileName)
	{
		//m_file = fopen(a_fileName, "rb");
		fopen_s(&m_file, a_fileName, "rb");
		if (!m_file)
		{
			return false;
		}
		return true;
	}

	bool OpenWrite(const char* a_fileName)
	{
		//m_file = fopen(a_fileName, "wb");
		fopen_s(&m_file, a_fileName, "wb");
		if (!m_file)
		{
			return false;
		}
		return true;
	}

	void Close()
	{
		if (m_file)
		{
			fclose(m_file);
			m_file = NULL;
		}
	}

	template< typename TYPE >
	size_t Write(const TYPE& a_value)
	{
		ASSERT(m_file);
		return fwrite((void*)&a_value, sizeof(a_value), 1, m_file);
	}

	template< typename TYPE >
	size_t WriteArray(const TYPE* a_array, int a_count)
	{
		ASSERT(m_file);
		return fwrite((void*)a_array, sizeof(TYPE) * a_count, 1, m_file);
	}

	template< typename TYPE >
	size_t Read(TYPE& a_value)
	{
		ASSERT(m_file);
		return fread((void*)&a_value, sizeof(a_value), 1, m_file);
	}

	template< typename TYPE >
	size_t ReadArray(TYPE* a_array, int a_count)
	{
		ASSERT(m_file);
		return fread((void*)a_array, sizeof(TYPE) * a_count, 1, m_file);
	}
};


/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&           Implement Tree Class         &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/


RTREE_TEMPLATE
template<typename T>
struct RTREE_PARTITION_QUAL::RECTANGLE_ENTRY_STRUCT {
public:
	shared_ptr<vector<long double>> original_time_series_vector_pointer;
	//shared_ptr<vector<long double>> reconstruct_time_series_vector_pointer;
	shared_ptr<DoublyLinkedList<T>> list_original_pointer;//
	shared_ptr<DoublyLinkedList<T>> list_partition_pointer;
};

//211008 min max point
RTREE_TEMPLATE
template<typename T>
struct RTREE_PARTITION_QUAL::MINMAX_POINT_STRUCT {
public:
	T volume_minmax_point = -std::numeric_limits<T>::infinity();// max reconstruct point difference
	T volume_minmax_point_sqrt = -std::numeric_limits<T>::infinity(); // max reconstruct point difference square root
	vector<T> value_point_min_vector;
	vector<T> value_point_max_vector;
};

//211008 min max global
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
struct RTREE_PARTITION_QUAL::MINMAX_GLOBAL_STRUCT {
public:
	shared_ptr<DoublyLinkedList<T>> list_partition_global_min_pointer;//AREA_COEFFICIENT_RTREE_NODE
	shared_ptr<DoublyLinkedList<T>> list_partition_global_max_pointer;

	/*211225*/
	shared_ptr<vector<Y>> original_time_series_vector_min_pointer;
	shared_ptr<DoublyLinkedList<T>> list_original_min_pointer;//AREA_COEFFICIENT_RTREE_NODE
	shared_ptr<vector<Y>> original_time_series_vector_max_pointer;
	shared_ptr<DoublyLinkedList<T>> list_original_max_pointer;//AREA_COEFFICIENT_RTREE_NODE
	U seed_branch = -1;
	/**/

	Y volume_minmax_whole = -std::numeric_limits<Y>::infinity();// max whole difference
	Y volume_minmax_whole_sqrt = -std::numeric_limits<Y>::infinity();// max whole difference square root
	U seed_0 = (std::numeric_limits<U>::min)();
	U seed_1 = (std::numeric_limits<U>::min)();
};

//211008 min max local
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
struct RTREE_PARTITION_QUAL::MINMAX_LOCAL_STRUCT {
public:
	shared_ptr<DoublyLinkedList<T>> list_partition_local_min_pointer;//AREA_COEFFICIENT_RTREE_NODE
	shared_ptr<DoublyLinkedList<T>> list_partition_local_max_pointer;
	Y volume_minmax_segment = -std::numeric_limits<Y>::infinity();// max semgent difference
	Y volume_minmax_segment_sqrt = -std::numeric_limits<Y>::infinity();// max segment difference square root
	vector<Y> difference_each_segment_max_vector;
	vector<Y> difference_power_each_segment_max_vector;
	vector<U> id_branch_segment_min_vector;
	vector<U> id_branch_segment_max_vector;

};

//211011 min max MBR
RTREE_TEMPLATE
template<typename T, typename Y>
struct RTREE_PARTITION_QUAL::MINMAX_MBR_STRUCT {
public:
	vector<DoublyLinkedList<T>> MBR_list_max_vector;//211011
	vector<DoublyLinkedList<T>> MBR_list_min_vector;
	Y volume_MBR_segment = -std::numeric_limits<Y>::infinity();// max semgent difference
	Y volume_MBR_segment_sqrt = -std::numeric_limits<Y>::infinity();// max segment difference square root
};

//211008 vector of vector
RTREE_TEMPLATE
template<typename T, typename Y>
struct RTREE_PARTITION_QUAL::VECTOR_VECTOR_STRUCT {
public:

	/*------------------------   210803   -----------------------------*/
	inline auto size_entry() { return (list_original_pointer_vector.size()); }
	/*-----------------------------------------------------------------*/
	/*------------------------------------------------    210710 shared pointer    ----------------------------------------------------------------*/
	//AREA_COEFFICIENT_RTREE_NODE
	vector<shared_ptr<DoublyLinkedList<T>>> list_original_pointer_vector;//for each partitioned approximation point
	vector<shared_ptr<DoublyLinkedList<T>>> list_partition_pointer_vector;//for each partitioned approximation point
	vector<shared_ptr<vector<Y>>> original_time_series_pointer_vector;//210706
	//vector<shared_ptr<vector<long double>>> reconstruct_time_series_pointer_vector;//210806
	/*----------------------------------------------------------------------------------------------------------------------------------------------*/
};

//211008 parameter in internal node
RTREE_TEMPLATE
struct RTREE_PARTITION_QUAL::RECTANGLE_INTER_NODE_STRUCT {
public:
	/*#############    210718 Internal Node, Rectangle in node: Convex hull & Volume    #################*/

	/*--------- max(at+b)-min(at+b) ---------------------*/
	MINMAX_POINT_STRUCT<long double> minmax_point_struct;
	/*---------------------------------------------------*/

	/*--------------------------- for one whole approximation -----------------------------*/
	MINMAX_GLOBAL_STRUCT<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX, long double, int> minmax_global_struct;
	/*-------------------------------------------------------------------------------------*/

	/*---------------------211011 for MBR struct ---------------------------------*/
	MINMAX_MBR_STRUCT<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX, long double> minmax_MBR_struct;
	/*----------------------------------------------------------------------------*/

	/*--------- for each segment approximation that has max differnce --------------------*/
	MINMAX_LOCAL_STRUCT<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX, long double, int> minmax_local_struct;
	/*------------------------------------------------------------------------------------*/

	/*###################################################################################################*/

	/*###################         vector of vector / pointer       #####################*/
	VECTOR_VECTOR_STRUCT<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX, long double> vector_vector_struct;
	/*##################################################################################*/

};


// Minimal bounding rectangle (n-dimensional)
RTREE_TEMPLATE
struct RTREE_PARTITION_QUAL::Rect {
public:

	/*######### origianl variantes ##################*/
	//ELEMTYPE m _min[NUMDIMS];   //< Min dimensions of bounding box 
	//ELEMTYPE* m _min = nullptr;
	//ELEMTYPE m_max[NUMDIMS];   //< Max dimensions of bounding box 
	//ELEMTYPE* m_max = nullptr;

	/*-------210816---------*/
	vector<ELEMTYPE> m_min;//instead pointer
	vector<ELEMTYPE> m_max;
	/*----------------------*/
	/*##################################*/

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&           210618 For SAPLA approximation       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	/*###############################      210716 Leaf Node   ###############################*/
	//vector<long double> id_data_vector;//210709
	//vector<long double> m_min_vector;//210709
	//vector<long double> m_max_vector;//210709
	//vector<long double> original_time_series_vector_pointer;
	//vector<long double> reconstruct_time_series_vector_pointer;//210718
	/*=====================  Instead struct  ====================================*/
	//shared_ptr<vector<long double>> original_time_series_vector_pointer;
	////shared_ptr<vector<long double>> reconstruct_time_series_vector_pointer;
	//shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> list_original_pointer;
	//shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> list_partition_pointer;
	RECTANGLE_ENTRY_STRUCT<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX> rectangle_entry_struct;
	/*===========================================================================*/
	/*#########################################################################################*/

	/*###############################    210718 Internal Node, Rectangle in node: Convex hull & Volume      ###############################*/
	int seed_0 = (std::numeric_limits<int>::min)();
	int seed_1 = (std::numeric_limits<int>::min)();
	RECTANGLE_INTER_NODE_STRUCT rectangle_inter_node_struct;
	/*--------------------------- max(at+b)-min(at+b) ---------------------------------*/
	//long double volume_minmax_point = -std::numeric_limits<long double>::infinity();// max reconstruct point difference
	//long double volume_minmax_point_sqrt = -std::numeric_limits<long double>::infinity(); // max reconstruct point difference square root
	//vector<long double> value_point_min_vector;
	//vector<long double> value_point_max_vector;
	//MINMAX_POINT_STRUCT<long double> minmax_point_struct;
	/*----------------------------------------------------------------------------------*/

	/*--------------------------- for one whole approximation ---------------------------------*/
	//shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> list_partition_global_min_pointer;
	//shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> list_partition_global_max_pointer;
	//long double volume_minmax_whole = -std::numeric_limits<long double>::infinity();// max whole difference
	//long double volume_minmax_whole_sqrt = -std::numeric_limits<long double>::infinity();// max whole difference square root
	//MINMAX_GLOBAL_STRUCT<AREA_COEFFICIENT_RTREE_NODE, long double, int> minmax_global_struct;
	/*-----------------------------------------------------------------------------------------*/

	/*-------------- for each segment approximation that has max differnce ---------------------*/
	//shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> list_partition_local_min_pointer;
	//shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> list_partition_local_max_pointer;
	//long double volume_minmax_segment = -std::numeric_limits<long double>::infinity();// max semgent difference
	//long double volume_minmax_segment_sqrt = -std::numeric_limits<long double>::infinity();// max segment difference square root
	//vector<long double> difference_each_segment_max_vector;
	//vector<long double> difference_power_each_segment_max_vector;
	//vector<int> id_branch_segment_min_vector;
	//vector<int> id_branch_segment_max_vector;

	//MINMAX_LOCAL_STRUCT<AREA_COEFFICIENT_RTREE_NODE, long double, int> minmax_local_struct;
	/*------------------------------------------------------------------------------------------*/

	/*#################################################################################################################################################*/

	/*######################################################         vector of vector / pointer       #################################################*/
	//vector<vector<long double>> original_time_series_vector_vector;//210706
	//vector<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> linked_list_original_vector;//for each partitioned approximation point
	//vector<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>> linked_list_partition_vector;//for each partitioned approximation point

	/*------------------------------------------------    210710 shared pointer    ----------------------------------------------------------------*/
	///*------------------------   210803   -----------------------------*/
	//inline auto size_entry() { return (list_original_pointer_vector.size()); }
	///*-----------------------------------------------------------------*/
	//vector<shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>>> list_original_pointer_vector;//for each partitioned approximation point
	//vector<shared_ptr<DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>>> list_partition_pointer_vector;//for each partitioned approximation point
	//vector<shared_ptr<vector<long double>>> original_time_series_pointer_vector;//210706
	////vector<shared_ptr<vector<long double>>> reconstruct_time_series_pointer_vector;//210806
	//VECTOR_VECTOR_STRUCT<AREA_COEFFICIENT_RTREE_NODE, long double> vector_vector_struct;
	/*----------------------------------------------------------------------------------------------------------------------------------------------*/
	/*##################################################################################################################################################*/

	
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	/*~Rect() {
		delete[] m_min;
		m_min = nullptr;
		delete[] m_max;
		m_max = nullptr;
	}*/
};

RTREE_TEMPLATE
RTREE_PARTITION_QUAL::RTree_partition()
{
	//assert(0);
	//this->NUMDIMS = 20;
	//ASSERT(TMAXNODES > TMINNODES);
	//ASSERT(TMINNODES > 0);

	// Precomputed volumes of the unit spheres for the first few dimensions
	const float UNIT_SPHERE_VOLUMES[] = {
	  0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
	  4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
	  5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
	  3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
	  1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
	  0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
	  0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20
	  0.014060f, // Dimension  21 
	};

	//211010 initial all Options in Parition tree
	TOOL::initial_option_tree(1, 0, 0, 0, 0, this->option_tree_struct);

	m_root = AllocNode();
	//pointer[0] = m_root;
	m_root->m_level = 0;
	m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[NUMDIMS];
}

RTREE_TEMPLATE
RTREE_PARTITION_QUAL::RTree_partition(const RTree_partition& other) : RTree_partition()
{
	CopyRec(m_root, other.m_root);
}


RTREE_TEMPLATE
RTREE_PARTITION_QUAL::~RTree_partition()
{
#ifdef _DEBUG
	//cout << "RTree_partition clear memory \n";
#endif
	RemoveAll();// 180917 auto free memeory
	Reset(); // Free, or reset node memory
}


RTREE_TEMPLATE
RTREE_PARTITION_QUAL::RTree_partition(int NUMDIMS, int& TMAXNODESs) {//Defualt TMAXNODES=8 TMINNODES = TMAXNODES / 2
	this->NUMDIMS = NUMDIMS;
	TMAXNODES = TMAXNODESs;
	TMINNODES = 2 ;// TMAXNODESs / 2;
	ASSERT(TMAXNODES > TMINNODES);
	ASSERT(TMINNODES > 0);
	//cout << "RTREE: NUMDIMS = " << NUMDIMS << ", TMAXNODES = " << TMAXNODES << ", TMINNODES =" << TMINNODES << endl;
	// Precomputed volumes of the unit spheres for the first few dimensions
	m_root = AllocNode();
	//pointer[1] = m_root;
	m_root->m_level = 0;

	// Precomputed volumes of the unit spheres for the first few dimensions
	const double UNIT_SPHERE_VOLUMES[] = {
		0.000000, 2.000000, 3.141593, // Dimension  0,1,2
		4.188790, 4.934802, 5.263789, // Dimension  3,4,5
		5.167713, 4.724766, 4.058712, // Dimension  6,7,8
		3.298509, 2.550164, 1.884104, // Dimension  9,10,11
		1.335263, 0.910629, 0.599265, // Dimension  12,13,14
		0.381443, 0.235331, 0.140981, // Dimension  15,16,17
		0.082146, 0.046622, 0.025807, // Dimension  18,19,20
		0.014060, 0.007426, 0.003838,  // Dimension  21,22,23
		0.001943, 0.000964, 0.000469,  // Dimension  24,25,26
		0.000224, 0.000105, 0.000469,  // Dimension  27,28,29
		0.0000220374261236442, 0.00000983989335030794, 0.0000043255369431949,  // Dimension  30,31,32
		0.00000187290118762488, 0.000000799112584179925, 0.000000336125448655506,  // Dimension  33,34,35
		0.000000139433492032519, 0.0000000570647325568745, 0.0000000230492917233849,  // Dimension  36,37,38
		0.00000000919142299825764, 0.000000003619781, 0.000000001408278,  // Dimension  39,40,41
		0.000000000541411, 0.000000000205739, 0.000000000077299,  // Dimension  42,43,44
		0.000000000028722, 0.000000000010557, 0.000000000003839,  // Dimension  45,46,47
		0.000000000001382, 0.000000000000492, 0.000000000000174 , // Dimension  48,49,50
		0.000000000000061, 0.000000000000021, 0.000000000000007 , // Dimension  51,52,53
		0.000000000000002, 0.000000000000001 //Dimension  54,55,
	};
	m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[NUMDIMS];
};


//************************************
// Method:RTree_partition
// Qualifier:
// date:180916
// author: 
//************************************
RTREE_TEMPLATE
RTREE_PARTITION_QUAL::RTree_partition(int NUMDIMS, int& TMAXNODESs, const int& representation_option) {//Defualt TMAXNODES=8 TMINNODES = TMAXNODES / 2
	this->NUMDIMS = NUMDIMS;
	TMAXNODES = TMAXNODESs;
	TMINNODES = 2; //TMAXNODESs / 2;
	this->representation_option = representation_option;
	ASSERT(TMAXNODES > TMINNODES);
	ASSERT(TMINNODES > 0);
	//cout << "RTREE: NUMDIMS = " << NUMDIMS << ", TMAXNODES = " << TMAXNODES << ", TMINNODES =" << TMINNODES << endl;
	// Precomputed volumes of the unit spheres for the first few dimensions
	m_root = AllocNode();
	//pointer[1] = m_root;
	m_root->m_level = 0;

	// Precomputed volumes of the unit spheres for the first few dimensions
	const double UNIT_SPHERE_VOLUMES[] = {
		0.000000, 2.000000, 3.141593, // Dimension  0,1,2
		4.188790, 4.934802, 5.263789, // Dimension  3,4,5
		5.167713, 4.724766, 4.058712, // Dimension  6,7,8
		3.298509, 2.550164, 1.884104, // Dimension  9,10,11
		1.335263, 0.910629, 0.599265, // Dimension  12,13,14
		0.381443, 0.235331, 0.140981, // Dimension  15,16,17
		0.082146, 0.046622, 0.025807, // Dimension  18,19,20
		0.014060, 0.007426, 0.003838,  // Dimension  21,22,23
		0.001943, 0.000964, 0.000469,  // Dimension  24,25,26
		0.000224, 0.000105, 0.000469,  // Dimension  27,28,29
		0.0000220374261236442, 0.00000983989335030794, 0.0000043255369431949,  // Dimension  30,31,32
		0.00000187290118762488, 0.000000799112584179925, 0.000000336125448655506,  // Dimension  33,34,35
		0.000000139433492032519, 0.0000000570647325568745, 0.0000000230492917233849,  // Dimension  36,37,38
		0.00000000919142299825764, 0.000000003619781, 0.000000001408278,  // Dimension  39,40,41
		0.000000000541411, 0.000000000205739, 0.000000000077299,  // Dimension  42,43,44
		0.000000000028722, 0.000000000010557, 0.000000000003839,  // Dimension  45,46,47
		0.000000000001382, 0.000000000000492, 0.000000000000174 , // Dimension  48,49,50
		0.000000000000061, 0.000000000000021, 0.000000000000007 , // Dimension  51,52,53
		0.000000000000002, 0.000000000000001 //Dimension  54,55,
	};
	m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[NUMDIMS];
};

//210618 For multi_dimension class to use instead construct class
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::initialRTree(int NUMDIMS, int& TMAXNODESs, const int& representation_option) {// 180917 for other class to use
	this->NUMDIMS = NUMDIMS;
	TMAXNODES = TMAXNODESs;
	TMINNODES = 2;// TMAXNODESs / 2;
	this->representation_option = representation_option;
	ASSERT(TMAXNODES > TMINNODES);
	ASSERT(TMINNODES > 0);

	//cout << "RTREE: NUMDIMS = " << NUMDIMS << ", TMAXNODES = " << TMAXNODES << ", TMINNODES =" << TMINNODES << endl;
	// Precomputed volumes of the unit spheres for the first few dimensions
	m_root = AllocNode();
	//pointer[1] = m_root;
	m_root->m_level = 0;

	/*---- 210618 partition -------*/
	//m_root->node_to_parent_branch_pointer = &branch_root;// get parent branch of root node.
	/*-----------------------------*/

	// Precomputed volumes of the unit spheres for the first few dimensions
	const double UNIT_SPHERE_VOLUMES[] = {
		0.000000, 2.000000, 3.141593, // Dimension  0,1,2
		4.188790, 4.934802, 5.263789, // Dimension  3,4,5
		5.167713, 4.724766, 4.058712, // Dimension  6,7,8
		3.298509, 2.550164, 1.884104, // Dimension  9,10,11
		1.335263, 0.910629, 0.599265, // Dimension  12,13,14
		0.381443, 0.235331, 0.140981, // Dimension  15,16,17
		0.082146, 0.046622, 0.025807, // Dimension  18,19,20
		0.014060, 0.007426, 0.003838,  // Dimension  21,22,23
		0.001943, 0.000964, 0.000469,  // Dimension  24,25,26
		0.000224, 0.000105, 0.000469,  // Dimension  27,28,29
		0.0000220374261236442, 0.00000983989335030794, 0.0000043255369431949,  // Dimension  30,31,32
		0.00000187290118762488, 0.000000799112584179925, 0.000000336125448655506,  // Dimension  33,34,35
		0.000000139433492032519, 0.0000000570647325568745, 0.0000000230492917233849,  // Dimension  36,37,38
		0.00000000919142299825764, 0.000000003619781, 0.000000001408278,  // Dimension  39,40,41
		0.000000000541411, 0.000000000205739, 0.000000000077299,  // Dimension  42,43,44
		0.000000000028722, 0.000000000010557, 0.000000000003839,  // Dimension  45,46,47
		0.000000000001382, 0.000000000000492, 0.000000000000174 , // Dimension  48,49,50
		0.000000000000061, 0.000000000000021, 0.000000000000007 , // Dimension  51,52,53
		0.000000000000002, 0.000000000000001 //Dimension  54,55,
	};
	m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[NUMDIMS];
}

// 211008 add option tree.
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
void RTREE_PARTITION_QUAL::initialRTree(const T& const NUMDIMS, const Y& const TMAXNODESs, const U& const option_tree_struct) {
	this->NUMDIMS = NUMDIMS;
	TMAXNODES = TMAXNODESs;
	TMINNODES = 2;// TMAXNODESs / 2;
	this->option_tree_struct = option_tree_struct;
	this->representation_option = this->option_tree_struct.type_representation;

#ifdef _DEBUG
	ASSERT(TMAXNODES > TMINNODES);
	ASSERT(TMINNODES > 0);
	assert_option_tree(option_tree_struct);
#endif
	
	//cout << "RTREE: NUMDIMS = " << NUMDIMS << ", TMAXNODES = " << TMAXNODES << ", TMINNODES =" << TMINNODES << endl;
	// Precomputed volumes of the unit spheres for the first few dimensions
	m_root = AllocNode();
	//pointer[1] = m_root;
	m_root->m_level = 0;

	/*---- 210618 partition -------*/
	//m_root->node_to_parent_branch_pointer = &branch_root;// get parent branch of root node.
	/*-----------------------------*/

	// Precomputed volumes of the unit spheres for the first few dimensions
	const double UNIT_SPHERE_VOLUMES[] = {
		0.000000, 2.000000, 3.141593, // Dimension  0,1,2
		4.188790, 4.934802, 5.263789, // Dimension  3,4,5
		5.167713, 4.724766, 4.058712, // Dimension  6,7,8
		3.298509, 2.550164, 1.884104, // Dimension  9,10,11
		1.335263, 0.910629, 0.599265, // Dimension  12,13,14
		0.381443, 0.235331, 0.140981, // Dimension  15,16,17
		0.082146, 0.046622, 0.025807, // Dimension  18,19,20
		0.014060, 0.007426, 0.003838,  // Dimension  21,22,23
		0.001943, 0.000964, 0.000469,  // Dimension  24,25,26
		0.000224, 0.000105, 0.000469,  // Dimension  27,28,29
		0.0000220374261236442, 0.00000983989335030794, 0.0000043255369431949,  // Dimension  30,31,32
		0.00000187290118762488, 0.000000799112584179925, 0.000000336125448655506,  // Dimension  33,34,35
		0.000000139433492032519, 0.0000000570647325568745, 0.0000000230492917233849,  // Dimension  36,37,38
		0.00000000919142299825764, 0.000000003619781, 0.000000001408278,  // Dimension  39,40,41
		0.000000000541411, 0.000000000205739, 0.000000000077299,  // Dimension  42,43,44
		0.000000000028722, 0.000000000010557, 0.000000000003839,  // Dimension  45,46,47
		0.000000000001382, 0.000000000000492, 0.000000000000174 , // Dimension  48,49,50
		0.000000000000061, 0.000000000000021, 0.000000000000007 , // Dimension  51,52,53
		0.000000000000002, 0.000000000000001 //Dimension  54,55,
	};
	m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[NUMDIMS];
}

RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::Insert(const ELEMTYPE a_min[], const ELEMTYPE a_max[], const DATATYPE& a_dataId)
{
#ifdef _DEBUG
	//cout <<"RTree_partition::Insert() = "<< NUMDIMS << endl;
	for (int index = 0; index < NUMDIMS; ++index)
	{
		if (index & 1) {//odd point value
			ASSERT(a_min[index] <= a_max[index]);
		}
		else {//even point id
			ASSERT(a_min[index] <= a_max[index]);
		}

		a_min[index] != INF;
		a_max[index] != INF;

		switch (representation_option) {
		case 1://MSPLA
		case 2://PLA
		case 3:
		case 4:
		case 5:
		case 6:
			if (index & 1 == 0) {//even
				ASSERT(a_min[index] == a_max[index]);
			}
			break;
		default:
			assert(0);
			break;
		}
	}
#endif

	Branch branch;
	//branch.m_rect.m_min = new ELEMTYPE[NUMDIMS];
	//branch.m_rect.m_max = new ELEMTYPE[NUMDIMS];
	allocate_minmax_vector_in_rectangle(branch.m_rect);

	branch.m_data = a_dataId;
	branch.m_child = NULL;

	for (int axis = 0; axis < NUMDIMS; ++axis)
	{
		//cout <<"(a_min : "<< a_min[axis] << ", a_max : " << a_max[axis] <<")  ";
		branch.m_rect.m_min[axis] = a_min[axis];
		branch.m_rect.m_max[axis] = a_max[axis];
	}
	//cout << endl;
	InsertRect(branch, &m_root, 0);
}

//***************************************************************
// Method:Insert
// Qualifier: insertion of vector
// Input:
// Output: 
// Note: 
// date:191129
// author:
//***************************************************************
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::Insert(const vector<ELEMTYPE>& const a_min, const vector<ELEMTYPE>& const a_max, const DATATYPE& const a_dataId) {

	/*..............................................................*/
#ifdef _DEBUG
	//cout <<"RTree_partition::Insert() = "<< NUMDIMS << endl;
	for (int index = 0; index < NUMDIMS; ++index)
	{
		if (index & 1) {//odd point value
			ASSERT(a_min[index] <= a_max[index]);
		}
		else {//even point id
			ASSERT(a_min[index] <= a_max[index]);
		}

		a_min[index] != INF;
		a_max[index] != INF;

		switch (representation_option) {
		case 1://PAA
		case 2://PLA
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:// Initial 200706
		case 9:// SAX
			if (index & 1 == 0) {//even
				ASSERT(a_min[index] == a_max[index]);
			}
			break;
		default:
			assert(0);
			break;
		}
	}
#endif
	/*..............................................................*/

	Branch branch;
	//branch.m_rect.m_min = new ELEMTYPE[NUMDIMS];
	//branch.m_rect.m_max = new ELEMTYPE[NUMDIMS];
	allocate_minmax_vector_in_rectangle(branch.m_rect);

	branch.m_data = a_dataId;
	branch.m_child = NULL;

	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, INF);

	for (int axis = 0; axis < NUMDIMS; ++axis)
	{
		//cout <<"(a_min : "<< a_min[axis] << ", a_max : " << a_max[axis] <<")  ";
		assert(a_min[axis] != INF && a_max[axis] != INF);
		branch.m_rect.m_min[axis] = a_min[axis];
		branch.m_rect.m_max[axis] = a_max[axis];
	}
	//cout << endl;
	InsertRect(branch, &m_root, 0);
}

////210906 Add approximation point. Use pointer instead vector
//RTREE_TEMPLATE
//template<typename T, typename Y, typename U, typename T1>
//void RTREE_PARTITION_QUAL::Insert(const T a_min[], const T a_max[], const Y& a_dataId, const DoublyLinkedList<U>& const doubly_linked_list, const vector<T1>& const normalized_time_series_vector) {
//
//}

//***************************************************************
// Method:Insert
// Qualifier: insertion of vector. 
// Input:
// Output: 
// Note: Distance Tree. Add approximation point. Add template
// date:210618
// author:
//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y, typename U, typename T1>
void RTREE_PARTITION_QUAL::Insert(const vector<T>& const a_min, const vector<T>& const a_max, const Y& a_dataId, const DoublyLinkedList<U>& const doubly_linked_list, const vector<T1>& const normalized_time_series_vector) {
	/*............................................................................................*/
#ifdef _DEBUG
	assert_option_tree(this->option_tree_struct);
	//cout <<"RTree_partition::Insert() = "<< NUMDIMS << endl;
	for (int index = 0; index < NUMDIMS; ++index)
	{
		if (index & 1) {//odd point value
			ASSERT(a_min[index] <= a_max[index]);
		}
		else {//even point id
			ASSERT(a_min[index] <= a_max[index]);
		}

		a_min[index] != INF;
		a_max[index] != INF;

		switch (representation_option) {
		case 1://PAA
		case 2://PLA
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:// Initial 200706
		case 9:// SAX
			if (index & 1 == 0) {//even
				ASSERT(a_min[index] == a_max[index]);
			}
			break;
		default:
			assert(0);
			break;
		}
	}

	switch (representation_option) {
	case 3://APCA
		break;
	case 5://CHEBY
		assert(NUMDIMS == doubly_linked_list.size());
		break;
	case 2: // PLA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 6: // ICDE07
	case 8:
		assert(NUMDIMS / 2 == doubly_linked_list.size());
		break;
	default:
		break;
	}

	//update_branch_pointer_vector(*m_root);
#endif
	/*......................................................................................*/

	Branch branch;

	/*---------210705 -----*/
	allocate_rectangle(1, branch.m_rect);
	//AllocBranch(branch);
	/*---------------------*/
	/*---------210705 -----*/
	/*branch.m_rect.m_ min = new ELEMTYPE[NUMDIMS];
	branch.m_rect.m_max = new ELEMTYPE[NUMDIMS];*/
	/*---------------------*/

	branch.m_data = a_dataId;
	branch.m_child = NULL;

	for (int axis = 0; axis < NUMDIMS; ++axis)
	{
		//cout <<"(a_min : "<< a_min[axis] << ", a_max : " << a_max[axis] <<")  ";
		assert(a_min[axis] != INF && a_max[axis] != INF);
		branch.m_rect.m_min[axis] = a_min[axis];
		branch.m_rect.m_max[axis] = a_max[axis];
	}
	//cout << endl;

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     210618 partition     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	/*------------- 210705 AllocBranch instead ------------------*/
	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, INF);
	//branch.m_rect.original_time_series_vector_pointer.assign(1024, INF);
	//branch.m_rect.linked_list_original_vector.resize(TMAXNODES, DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>());
	//branch.m_rect.rectangle_entry_struct.list_original_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	//branch.m_rect.rectangle_entry_struct.list_partition_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	/*----------------------------------------------------------*/

	//210618 copy new original approximation point to 2
	copy_approximation_point_original(doubly_linked_list, branch.m_rect);
	//210701 copy originl time series
	branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->assign(normalized_time_series_vector.begin(), normalized_time_series_vector.end());
	//branch.m_rect.reconstruct_time_series_vector_pointer->assign(1024, INF);

	//copy_n(normalized_time_series_vector.begin(), normalized_time_series_vector.size(), branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->begin());
	/*+++++++++++++++++++++++++++++++++++++++   210709 Rect in Node  +++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*========================  210710 vector of vector  =====================================*/
	/*m_root->rectangle_in_node.original_time_series_vector_vector.emplace_back(branch.m_rect.original_time_series_vector_pointer);
	m_root->rectangle_in_node.linked_list_original_vector.emplace_back(*branch.m_rect.rectangle_entry_struct.list_original_pointer);
	m_root->rectangle_in_node.linked_list_partition_vector.emplace_back(*branch.m_rect.rectangle_entry_struct.list_partition_pointer);*/
	/*=========================================================================================*/
	/*========================  210711 vector of pointer ======================================*/
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.emplace_back(&branch.m_rect.original_time_series_vector_pointer);//210706
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.emplace_back(branch.m_rect.rectangle_entry_struct.list_original_pointer);//for each partitioned approximation point
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.emplace_back(branch.m_rect.rectangle_entry_struct.list_partition_pointer);//for each partitioned approximation point
	/*=========================================================================================*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	InsertRect(branch, &m_root, 0);
}


//211010 For SAPLA MBR Add approximation point, no minmax point
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
void RTREE_PARTITION_QUAL::Insert(const T& a_dataId, const DoublyLinkedList<Y>& const doubly_linked_list, const vector<U>& const normalized_time_series_vector) {

	/*............................................................................................*/
#ifdef _DEBUG
	assert_option_tree(this->option_tree_struct);
	switch (representation_option) {
	case 3://APCA
		break;
	case 5://CHEBY
		assert(NUMDIMS == doubly_linked_list.size());
		break;
	case 2: // PLA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 6: // ICDE07
	case 8:
		assert(NUMDIMS / 2 == doubly_linked_list.size());
		break;
	default:
		break;
	}

	//update_branch_pointer_vector(*m_root);
#endif
	/*......................................................................................*/

	Branch branch;

	/*------210705 branch rectangle -------*/
	allocate_rectangle(1, branch.m_rect);// 1 entry rectangle
	//AllocBranch(branch);
	/*-------------------------------------*/

	branch.m_data = a_dataId;
	branch.m_child = NULL;

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     210618 partition     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	/*------------- 210705 AllocBranch instead ------------------*/
	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, INF);
	//branch.m_rect.original_time_series_vector_pointer.assign(1024, INF);
	//branch.m_rect.linked_list_original_vector.resize(TMAXNODES, DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>());
	//branch.m_rect.rectangle_entry_struct.list_original_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	//branch.m_rect.rectangle_entry_struct.list_partition_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	/*----------------------------------------------------------*/


	/*.......*/
	//210618 copy new original approximation point to 2
	//copy_approximation_point_original(doubly_linked_list, branch.m_rect);

	branch.m_rect.rectangle_entry_struct.list_original_pointer ->copy(doubly_linked_list);
	branch.m_rect.rectangle_entry_struct.list_partition_pointer->copy(doubly_linked_list);
	/*.......*/

	//210701 copy originl time series
	//branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->assign(normalized_time_series_vector.size(), INF);
	//branch.m_rect.reconstruct_time_series_vector_pointer->assign(1024, INF);
	//copy_n(normalized_time_series_vector.begin(), normalized_time_series_vector.size(), branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->begin());
	branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->assign(normalized_time_series_vector.begin(), normalized_time_series_vector.end());
	/*+++++++++++++++++++++++++++++++++++++++   210709 Rect in Node  +++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*========================  210710 vector of vector  =====================================*/
	/*m_root->rectangle_in_node.original_time_series_vector_vector.emplace_back(branch.m_rect.original_time_series_vector_pointer);
	m_root->rectangle_in_node.linked_list_original_vector.emplace_back(*branch.m_rect.rectangle_entry_struct.list_original_pointer);
	m_root->rectangle_in_node.linked_list_partition_vector.emplace_back(*branch.m_rect.rectangle_entry_struct.list_partition_pointer);*/
	/*=========================================================================================*/
	/*========================  210711 vector of pointer ======================================*/
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.emplace_back(&branch.m_rect.original_time_series_vector_pointer);//210706
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.emplace_back(branch.m_rect.rectangle_entry_struct.list_original_pointer);//for each partitioned approximation point
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.emplace_back(branch.m_rect.rectangle_entry_struct.list_partition_pointer);//for each partitioned approximation point
	/*=========================================================================================*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	
	InsertRect(branch, &m_root, 0);
}

//211010 Add approximation point, no minmax point
RTREE_TEMPLATE
template<typename T, typename Y, typename U, typename T1>
void RTREE_PARTITION_QUAL::Insert(const T& a_dataId, const DoublyLinkedList<Y>& const doubly_linked_list, const vector<U>& const normalized_time_series_vector, T1& const output_struct) {
	/*............................................................................................*/
#ifdef _DEBUG
	assert_option_tree(this->option_tree_struct);
	switch (representation_option) {
	case 3://APCA
		break;
	case 5://CHEBY
		assert(NUMDIMS == doubly_linked_list.size());
		break;
	case 2: // PLA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 6: // ICDE07
	case 8:
		assert(NUMDIMS / 2 == doubly_linked_list.size());
		break;
	default:
		break;
	}

	//update_branch_pointer_vector(*m_root);
#endif
	/*......................................................................................*/

	Branch branch;

	/*------210705 branch rectangle -------*/
	allocate_rectangle(1, branch.m_rect);// 1 entry rectangle
	//AllocBranch(branch);
	/*-------------------------------------*/

	branch.m_data = a_dataId;
	branch.m_child = NULL;

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     210618 partition     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	/*------------- 210705 AllocBranch instead ------------------*/
	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	//branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, INF);
	//branch.m_rect.original_time_series_vector_pointer.assign(1024, INF);
	//branch.m_rect.linked_list_original_vector.resize(TMAXNODES, DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>());
	//branch.m_rect.rectangle_entry_struct.list_original_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	//branch.m_rect.rectangle_entry_struct.list_partition_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	/*----------------------------------------------------------*/


	/*.......*/
	//210618 copy new original approximation point to 2
	//copy_approximation_point_original(doubly_linked_list, branch.m_rect);

	branch.m_rect.rectangle_entry_struct.list_original_pointer->copy(doubly_linked_list);
	branch.m_rect.rectangle_entry_struct.list_partition_pointer->copy(doubly_linked_list);
	/*.......*/

	//210701 copy originl time series
	//branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->assign(normalized_time_series_vector.size(), INF);
	//branch.m_rect.reconstruct_time_series_vector_pointer->assign(1024, INF);
	//copy_n(normalized_time_series_vector.begin(), normalized_time_series_vector.size(), branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->begin());
	branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer->assign(normalized_time_series_vector.begin(), normalized_time_series_vector.end());
	/*+++++++++++++++++++++++++++++++++++++++   210709 Rect in Node  +++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*========================  210710 vector of vector  =====================================*/
	/*m_root->rectangle_in_node.original_time_series_vector_vector.emplace_back(branch.m_rect.original_time_series_vector_pointer);
	m_root->rectangle_in_node.linked_list_original_vector.emplace_back(*branch.m_rect.rectangle_entry_struct.list_original_pointer);
	m_root->rectangle_in_node.linked_list_partition_vector.emplace_back(*branch.m_rect.rectangle_entry_struct.list_partition_pointer);*/
	/*=========================================================================================*/
	/*========================  210711 vector of pointer ======================================*/
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.emplace_back(&branch.m_rect.original_time_series_vector_pointer);//210706
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.emplace_back(branch.m_rect.rectangle_entry_struct.list_original_pointer);//for each partitioned approximation point
	//m_root->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.emplace_back(branch.m_rect.rectangle_entry_struct.list_partition_pointer);//for each partitioned approximation point
	/*=========================================================================================*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	TOOL::recordStartTime(TOOL::time_record[10]);
	InsertRect(branch, &m_root, 0);
	output_struct = TOOL::recordFinishTime(TOOL::time_record[10]);
}

RTREE_TEMPLATE
template<typename T, typename Y>//210815 assign, not copy
inline void RTREE_PARTITION_QUAL::assign_node_rectangle_to_branch_rectangle(T& const node, Y& const rect_in_branch) {

	/*==============        210802     =========================*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		/*..............................................................................*/
#ifdef _DEBUG
		assert_node(node);
#endif
		/*..............................................................................*/
		rect_in_branch = node.rectangle_in_node;
		break;
	}
	default:
		break;
	}
	/*==========================================================*/

}

//210618 branch copy origianl approxiamtion point.
RTREE_TEMPLATE
template<typename T, typename Y>
void RTREE_PARTITION_QUAL::copy_approximation_point_original(const DoublyLinkedList<T>& const doubly_linked_list, Y& const m_rect) {

	//	/*..............................................................*/
	//#ifdef _DEBUG
	//	
	//#endif
	//	/*..............................................................*/
	assert(0);
	//m_rect.rectangle_entry_struct.list_original_pointer->clear();
	m_rect.rectangle_entry_struct.list_partition_pointer->clear();

	/*---*/
	//m_rect.rectangle_entry_struct.list_original_pointer.copy(doubly_linked_list);
	//m_rect.rectangle_entry_struct.list_partition_pointer.copy(doubly_linked_list);
	/*---*/

	/*---*/
	APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX temp_segment_coefficients;
	for (int id_segment = 0; id_segment < doubly_linked_list.size(); id_segment++) {
		temp_segment_coefficients.right_endpoint = doubly_linked_list[id_segment].right_endpoint;
		temp_segment_coefficients.rectangle_width = doubly_linked_list[id_segment].rectangle_width;
		temp_segment_coefficients.apla = doubly_linked_list[id_segment].apla;
		m_rect.rectangle_entry_struct.list_original_pointer->emplace_back(temp_segment_coefficients);
		m_rect.rectangle_entry_struct.list_partition_pointer->emplace_back(temp_segment_coefficients);
	}
	/*--*/

	//m_rect.linked_list_original_vector.emplace_back(m_rect.rectangle_entry_struct.list_original_pointer);// 210623 delete crash
}

//210711 rectangle 1, vector emplace original time series and approximation. 2, shared_ptr point to original time series and approximation.
RTREE_TEMPLATE
template<typename T, typename Y>
inline void RTREE_PARTITION_QUAL::copy_rect_original_series_approximation_vector_pointer(const T& const rectangle_original, Y& const m_rect) {
	/*............................................................................................*/
#ifdef _DEBUG
	//assert_rectangle(rectangle_original);
	//assert_rectangle(m_rect);
#endif
	/*............................................................................................*/

	/*======================================================  210710 vector of vector  ============================================================*/
	// rectangle_original in Leaf Node (is an entry)
	if (rectangle_original.rectangle_entry_struct.original_time_series_vector_pointer && (*rectangle_original.rectangle_entry_struct.original_time_series_vector_pointer)[0] != INF) {
		/*............................................................................................*/
#ifdef _DEBUG
		assert(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.empty());
		assert(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.empty() && rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.empty());
#endif
		/*............................................................................................*/

		//m_rect.original_time_series_vector_vector.emplace_back(rectangle_original.original_time_series_vector_pointer);
		m_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.emplace_back(rectangle_original.rectangle_entry_struct.original_time_series_vector_pointer);//210706
		//m_rect.reconstruct_time_series_pointer_vector.emplace_back(rectangle_original.reconstruct_time_series_vector_pointer);//210806
		/*========================            210711 vector of pointer        ======================================*/
		m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.emplace_back(rectangle_original.rectangle_entry_struct.list_original_pointer);//for each partitioned approximation point
		m_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.emplace_back(rectangle_original.rectangle_entry_struct.list_partition_pointer);//for each partitioned approximation point
		/*==========================================================================================================*/
	}
	else {// rectangle_original in Internal Node, copy all its entries to parent node

		/*............................................................................................*/
#ifdef _DEBUG
		assert(!rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.empty());
		assert(!rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.empty() && !rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.empty());
		assert(!rectangle_original.rectangle_entry_struct.original_time_series_vector_pointer && !m_rect.rectangle_entry_struct.original_time_series_vector_pointer);
		assert_rectangle(rectangle_original);
		assert_rectangle(m_rect);
#endif
		/*............................................................................................*/
		switch (option_tree_struct.type_tree) {
		case 1: {
			for (int id_pointer = 0; id_pointer < rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size(); id_pointer++) {
				//Only copy non duplicate entry from rectangle_original
				if (find(m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.begin(), m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.end(), rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_pointer]) == m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.end()) {
					m_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.emplace_back(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_pointer]);
					//m_rect.reconstruct_time_series_pointer_vector.emplace_back(rectangle_original.reconstruct_time_series_pointer_vector[id_pointer]);
					m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.emplace_back(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_pointer]);
					m_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.emplace_back(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_pointer]);
				}
			}
			break;
		}
		case 2:
			//assert(0);
			/*--------- option_tree_struct.type_tree is 1 ----------*/
			for (int id_pointer = 0; id_pointer < rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size(); id_pointer++) {
				//Only copy non duplicate entry from rectangle_original
				if (find(m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.begin(), m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.end(), rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_pointer]) == m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.end()) {
					m_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.emplace_back(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_pointer]);
					//m_rect.reconstruct_time_series_pointer_vector.emplace_back(rectangle_original.reconstruct_time_series_pointer_vector[id_pointer]);
					m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.emplace_back(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_pointer]);
					m_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.emplace_back(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_pointer]);
				}
			}
			/*-------------------*/
			break;
		default:
			assert(0);
		}
		
		/*========================  210711 vector of pointer ======================================*/
		//copy_ n(rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.begin(), rectangle_original.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size(), m_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.end());
		//for each partitioned approximation point
		//copy_n (rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.begin(), rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(), m_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.end());
		//for each partitioned approximation point
		//copy_n (rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.begin(), rectangle_original.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size(), m_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.end());
		/*=========================================================================================*/
		/*............................................................................................*/
#ifdef _DEBUG
		assert_rectangle(m_rect);
#endif
		/*............................................................................................*/
	}
	//m_rect.linked_list_original_vector.emplace_back(*rectangle_original.rectangle_entry_struct.list_original_pointer);
	//m_rect.linked_list_partition_vector.emplace_back(*rectangle_original.rectangle_entry_struct.list_partition_pointer);
	/*=========================================================================================*/

	/*.....................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_original);
	//assert_rectangle(m_rect);
#endif
	/*......................................*/
}

//210814 erase sub rectangle vector from one rectangle
RTREE_TEMPLATE
template<typename T, typename Y>
void RTREE_PARTITION_QUAL::erase_rectangle_from_sub_rectangle(T& const rectangle_sub, Y& const rectangle) {
	/*......................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_sub);
	assert_rectangle(rectangle);
	//assert(rectangle_sub.rectangle_inter_node_struct.vector_vector_struct.size_entry()<= rectangle.rectangle_inter_node_struct.vector_vector_struct.size_entry());
#endif
	/*......................................*/

	for (int id_entry = 0; id_entry < rectangle_sub.rectangle_inter_node_struct.vector_vector_struct.size_entry(); id_entry++) {
		rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.erase(std::remove(rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.begin(), rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.end(), rectangle_sub.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_entry]), rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.end());
		//rectangle.reconstruct_time_series_pointer_vector.erase(std::remove(rectangle.reconstruct_time_series_pointer_vector.begin(), rectangle.reconstruct_time_series_pointer_vector.end(), rectangle_sub.reconstruct_time_series_pointer_vector[id_entry]), rectangle.reconstruct_time_series_pointer_vector.end());
		rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.erase(std::remove(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.begin(), rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.end(), rectangle_sub.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_entry]), rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.end());
		rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.erase(std::remove(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.begin(), rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.end(), rectangle_sub.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_entry]), rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.end());
	}

	get_convex_hull_by_rectangle(rectangle);

	/*......................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle);
#endif
	/*......................................*/

}

//210814 update rectangle, delete rectangle_sub_old and insert rectangle_sub_new
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
void RTREE_PARTITION_QUAL::update_rectangle_from_splited_rectangles(T& const rectangle_sub_old, Y& const rectangle_sub_new, U& const rectangle) {

	/*~~~~~~~~~~~~~~~~~  partition 210618 ~~~~~~~~~~~~~~~~~~~~~~~~~*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		/*......................................*/
#ifdef _DEBUG
		assert_rectangle(rectangle_sub_old);
		assert_rectangle(rectangle_sub_new);
		assert_rectangle(rectangle);
#endif
		/*......................................*/

		erase_rectangle_from_sub_rectangle(rectangle_sub_old, rectangle);
		erase_rectangle_from_sub_rectangle(rectangle_sub_new, rectangle);
		copy_rect_original_series_approximation_vector_pointer(rectangle_sub_old, rectangle);

		/*......................................*/
#ifdef _DEBUG
		assert_rectangle(rectangle_sub_old);
		assert_rectangle(rectangle_sub_new);
		assert_rectangle(rectangle);
#endif
		/*......................................*/
		break;
	}
	default: {
		break;
	}
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
}

RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::Remove(const ELEMTYPE a_min[], const  ELEMTYPE a_max[], const DATATYPE& a_dataId)
{
#ifdef _DEBUG
	for (int index = 0; index < NUMDIMS; ++index)
	{
		if (index & 1) {//odd point value
			ASSERT(a_min[index] <= a_max[index]);
		}
		else {//even point id
			ASSERT(a_min[index] < a_max[index]);
		}
	}
#endif 

	Rect rect;

	for (int axis = 0; axis < NUMDIMS; ++axis)
	{
		rect.m_min[axis] = a_min[axis];
		rect.m_max[axis] = a_max[axis];
	}

	RemoveRect(&rect, a_dataId, &m_root);
}


RTREE_TEMPLATE
int RTREE_PARTITION_QUAL::Search(const ELEMTYPE a_min[], const ELEMTYPE a_max[], std::function<bool(const DATATYPE&)> callback) const
{
#ifdef _DEBUG
	for (int index = 0; index < NUMDIMS; ++index)
	{
		ASSERT(a_min[index] <= a_max[index]);
	}
#endif //_DEBUG

	Rect rect;

	for (int axis = 0; axis < NUMDIMS; ++axis)
	{
		rect.m_min[axis] = a_min[axis];
		rect.m_max[axis] = a_max[axis];
	}

	// NOTE: May want to return search result another way, perhaps returning the number of found elements here.

	int foundCount = 0;
	Search(m_root, &rect, foundCount, callback);

	return foundCount;
}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::PrintBranchRect(const Branch& m_branch) {
	for (int i = 0; i < NUMDIMS; i++) {
		cout << "(" << m_branch.m_rect.m_min[i] << "," << m_branch.m_rect.m_max[i] << ") ";
	}
}

RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::PrintRect(const Rect& m_rect) {
	for (int i = 0; i < NUMDIMS; i++) {
		cout << "(" << m_rect.m_min[i] << "," << m_rect.m_max[i] << ") ";
	}
}

//************************************
// Method:printRTree()
// Qualifier:
// date: 180921
// author:
//************************************
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::printRTree() {// 180921
	cout << "Root Node : sub node number = " << m_root->m_count << " Root level = : " << m_root->m_level << "\n\nBegin to build a RTree_partition:\n";
	cout << "\n RTree_partition conclusion\n The number of RTree_partition Data Point = : " << Count() << endl;
}

RTREE_TEMPLATE
template<typename T>
void RTREE_PARTITION_QUAL::print_tree(T& const node) {
	if (node.IsInternalNode()) {
		cout << "Internal Node. Level: " << node.m_level << endl;
		for (int id_branch = 0; id_branch < node.m_count; ++id_branch) {
			cout << " Level: " << node.m_level << " sub branch " << id_branch << endl;
			print_tree(*node.m_branch[id_branch].m_child);
		}
	}
	else {
		cout << "Leaf Node. Level: " << node.m_level << endl;
		for (int id_branch = 0; id_branch < node.m_count; id_branch++) {
			cout <<"  branch id " << id_branch << " id:" << node.m_branch[id_branch].m_data;
		}
		cout << endl;
	}

}

RTREE_TEMPLATE
int RTREE_PARTITION_QUAL::Count()
{
	int count = 0;
	CountRec(m_root, count);
	return count;
}

RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::CountRec(Node* a_node, int& a_count)           //  count the data of the subtree, not Node???
{

	if (a_node->IsInternalNode())  // not a leaf node
	{
		cout << "    Internal Node Number: " << a_node->m_count << ",  Node level =" << a_node->m_level << endl;

		for (int index = 0; index < a_node->m_count; ++index)
		{
			cout << "       Internal Node Branch rectangle : ";
			PrintBranchRect(a_node->m_branch[index]);
			cout << " Internal Node Rectangle Volume: " << CalcRectVolume(&a_node->m_branch[index].m_rect) << endl;

			CountRec(a_node->m_branch[index].m_child, a_count);
		}
	}
	else // A leaf node
	{

		a_count += a_node->m_count;

		cout << "        Leaf Data Number = " << a_node->m_count << ",  Node level =" << a_node->m_level << endl;

		for (int index = 0; index < a_node->m_count; ++index) {
			cout << "            id: " << a_node->m_branch[index].m_data << " Leaf Node branch rectangle : ";
			PrintBranchRect(a_node->m_branch[index]);
			cout << ". Leaf Rectangle Volume: " << CalcRectVolume(&a_node->m_branch[index].m_rect) << endl;
			cout << a_count << endl;
		}
	}
}

//211215 Rtree size: internal node size. leaf node number.
RTREE_TEMPLATE
template<typename T, typename Y>
void RTREE_PARTITION_QUAL::count_node_size(Node* a_node, T& const count_node_internal, Y& const count_node_leaf) {

	if (a_node->IsInternalNode()) {  // not a leaf node
		assert(a_node->m_level > 0);

		if (a_node->m_level > 1) {
			count_node_internal += a_node->m_count;
		}

		for (int index = 0; index < a_node->m_count; ++index) {
			count_node_size(a_node->m_branch[index].m_child, count_node_internal, count_node_leaf);
		}
	}
	else { // A leaf node
		count_node_leaf++;
	}

}

//RTREE_TEMPLATE
//void RTREE_PARTITION_QUAL::CountRec(Node* a_node, int& a_count)           //  count the data of the subtree, not Node???
//{
//	if (a_node->IsInternalNode())  // not a leaf node
//	{
//		for (int index = 0; index < a_node->m_count; ++index)
//		{
//			CountRec(a_node->m_branch[index].m_child, a_count);
//		}
//	}
//	else // A leaf node
//	{
//		a_count += a_node->m_count;
//	}
//}



RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::Load(const char* a_fileName)
{
	RemoveAll(); // Clear existing tree

	RTFileStreamPartition stream;
	if (!stream.OpenRead(a_fileName))
	{
		return false;
	}

	bool result = Load(stream);

	stream.Close();

	return result;
}



RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::Load(RTFileStreamPartition& a_stream)
{
	// Write some kind of header
	int _dataFileId = ('R' << 0) | ('T' << 8) | ('R' << 16) | ('E' << 24);
	int _dataSize = sizeof(DATATYPE);
	int _dataNumDims = NUMDIMS;
	int _dataElemSize = sizeof(ELEMTYPE);
	int _dataElemRealSize = sizeof(ELEMTYPEREAL);
	int _dataMaxNodes = TMAXNODES;
	int _dataMinNodes = TMINNODES;

	int dataFileId = 0;
	int dataSize = 0;
	int dataNumDims = 0;
	int dataElemSize = 0;
	int dataElemRealSize = 0;
	int dataMaxNodes = 0;
	int dataMinNodes = 0;

	a_stream.Read(dataFileId);
	a_stream.Read(dataSize);
	a_stream.Read(dataNumDims);
	a_stream.Read(dataElemSize);
	a_stream.Read(dataElemRealSize);
	a_stream.Read(dataMaxNodes);
	a_stream.Read(dataMinNodes);

	bool result = false;

	// Test if header was valid and compatible
	if ((dataFileId == _dataFileId)
		&& (dataSize == _dataSize)
		&& (dataNumDims == _dataNumDims)
		&& (dataElemSize == _dataElemSize)
		&& (dataElemRealSize == _dataElemRealSize)
		&& (dataMaxNodes == _dataMaxNodes)
		&& (dataMinNodes == _dataMinNodes)
		)
	{
		// Recursively load tree
		result = LoadRec(m_root, a_stream);
	}

	return result;
}


RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::LoadRec(Node* a_node, RTFileStreamPartition& a_stream)
{
	a_stream.Read(a_node->m_level);
	a_stream.Read(a_node->m_count);

	if (a_node->IsInternalNode())  // not a leaf node
	{
		for (int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

			curBranch->m_child = AllocNode();
			LoadRec(curBranch->m_child, a_stream);
		}
	}
	else // A leaf node
	{
		for (int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.ReadArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.ReadArray(curBranch->m_rect.m_max, NUMDIMS);

			a_stream.Read(curBranch->m_data);
		}
	}

	return true; // Should do more error checking on I/O operations
}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::CopyRec(Node* current, Node* other)
{
	current->m_level = other->m_level;
	current->m_count = other->m_count;

	if (current->IsInternalNode())  // not a leaf node
	{
		for (int index = 0; index < current->m_count; ++index)
		{
			Branch* currentBranch = &current->m_branch[index];
			Branch* otherBranch = &other->m_branch[index];

			std::copy(otherBranch->m_rect.m_min,
				otherBranch->m_rect.m_min + NUMDIMS,
				currentBranch->m_rect.m_min);

			std::copy(otherBranch->m_rect.m_max,
				otherBranch->m_rect.m_max + NUMDIMS,
				currentBranch->m_rect.m_max);

			currentBranch->m_child = AllocNode();
			CopyRec(currentBranch->m_child, otherBranch->m_child);
		}
	}
	else // A leaf node
	{
		for (int index = 0; index < current->m_count; ++index)
		{
			Branch* currentBranch = &current->m_branch[index];
			Branch* otherBranch = &other->m_branch[index];

			std::copy(otherBranch->m_rect.m_min,
				otherBranch->m_rect.m_min + NUMDIMS,
				currentBranch->m_rect.m_min);

			std::copy(otherBranch->m_rect.m_max,
				otherBranch->m_rect.m_max + NUMDIMS,
				currentBranch->m_rect.m_max);

			currentBranch->m_data = otherBranch->m_data;
		}
	}
}


RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::Save(const char* a_fileName)
{
	RTFileStreamPartition stream;
	if (!stream.OpenWrite(a_fileName))
	{
		return false;
	}

	bool result = Save(stream);

	stream.Close();

	return result;
}


RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::Save(RTFileStreamPartition& a_stream)
{
	// Write some kind of header
	int dataFileId = ('R' << 0) | ('T' << 8) | ('R' << 16) | ('E' << 24);
	int dataSize = sizeof(DATATYPE);
	int dataNumDims = NUMDIMS;
	int dataElemSize = sizeof(ELEMTYPE);
	int dataElemRealSize = sizeof(ELEMTYPEREAL);
	int dataMaxNodes = TMAXNODES;
	int dataMinNodes = TMINNODES;

	a_stream.Write(dataFileId);
	a_stream.Write(dataSize);
	a_stream.Write(dataNumDims);
	a_stream.Write(dataElemSize);
	a_stream.Write(dataElemRealSize);
	a_stream.Write(dataMaxNodes);
	a_stream.Write(dataMinNodes);

	// Recursively save tree
	bool result = SaveRec(m_root, a_stream);

	return result;
}


RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::SaveRec(Node* a_node, RTFileStreamPartition& a_stream)
{
	a_stream.Write(a_node->m_level);
	a_stream.Write(a_node->m_count);

	if (a_node->IsInternalNode())  // not a leaf node
	{
		for (int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);

			SaveRec(curBranch->m_child, a_stream);
		}
	}
	else // A leaf node
	{
		for (int index = 0; index < a_node->m_count; ++index)
		{
			Branch* curBranch = &a_node->m_branch[index];

			a_stream.WriteArray(curBranch->m_rect.m_min, NUMDIMS);
			a_stream.WriteArray(curBranch->m_rect.m_max, NUMDIMS);

			a_stream.Write(curBranch->m_data);
		}
	}

	return true; // Should do more error checking on I/O operations
}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::RemoveAll()
{
	// Delete all existing nodes
	Reset();

	
	m_root = new Node;// AllocNode();
	m_root->m_level = 0;
}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::Reset()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	// Delete all existing nodes

	//reset_root(m_root);

	RemoveAllRec(m_root);
#else // RTREE_DONT_USE_MEMPOOLS
	// Just reset memory pools.  We are not using complex types
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::reset_root(Node* a_node) {
	if (a_node->IsInternalNode()) // This is an internal node in the tree
	{
		//if (a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer) {
			//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer->clear();
			//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer->~DoublyLinkedList();
			a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = nullptr;
		//}
		//if (a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer) {
			//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer->clear();
			//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer->~DoublyLinkedList();
			a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = nullptr;
		//}
		
		reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);// = nullptr;
		reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);// = nullptr;

		for (int index = 0; index < a_node->m_count; ++index) {
			reset_root(a_node->m_branch[index].m_child);
		}
	}
	//if (a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer) {
		//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer->clear();
		//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer->~DoublyLinkedList();
		a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = nullptr;
	//}
	//if (a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer) {
		//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer->clear();
		//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer->~DoublyLinkedList();
		a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = nullptr;
	//}
	reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);// = nullptr;
	reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);// = nullptr;
}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::RemoveAllRec(Node* a_node)
{
	ASSERT(a_node);
	ASSERT(a_node->m_level >= 0);

	if (a_node->IsInternalNode()) // This is an internal node in the tree
	{
		//reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);// = nullptr;
		//reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);// = nullptr;

		for (int index = 0; index < a_node->m_count; ++index) {
			RemoveAllRec(a_node->m_branch[index].m_child);
		}
	}
	//reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);// = nullptr;
	//reset_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);// = nullptr;
	FreeNode(a_node);
}


RTREE_TEMPLATE
typename RTREE_PARTITION_QUAL::Node* RTREE_PARTITION_QUAL::AllocNode()
{
	Node* newNode;
#ifdef RTREE_DONT_USE_MEMPOOLS
	newNode = new Node;
	//newNode->m_branch = new Branch[TMAXNODES];
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
	InitNode(newNode);
	return newNode;
}

RTREE_TEMPLATE
template<typename T>
inline void RTREE_PARTITION_QUAL::free_shared_ptr_linked_list(shared_ptr<T>& const smart_pointer) {
	if (smart_pointer != nullptr) {
		smart_pointer->clear();
		smart_pointer.reset();
		//a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer = nullptr;
		/* 210709 delete a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer;
		a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer = nullptr;*/
	}
	smart_pointer = nullptr;
}

RTREE_TEMPLATE
template<typename T>
inline void RTREE_PARTITION_QUAL::reset_shared_ptr_linked_list(shared_ptr<T>& const smart_pointer) {
	if (smart_pointer != nullptr) {
		//smart_pointer->clear();
		//smart_pointer.reset();
	}
	smart_pointer = nullptr;
}

//210730
RTREE_TEMPLATE
template<typename T>
void  RTREE_PARTITION_QUAL::free_pointer_vector(vector<T>& const ponter_vector) {
	/*for (auto au : ponter_vector) {
		free_shared_ptr_linked_list(au);
	}*/

	for (int id_pointer = 0; id_pointer < ponter_vector.size(); id_pointer++) {
		free_shared_ptr_linked_list(ponter_vector[id_pointer]);
	}

	ponter_vector.clear();
	ponter_vector.shrink_to_fit();
}

RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::FreeNode(Node* a_node)
{
	ASSERT(a_node);

#ifdef RTREE_DONT_USE_MEMPOOLS

	//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = nullptr;
	//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = nullptr;

	/*=================================================== Rect in Node ========================================================*/
	free_rectangle(a_node->rectangle_in_node);

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~    210718 for rectangle in node     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*--------------------------- for one whole approximation ---------------------------------*/
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = nullptr;
	//a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = nullptr;
	/*free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);*/
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer);
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer);
	/*-----------------------------------------------------------------------------------------*/

	/*-------------- for each segment approximation that has max differnce ---------------------*/
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
	/*------------------------------------------------------------------------------------------*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~ partition 210701 clear pointer of approxiation points ~~~~~~~~~~~~~~~~~~~~~*/
	//free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_entry_struct.list_original_pointer);
	//a_node->rectangle_in_node.rectangle_entry_struct.list_original_pointer.reset();
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_entry_struct.list_original_pointer);
	free_shared_ptr_linked_list(a_node->rectangle_in_node.rectangle_entry_struct.list_partition_pointer);
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*=====================================================================================================================*/

	a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = nullptr;
	/*=================================================== Rect in Branch ===================================================*/
	for (int index = 0; index < a_node->m_count; ++index) {

		//free_rectangle(a_node->m_branch[index].m_rect);//210730

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~    210718 for rectangle in node     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		/*--------------------------- for one whole approximation ---------------------------------*/
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
		//free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);
		//free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer);
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer);
		/*-----------------------------------------------------------------------------------------*/

		/*-------------- for each segment approximation that has max differnce ---------------------*/
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
		/*------------------------------------------------------------------------------------------*/

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~ partition 210701 clear pointer of approxiation points ~~~~~~~~~~~~~~~~~~~~~*/
		//a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer = nullptr;
		// a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer
		//a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
		//a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer = nullptr;
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer);
		free_shared_ptr_linked_list(a_node->m_branch[index].m_rect.rectangle_entry_struct.list_partition_pointer);

		//if (a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer != nullptr) {
		//	a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer->clear();
		//	a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer.reset();
		//	//a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer = nullptr;
		//	/* 210709 delete a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer;
		//	a_node->m_branch[index].m_rect.rectangle_entry_struct.list_original_pointer = nullptr;*/
		//}
		//
		//if (a_node->m_branch[index].m_rect.rectangle_entry_struct.list_partition_pointer != nullptr) {
		//	a_node->m_branch[index].m_rect.rectangle_entry_struct.list_partition_pointer->clear();
		//	a_node->m_branch[index].m_rect.rectangle_entry_struct.list_partition_pointer.reset();
		//	//a_node->m_branch[index].m_rect.rectangle_entry_struct.list_partition_pointer = nullptr;
		//	/* 210709delete a_node->m_branch[index].m_rect.rectangle_entry_struct.list_partition_pointer;
		//	a_node->m_branch[index].m_rect.rectangle_entry_struct.list_partition_pointer = nullptr;*/
		//}
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	}
	/*=====================================================================================================================*/

	delete[] a_node->m_branch;
	a_node->m_branch = nullptr;

	delete a_node;
	a_node = nullptr;

#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


// Allocate space for a node in the list used in DeletRect to
// store Nodes that are too empty.
RTREE_TEMPLATE
typename RTREE_PARTITION_QUAL::ListNode* RTREE_PARTITION_QUAL::AllocListNode()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	return new ListNode;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

RTREE_TEMPLATE
template<typename T>
inline void RTREE_PARTITION_QUAL::allocate_minmax_vector_in_rectangle(T& const a_rect) {//210707 allocatge new memory to vector for rectangle
	a_rect.m_min.assign(NUMDIMS, INF);
	a_rect.m_max.assign(NUMDIMS, -INF);
}

RTREE_TEMPLATE
template<typename T>
inline void RTREE_PARTITION_QUAL::free_minmax_vector_in_rectangle(T& const a_rect) {//210816 free new memory to vector for rectangle
	a_rect.m_min.clear();
	a_rect.m_min.shrink_to_fit();
	a_rect.m_max.clear();
	a_rect.m_max.shrink_to_fit();
}

//************************************
// Method:allocate_rectangle_entry
// Qualifier: 1 allocatge new memory to pointer or vectotr. Because rect was changed from array to pointer. pointer needs new memory. So whenever apply new Branch, needs this frunction for new memory
//            
//              1 is entry rectangle. For approximation & time series.
// date:210707 personal function
// author:
//************************************
RTREE_TEMPLATE
template<typename T>
inline T& RTREE_PARTITION_QUAL::allocate_rectangle_entry(T& const a_rect) { //210804 allocatge new memory to pointer or vector for entry

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~    210718 for rectangle in branch     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//a_rect.original_time_series_vector_pointer.assign(1024, INF);
	//a_rect.reconstruct_time_series_vector_pointer.assign(1024, INF);//210718

	a_rect.rectangle_entry_struct.original_time_series_vector_pointer = std::make_shared<vector<long double>>();
	//a_rect.original_time_series_vector_pointer->assign(1024, INF);
	//a_rect.reconstruct_time_series_vector_pointer = std::make_shared<vector<long double>>();
	//a_rect.reconstruct_time_series_vector_pointer->assign(1024, INF);

	//a_rect.linked_list_original_vector.resize(TMAXNODES, DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>());
	//a_rect.linked_list_partition_vector.resize(TMAXNODES, DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>());
	a_rect.rectangle_entry_struct.list_original_pointer = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();// = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
	a_rect.rectangle_entry_struct.list_partition_pointer = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
	/*a_rect.rectangle_entry_struct.list_original_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	a_rect.rectangle_entry_struct.list_partition_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;*/
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	return a_rect;
}


//************************************
// Method:allocate_rectangle_node
// Qualifier: 1 allocatge new memory to pointer or vectotr. Because rect was changed from array to pointer. pointer needs new memory. So whenever apply new Branch, needs this frunction for new memory
//            2 0 is node rectangle. For cover convex hull. pointers vector. 
//             
// date:210707 personal function
// author:
//************************************
RTREE_TEMPLATE
template<typename T>
inline T& RTREE_PARTITION_QUAL::allocate_rectangle_node(T& const a_rect) { //210804 allocatge new memory to pointer or vector for node:
	
	if (option_tree_struct.type_tree == 2) {

		a_rect.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;// max whole difference
		a_rect.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = -INF;// max whole difference square root
		/*--------------------------- for one whole approximation ---------------------------------*/
		a_rect.seed_0 = (std::numeric_limits<int>::min)();
		a_rect.seed_1 = (std::numeric_limits<int>::min)();
		a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
		a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
		//a_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer;// = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
		//a_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer;// = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();

		a_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = std::make_shared<vector<long double>>();
		a_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = std::make_shared<vector<long double>>();
		/*-----------------------------------------------------------------------------------------*/

		return a_rect;
	}

	/*~~~~~~~~~~~~~~~~~~~   210716 Volume   ~~~~~~~~~~~~~~~~~~~~~~~~*/
	a_rect.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment = -INF;// max semgent difference
	a_rect.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = -INF;// max segment difference square root
	a_rect.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point = -INF;// max reconstruct point difference
	a_rect.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = -INF; // max reconstruct point difference square root
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~    210718 for rectangle in node     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*--------------------------- max(at+b)-min(at+b) ---------------------------------*/
	//a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	//a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, -INF);
	/*----------------------------------------------------------------------------------*/

	/*-------------- for each segment approximation that has max differnce ---------------------*/
	a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer = std::make_shared<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>();
	/*------------------------------------------------------------------------------------------*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	return a_rect;
}

//************************************
// Method:allocate_rectangle
// Qualifier: 1 allocatge new memory to pointer or vectotr. Because rect was changed from array to pointer. pointer needs new memory. So whenever apply new Branch, needs this frunction for new memory
//            2 0 is node rectangle. For cover convex hull. pointers vector. 
//              1 is entry rectangle. For approximation & time series.
// date:210707 personal function
// author:
//************************************
RTREE_TEMPLATE
template<typename T, typename Y>
inline Y& RTREE_PARTITION_QUAL::allocate_rectangle(const T& const option_rectangle, Y& const a_rect) {

	/*###############################   option_tree_struct   #######################################*/
	switch (option_tree_struct.type_tree) {

	case 0: {// R-tree
		assert(0);
		break;
	}
	case 1:{// distance tree
		/*============= partition 210705 ===================================*/
		switch (option_rectangle) {
		case 0: {// 0 is node rectangle. For cover convex hull. pointers vector, 
			allocate_minmax_vector_in_rectangle(a_rect);
			allocate_rectangle_node(a_rect);
			break;
		}
		case 1: {//1 is entry rectangle. For approximation & time series
			allocate_minmax_vector_in_rectangle(a_rect);
			allocate_rectangle_entry(a_rect);
			break;
		}
		case 2://no make_shared for pointer or assign for vector
			allocate_minmax_vector_in_rectangle(a_rect);
			break;
		default:
			assert(0);
		}
		/*=====================================================================*/
		break;
	}
	case 2: {//211008 distance tree with SAPLA_MBR
		/*============= SAPLA_MBR 211008 ===================================*/
		switch (option_rectangle) {
		case 0: {// 0 is node rectangle. For cover SAPLA_MBR. pointers vector,
			//allocate_minmax_vector_in_rectangle(a_rect);
			allocate_rectangle_node(a_rect);
			break;
		}
		case 1: {//1 is entry rectangle. For approximation & time series
			//allocate_minmax_vector_in_rectangle(a_rect);
			allocate_rectangle_entry(a_rect);
			break;
		}
		case 2:// redundant no make_shared for pointer or assign for vector
			//allocate_minmax_vector_in_rectangle(a_rect);
			break;
		default:
			assert(0);
		}
		/*=====================================================================*/
		break;
	}
	default:
		assert(0);
	}
	/*##################################################################################*/

	return a_rect;
}

//21802 No free memory, only clear vector.
RTREE_TEMPLATE
template<typename T>
inline T& RTREE_PARTITION_QUAL::clear_vector_in_rectangle(T& const a_rect) {

	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    210718 for rectangle in node     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	/*~~~~~~~~~~~~~~~~~~~   210716 Volume   ~~~~~~~~~~~~~~~~~~~~~~~~*/
	a_rect.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment = -INF;// max semgent difference
	a_rect.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = -INF;// max segment difference square root
	a_rect.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;// max whole difference
	a_rect.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = -INF;// max whole difference square root

	a_rect.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point = -INF;// max reconstruct point difference
	a_rect.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = -INF; // max reconstruct point difference square root
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*--------------------------- max(at+b)-min(at+b) ---------------------------------*/
	/*a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.shrink_to_fit();
	a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.shrink_to_fit();*/
	/*----------------------------------------------------------------------------------*/

	/*----- for one whole approximation ----*/
	a_rect.seed_0 = (std::numeric_limits<int>::min)();
	a_rect.seed_1 = (std::numeric_limits<int>::min)();
	/*--------------------------------------*/

	/*--------------------------- for one whole approximation ---------------------------------*/
	a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
	a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
	/*free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);*/
	/*-----------------------------------------------------------------------------------------*/

	/*-------------- for each segment approximation that has max differnce ---------------------*/
	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.shrink_to_fit();

	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.shrink_to_fit();

	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.shrink_to_fit();

	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.shrink_to_fit();

	a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->clear();
	/*free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);*/
	/*------------------------------------------------------------------------------------------*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~         vector of vector / pointer       ~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//for (auto&& au: a_rect.original_time_series_vector_vector) {
	//	au.clear();
	//	au.shrink_to_fit();
	//}
	//a_rect.original_time_series_vector_vector.clear();//210706
	//a_rect.original_time_series_vector_vector.shrink_to_fit();

	/*---------------    210710 shared pointer    ------------------------*/
	a_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.clear();
	a_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.shrink_to_fit();
	//a_rect.reconstruct_time_series_pointer_vector.clear();
	//a_rect.reconstruct_time_series_pointer_vector.shrink_to_fit();
	//a_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.clear();
	//a_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.shrink_to_fit();
	a_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.clear();
	a_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.shrink_to_fit();

	/*free_pointer_vector(a_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector);
	free_pointer_vector(a_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector);
	free_pointer_vector(a_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector);*/
	/*--------------------------------------------------------------------*/
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	return a_rect;
}

RTREE_TEMPLATE
template<typename T>
inline T& RTREE_PARTITION_QUAL::free_rectangle(T& const a_rect) {//210730 free memory to pointer or vectotr

	/*if (a_rect.m_ min != nullptr) {
		delete[] a_rect.m_ min;
		a_rect.m_min = nullptr;
	}

	if (a_rect.m_max != nullptr) {
		delete[] a_rect.m_max;
		a_rect.m_max = nullptr;
	}*/

	//free_minmax_vector_in_rectangle(a_rect);

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&           partition 210705          &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~    210718 for rectangle in branch     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//a_rect.rectangle_entry_struct.original_time_series_vector_pointer->reset();
	if (a_rect.rectangle_entry_struct.original_time_series_vector_pointer) {
		a_rect.rectangle_entry_struct.original_time_series_vector_pointer->clear();
		a_rect.rectangle_entry_struct.original_time_series_vector_pointer->shrink_to_fit();
		free_shared_ptr_linked_list(a_rect.rectangle_entry_struct.original_time_series_vector_pointer);
	}

	//if (a_rect.reconstruct_time_series_vector_pointer) {
	//	a_rect.reconstruct_time_series_vector_pointer->clear();//210718
	//	a_rect.reconstruct_time_series_vector_pointer->shrink_to_fit();
	//	free_shared_ptr_linked_list(a_rect.reconstruct_time_series_vector_pointer);
	//}
	//a_rect.rectangle_entry_struct.list_original_pointer.reset();
	//a_rect.rectangle_entry_struct.list_original_pointer = nullptr;
	free_shared_ptr_linked_list(a_rect.rectangle_entry_struct.list_original_pointer);
	free_shared_ptr_linked_list(a_rect.rectangle_entry_struct.list_partition_pointer);
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        210718 for rectangle in node        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	/*~~~~~~~~~~~~~~~~~~~   210716 Volume   ~~~~~~~~~~~~~~~~~~~~~~~~*/
	a_rect.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment = -INF;// max semgent difference
	a_rect.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = -INF;// max segment difference square root
	a_rect.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;// max whole difference
	a_rect.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = -INF;// max whole difference square root

	a_rect.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point = -INF;// max reconstruct point difference
	a_rect.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = -INF; // max reconstruct point difference square root
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*--------------------------- max(at+b)-min(at+b) ---------------------------------*/
	a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.shrink_to_fit();
	a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.shrink_to_fit();
	/*----------------------------------------------------------------------------------*/

	/*----- for one whole approximation ----*/
	a_rect.seed_0 = (std::numeric_limits<int>::min)();
	a_rect.seed_1 = (std::numeric_limits<int>::min)();
	/*--------------------------------------*/

	/*--------------------------- for one whole approximation ---------------------------------*/
	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	//a_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = nullptr;
	//a_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = nullptr;
	/*a_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer->reset();
	a_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer->reset();*/
	//free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);
	//free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);
	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer);
	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer);
	/*-----------------------------------------------------------------------------------------*/

	/*-------------- for each segment approximation that has max differnce ---------------------*/
	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.shrink_to_fit();

	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.shrink_to_fit();

	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.shrink_to_fit();

	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.clear();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.shrink_to_fit();

	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
	free_shared_ptr_linked_list(a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
	/*------------------------------------------------------------------------------------------*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~         vector of vector / pointer       ~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//for (auto&& au: a_rect.original_time_series_vector_vector) {
	//	au.clear();
	//	au.shrink_to_fit();
	//}
	//a_rect.original_time_series_vector_vector.clear();//210706
	//a_rect.original_time_series_vector_vector.shrink_to_fit();

	/*---------------    210710 shared pointer    ------------------------*/
	free_pointer_vector(a_rect.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector);
	//free_pointer_vector(a_rect.reconstruct_time_series_pointer_vector);
	assert(a_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.empty());
	//free_pointer_vector(a_rect.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector);
	/*for (int id_pointer = 0; id_pointer < ponter_vector.size(); id_pointer++) {
		free_shared_ptr_linked_list(ponter_vector[id_pointer]);
	}*/
	free_pointer_vector(a_rect.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector);
	/*--------------------------------------------------------------------*/
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

	return a_rect;

}

RTREE_TEMPLATE
template<typename T>//210804  For rectangle in struct PartitionVars
inline T& RTREE_PARTITION_QUAL::clear_rectangle_PartitionVars(T& const a_rect) {
	clear_vector_in_rectangle(a_rect);
	a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer.reset();
	a_rect.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer.reset();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer.reset();
	a_rect.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer.reset();

	return a_rect;
}

//************************************
// Method:AllocBranch
// Qualifier: Because rect was changed from array to pointer. pointer needs new memory. So whenever apply new Branch, needs this frunction for new memory
// date: personal function
// author:
//************************************
RTREE_TEMPLATE
template<typename T>
inline void RTREE_PARTITION_QUAL::AllocBranch(T& const a_branch) {
	assert(0);
	//allocate_rectangle(a_branch.m_rect);

	//a_branch.m_rect.m_ min = new ELEMTYPE[NUMDIMS];
	//a_branch.m_rect.m_ max = new ELEMTYPE[NUMDIMS];
	///*=================================== partition 210705 =========================================================*/
	//a_branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	//a_branch.m_rect.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, INF);
	//a_branch.m_rect.original_time_series_vector_pointer.assign(1024, INF);
	//a_branch.m_rect.linked_list_original_vector.resize(TMAXNODES, DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>());
	//a_branch.m_rect.linked_list_partition_vector.resize(TMAXNODES, DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>());
	//a_branch.m_rect.rectangle_entry_struct.list_original_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	//a_branch.m_rect.rectangle_entry_struct.list_partition_pointer = new DoublyLinkedList<AREA_COEFFICIENT_RTREE_NODE>;
	///*===========================================================================================================*/
}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::FreeListNode(ListNode* a_listNode)
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	delete a_listNode;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::InitNode(Node* a_node)
{
	a_node->m_count = 0;
	a_node->m_level = -1;
	a_node->m_branch = new Branch[TMAXNODES];

	/*++++++++++++++++++   partition 210618   +++++++++++++++++++++*/
	for (int id_branch = 0; id_branch < TMAXNODES; id_branch++) {
		a_node->m_branch[id_branch].branch_to_parent_node_pointer = a_node;
	}
	/*---------210709 rectangle in Node-------------*/
	allocate_rectangle(0, a_node->rectangle_in_node);
	/*----------------------------------------------*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
}

RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::InitRect(Rect* a_rect)
{

	for (int index = 0; index < NUMDIMS; ++index)
	{
		a_rect->m_min[index] = (ELEMTYPE)0;
		a_rect->m_max[index] = (ELEMTYPE)0;
	}
}

// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// Used in  InsertRect(), InsertRectRec()
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::InsertRectRec(Branch& a_branch, Node* a_node, Node** a_newNode, int a_level) {
	
	/*....................................................................................*/
#ifdef _DEBUG
	ASSERT(a_node && a_newNode);
	ASSERT(a_level >= 0 && a_level <= a_node->m_level);
	/*~~~~~~~~~~~~~~~~~  partition 210618 ~~~~~~~~~~~~~~~~~~~~~~~~~*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		assert(is_rectangle_entry(a_branch.m_rect));
		//update_branch_pointer_vector(*a_node);
		///update_branch_pointer_vector(**a_newNode);
		break;
	}
	default:
		break;
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#endif
	/*....................................................................................*/

	// recurse until we reach the correct level for the new record. data records, will always be called with a_level == 0 (leaf)
	if (a_node->m_level > a_level) {//a_node is internal node,  In insert function a_level == 0
		// Still above level for insertion, go down tree recursively
		Node* otherNode = nullptr;// = new Node;

		//memory_account[5] += sizeof(otherNode);
		//otherNode->m_branch = new Branch[TMAXNODES];

		// find the optimal branch in internal node for this record. a_branch is entry or node?
		int index = PickBranch(&a_branch.m_rect, a_node);

		//cout <<"a_node->m_level : " <<a_node->m_level << " Pick Branch :" << index << " a_level : "<< a_level << endl;

		// recursively insert this record into the picked branch
		bool childWasSplit = InsertRectRec(a_branch, a_node->m_branch[index].m_child, &otherNode, a_level);

		if (!childWasSplit) {//Not split. a_node is internal node 
			// Child was not split. Merge the bounding box of the new record with the existing bounding box
			//cout << "Child was not splited : \n";
			//a_node->m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));

			/*~~~~~~~~~~~~~~~~~  partition 210618 ~~~~~~~~~~~~~~~~~~~~~~~~~*/
			switch (representation_option) {
			case 2:// PLA
			case 3://APCA
			case 4://PAA
			case 7://PAALM
			case 9://SAX
			case 5://CHEBY
			case 6: // ICDE07
			case 8: {// Initial 200706
				switch (option_tree_struct.type_tree) {
				case 1: {// distance R-tree
					/*~~~~~~~~~~~210815  ~~~~~~~~~~~~*/
					assign_node_rectangle_to_branch_rectangle(*a_node->m_branch[index].m_child, a_node->m_branch[index].m_rect);
					/*------------------------------- Original --------------------------------------------------------*/
					CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect), &(a_node->m_branch[index].m_rect));
					/*-------------------------------------------------------------------------------------------------*/
					a_node->m_branch[index].m_child->rectangle_in_node.m_min = a_node->m_branch[index].m_rect.m_min;
					a_node->m_branch[index].m_child->rectangle_in_node.m_max = a_node->m_branch[index].m_rect.m_max;

					copy_rect_original_series_approximation_vector_pointer(a_branch.m_rect, a_node->rectangle_in_node);
					get_convex_hull_by_rectangle(a_node->rectangle_in_node);

					/*.......................................................................*/
#ifdef _DEBUG
					assert_node(*a_node);
#endif
					/*.......................................................................*/
					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
					break;
				}
				case 2: {// convex hull
					if (a_node->m_branch[index].m_child->IsLeaf()) { // Compute convex hull for Leaf node
						get_convex_hull_volume_by_distance_approximation_global(a_node->m_branch[index].m_child->m_branch, a_node->m_branch[index].m_child->m_count,a_node->m_branch[index].m_child->rectangle_in_node);
					}
					else {//a_node->m_branch[index].m_child is internal node after add new branch
					
						get_convex_hull_volume_by_distance_internal_node_global(a_node->m_branch[index].m_child->m_branch, a_node->m_branch[index].m_child->m_count, a_node->m_branch[index].m_child->rectangle_in_node);
						/*...................*/
#ifdef _DEBUG
						//assert_rectangle_global_minmax(a_branch.m_child->rectangle_in_node);
						assert_rectangle_global_minmax(a_node->m_branch[index].m_child->rectangle_in_node);
#endif
						/*...................*/
						
					}
					break;
				}
				default:
					assert(0);
					break;
				}
				break;
			}
			default: {
				assert(0);
				/*------------------------------- Original --------------------------------------------------------*/
				CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect), &(a_node->m_branch[index].m_rect));
				/*-------------------------------------------------------------------------------------------------*/
				break;
			}
			}
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		
			return false;
		}
		else {// Splited a_node->m_branch[index].m_child, otherNode
			// Child was split. The old branches are now re-partitioned to two nodes so we have to re-calculate the bounding boxes of each node
			//cout << "Child was splited : \n";
			//a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
			Branch branch;
			branch.m_child = otherNode;

			/*--------210814---------*/
			otherNode->node_to_parent_branch_pointer = &branch;
			/*----------------------*/

			/*~~~~~~~~~~~~~~~~~  partition 210618 ~~~~~~~~~~~~~~~~~~~~~~~~~*/
			switch (representation_option) {
			case 2: // PLA
			case 3://APCA
			case 4://PAA
			case 7://PAALM
			case 9://SAX
			case 5://CHEBY
			case 6: //ICDE07
			case 8: {// Initial 200706

				switch (option_tree_struct.type_tree) {
				case 1: {// distance R-tree Compute new MBR
					NodeCover(a_node->m_branch[index].m_child, a_node->m_branch[index].m_rect);// min max MBR
					/*~~~~~~~~~~~210814 a_node is internal node and a_node->m_branch[index].m_child was split otherNode  ~~~~~~~~~~~~*/
					update_rectangle_from_splited_rectangles(a_node->m_branch[index].m_rect, otherNode->rectangle_in_node, a_node->rectangle_in_node);
					if (a_node->node_to_parent_branch_pointer) {
						a_node->node_to_parent_branch_pointer->m_rect = a_node->rectangle_in_node;
					}
					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
					/*++++++++++++++++  210814  +++++++++++++++++++++++*/
					allocate_rectangle(2, branch.m_rect);
					//AllocBranch(branch);
					//branch.m_rect = NodeCover(otherNode);
					NodeCover(otherNode, branch.m_rect);
					/*++++++++++++++++++++++++++++++++++++++++++++++++*/
					
					break;
				}
				case 2: {// Convex Tree slited new nodes already have minmax global

					if (a_node->m_branch[index].m_child->IsLeaf()) {
						//break;
					}
					else {
						get_convex_hull_volume_by_distance_internal_node_global(a_node->m_branch[index].m_child->m_branch, a_node->m_branch[index].m_child->m_count, a_node->m_branch[index].m_child->rectangle_in_node);
						get_convex_hull_volume_by_distance_internal_node_global(otherNode->m_branch, otherNode->m_count, otherNode->rectangle_in_node);
					}
					/*.......................................................................*/
#ifdef _DEBUG
					assert_rectangle_global_minmax(a_node->m_branch[index].m_child->rectangle_in_node);
					assert_rectangle_global_minmax(otherNode->rectangle_in_node);
#endif
					/*.......................................................................*/
					break;
				}
				default:
					assert(0);
					break;
				}

				
				break;
			}
			default: {
				assert(0);
				break;
			}
			}
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			//branch-> otherNode. otherNode is sub node of a_node
			// The old node is already a child of a_node. Now add the newly-created node to a_node as well. a_node might be split because of that.
			return AddBranch(&branch, a_node, a_newNode);
		}
	}
	else if (a_node->m_level == a_level) {// Ad branch entry in leaf node.  this is leaf node->contain entry.  const a_level == 0, a_node->m_level == 0 is leaf node.

	/*....................................................................................*/
#ifdef _DEBUG
	assert(a_level == 0);
#endif
	/*....................................................................................*/

		// We have reached level for insertion. Add rect, split if necessary
	    // return false;//0  a_node not split, return true; 1 a_node is splited
		return AddBranch(&a_branch, a_node, a_newNode);
	}
	else {
		// Should never occur
		ASSERT(0);
		return false;
	}
}


// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1(true) if root was split, 0(false) if it was not.
// The level argument specifies the number of steps up from the leaf level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
// used in RemoveRect(), 4 Insert()s
// In Insert function a_level == 0
RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::InsertRect(Branch& a_branch, Node** a_root, int a_level) {// a_level == 0
	/*.......................................................................*/
#ifdef _DEBUG
	ASSERT(a_root);
	ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);
	if (option_tree_struct.type_tree == 1) {
		for (int index = 0; index < NUMDIMS; ++index) {
			ASSERT(a_branch.m_rect.m_min[index] <= a_branch.m_rect.m_max[index]);
		}
	}
	/*~~~~~~~~~~~~~~~~~  partition 210618 ~~~~~~~~~~~~~~~~~~~~~~~~~*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6://ICDE07
	case 8: {// Initial 200706
		//update_branch_pointer_vector(**a_root);
		break;
	}
	default:
		break;
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#endif
	/*.......................................................................*/

	//cout <<"a_branch id : "<< a_branch.m_data<< ", a_lelve : " << a_level << endl;

	Node* newNode = nullptr;// = new Node;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!maybe influce result
	//newNode->m_branch = new Branch[TMAXNODES];

	//memory_account[6] += sizeof(newNode->m_branch)*TMAXNODES;

	// In Insert function, a_level == 0
	// return false a_root is not split, return true is split
	if (InsertRectRec(a_branch, *a_root, &newNode, a_level)) {
		/*&&&&   a_root was split to newNode, Need a new root node &&&&*/
		
		/*....................................................................*/
#ifdef _DEBUG
		assert_node(**a_root);
		assert_node(*newNode);
#endif
		/*....................................................................*/
		switch (option_tree_struct.type_tree) {
		case 1: {
			// Grow tree taller and new root. newRoot point to a_root and newNode
			Node* newRoot = AllocNode();
			//pointer[2] = newRoot;
			newRoot->m_level = (*a_root)->m_level + 1;

			/*~~~~ rectangle in branch cover a_root ~~~~~ */
			Branch branch;
			allocate_rectangle(2, branch.m_rect);
			//AllocBranch(branch);
			// add old root node as a child of the new root, rectangle in branch cover a_root.
			// branch.m_rect = (*a_root);
			NodeCover(*a_root, branch.m_rect);
			branch.m_child = *a_root;
			/*------------210814-----------*/
			(*a_root)->node_to_parent_branch_pointer = &branch;
			/*-----------------------------*/
			AddBranch(&branch, newRoot, NULL);
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

			/*~~~~ rectangle in branch cover new Node ~~~~~ */
			// branch.m_rect = NodeCover(newNode);
			Branch branch0;
			allocate_rectangle(2, branch0.m_rect);
			// AllocBranch(branch0);
			// add the split node as a child of the new root, rectangle in branch0 cover newNode.
			NodeCover(newNode, branch0.m_rect);
			branch0.m_child = newNode;
			/*------------210814-----------*/
			newNode->node_to_parent_branch_pointer = &branch0;
			/*-----------------------------*/
			AddBranch(&branch0, newRoot, NULL);
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

			// set the new root as the root node
			*a_root = newRoot;
			//cout << "true: Root not split" << endl;
			break;
		}
		case 2: {
			if ((*a_root)->IsLeaf()) {
				*a_root = (*a_root)->node_to_parent_branch_pointer->branch_to_parent_node_pointer;
			}
			else {
				// Grow tree taller and new root. newRoot point to a_root and newNode
				Node* newRoot = new Node;
				//pointer[2] = newRoot;
				newRoot->m_level = (*a_root)->m_level + 1;
				newRoot->m_count = 0;
				newRoot->m_branch = new Branch[TMAXNODES];

				
				for (int id_branch = 0; id_branch < TMAXNODES; id_branch++) {
					newRoot->m_branch[id_branch].branch_to_parent_node_pointer = newRoot;
				}
				/*---------210709 rectangle in Node-------------*/
				allocate_rectangle(0, newRoot->rectangle_in_node);
				/*----------------------------------------------*/

				/*~~~~ rectangle in branch cover a_root ~~~~~ */
				get_convex_hull_volume_by_distance_internal_node_global((*a_root)->m_branch, (*a_root)->m_count, (*a_root)->rectangle_in_node);
				newRoot->m_branch[newRoot->m_count].m_child = (*a_root);
				//node_sup_new->m_branch[node_sup_new->m_count].node_rect_pointer = make_shared<Rect>(a_node->rectangle_in_node);
				/*------------210814-----------*/
				(*a_root)->node_to_parent_branch_pointer = &newRoot->m_branch[newRoot->m_count];
				/*-----------------------------*/
				newRoot->m_branch[newRoot->m_count].branch_to_parent_node_pointer = newRoot;

				//node_sup_new->m_branch[node_sup_new->m_count] = branch_0;

				++newRoot->m_count;
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

				/*~~~~ rectangle in branch cover new Node ~~~~~ */
				get_convex_hull_volume_by_distance_internal_node_global(newNode->m_branch, newNode->m_count, newNode->rectangle_in_node);
				newRoot->m_branch[newRoot->m_count].m_child = newNode;
				//newRoot->m_branch[newRoot->m_count].node_rect_pointer = make_shared<Rect>((*a_newNode)->rectangle_in_node);
				/*------------210814-----------*/
				newNode->node_to_parent_branch_pointer = &newRoot->m_branch[newRoot->m_count];
				/*-----------------------------*/
				newRoot->m_branch[newRoot->m_count].branch_to_parent_node_pointer = newRoot;

				//newRoot->m_branch[newRoot->m_count] = branch_1;

				++newRoot->m_count;

				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

				*a_root = newRoot;
			}

			break;
		}
		default:
			assert(0);
			break;
		}


		
		return true;
	}

	//cout << "true: Root was split" << endl;
	return false;
}


// Find the smallest rectangle that includes all rectangles in branches of a node.
//RTREE_TEMPLATE
//typename RTREE_PARTITION_QUAL::Rect& RTREE_PARTITION_QUAL::NodeCover(Node* a_node)
//{
//	ASSERT(a_node);
//
//	Rect rect = a_node->m_branch[0].m_rect;                                          //: What the index=0 mean?????
//	for (int index = 1; index < a_node->m_count; ++index)
//	{
//		//cout << "CombineRect()" << endl;
//		rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
//	}
//	//cout << "rect: "<< rect.m_min[0]  << endl;
//	return rect;
//}

// Find the smallest rectangle that includes all rectangles in branches of a node.
// a_rect: is rectangle in branch that pointer to a_node
RTREE_TEMPLATE
typename RTREE_PARTITION_QUAL::Rect& RTREE_PARTITION_QUAL::NodeCover(Node* a_node, Rect& a_rect) {

	ASSERT(a_node);
	ASSERT(!a_rect.m_min.empty() && !a_rect.m_max.empty());

	/*==============        210802     =========================*/
	assign_node_rectangle_to_branch_rectangle(*a_node, a_rect);
	/*==========================================================*/

	//assertRect(a_rect);

	/*if (nullptr == a_rect.m _min&&nullptr == a_rect.m_max) {
		a_rect.m _min= new ELEMTYPE[NUMDIMS];
		a_rect.m_max = new ELEMTYPE[NUMDIMS];
	}*/

	//a_rect = a_node->m_branch[0].m_rect;      //: What the index=0 mean?????
	copy_n(a_node->m_branch[0].m_rect.m_min.begin(), NUMDIMS, a_rect.m_min.begin());
	copy_n(a_node->m_branch[0].m_rect.m_max.begin(), NUMDIMS, a_rect.m_max.begin());
	//	Rect rect = a_node->m_branch[0].m_rect;

	for (int index = 1; index < a_node->m_count; ++index) {
		//cout << "CombineRect()" << endl;
		//a_rect = CombineRect(&a_rect, &(a_node->m_branch[index].m_rect));
		CombineRect(&a_rect, &a_node->m_branch[index].m_rect, &a_rect);
	}
	//cout << "rect: "<< rect.m _min[0]  << endl;

	/*==============        210802     =========================*/
	a_node->rectangle_in_node.m_min = a_rect.m_min;
	a_node->rectangle_in_node.m_max = a_rect.m_max;
	/*==========================================================*/

	return a_rect;
}


// Add a branch to a node.  Split the node if necessary.
// Returns 0 false if node not split.  Old node updated.
// Returns 1 true if node split, sets *new_node to address of new node.
// Used in InsertRect(), LoadNodes(), InsertRectRec()
// Old node updated, becomes one of two.
RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode) {

	/*................................*/
#ifdef _DEBUG
	ASSERT(a_branch);
	ASSERT(a_node);
	if (representation_option == 8) {
		//assert_node(*a_node);
		assert_branch(*a_branch);
	}
	//update_branch_pointer_vector(*a_node);
	//update_branch_pointer_vector(**a_newNode);
	assert_option_tree(option_tree_struct);
#endif
	/*................................*/

	if (a_node->m_count < TMAXNODES) {  // Split a_node won't be necessary, a_newNode is nullptr

		/*..............................................................*/
#ifdef _DEBUG
		const Branch& test_branch0 = a_node->m_branch[a_node->m_count];
#endif
		/*...............................................................*/

		// parent node connect sub new branch
		a_node->m_branch[a_node->m_count] = *a_branch;
		++a_node->m_count;
		a_branch->m_rect.rectangle_entry_struct.list_original_pointer = nullptr;
		/*++++++++++++++++++++++  partition 210618 +++++++++++++++++++*/
		switch (representation_option) {
		case 2: //PLA
		case 3://APCA
		case 4://PAA
		case 7://PAALM
		case 9://SAX
		case 5://CHEBY
		case 6://ICDE07
		case 8: {// Initial 200706

			switch (option_tree_struct.type_tree) {
			case 1: {// distance R-tree

				/*--------- node point to parent branch, branch to parent node ----------*/
				a_node->m_branch[a_node->m_count - 1].branch_to_parent_node_pointer = a_node;// sub new branch connect to its parent node

				if (a_node->m_branch[a_node->m_count - 1].m_child) {// sub node of sub branch connect to sub branch
					a_node->m_branch[a_node->m_count - 1].m_child->node_to_parent_branch_pointer = &a_node->m_branch[a_node->m_count - 1];
				}
				/*-----------------------------------------------------------------------*/

				/*-------------    210711 node contain all sub rectanles, copy them of new branch to parent node   ----------------------------*/
				//a_node->rectangle_in_node.linked_list_partition_vector.clear();
				/*210711 a_node is Leaf Node or Internal Node*/
				copy_rect_original_series_approximation_vector_pointer(a_node->m_branch[a_node->m_count - 1].m_rect, a_node->rectangle_in_node);
				/*-----------------------------------------------------------------------------------------------------------------------------*/

				// partiton all entry with new in leaf node
				if (a_node->m_count > 1 && a_node->IsLeaf()) {

					for (int id_branch = 0; id_branch < a_node->m_count - 1; id_branch++) {
						APLA::get_same_partition_SAPLA_limit_only_move_1(*a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, option_tree_struct);
						//APLA::get_same_partition_SAPLA_210818(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer, number_points);
						//APLA::get_same_partition_SAPLA_debug(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer);

						/*--------------------210817 So far not need reconstructed time series---------*/
						//APLA::getAPLAReconstructSeries(*a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.reconstruct_time_series_vector_pointer);
						//APLA::getAPLAReconstructSeries(*a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.reconstruct_time_series_vector_pointer);
						/*------------------------------------------------------------------------------*/
						//partition_rect_volume(a_node->m_branch[id_branch].m_rect, a_node->rectangle_in_node);
					}

					/*...................*/
#ifdef _DEBUG
					assert_node(*a_node);
#endif
					/*...................*/

				}

				//get_convex_hull_by_node(*a_node);
				//partition_rect_volume(a_node->m_branch[a_node->m_count - 1].m_rect, a_node->rectangle_in_node);

				/*210711*/
				//a_node->rectangle_in_node.linked_list_partition_vector.emplace_back(*a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer);
				break;
			}
			case 2: {// convex hull

				/*--------- node point to parent branch, branch to parent node ----------*/
				a_node->m_branch[a_node->m_count - 1].branch_to_parent_node_pointer = a_node;// sub new branch connect to its parent node

				if (a_node->m_branch[a_node->m_count - 1].m_child) {// sub node of sub branch connect to sub branch
					a_node->m_branch[a_node->m_count - 1].m_child->node_to_parent_branch_pointer = &a_node->m_branch[a_node->m_count - 1];
				}
				/*------------------------------------------------------------------------*/

				/*-------------    210711 node contain all sub rectanles, copy them of new branch to parent node   ----------------------------*/
				//a_node->rectangle_in_node.linked_list_partition_vector.clear();
				/*210711 a_node is Leaf Node or Internal Node*/
				//copy_rect_original_series_approximation_vector_pointer(a_node->m_branch[a_node->m_count - 1].m_rect, a_node->rectangle_in_node);
				/*-----------------------------------------------------------------------------------------------------------------------------*/

				// partiton all entry with new entry in leaf node
				if (a_node->m_count > 1 && a_node->IsLeaf()) {

					for (int id_branch = 0; id_branch < a_node->m_count - 1; id_branch++) {
						APLA::get_same_partition_SAPLA_limit_only_move_1(*a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, option_tree_struct);
						//APLA::get_same_partition_SAPLA_210818(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer, number_points);
						//APLA::get_same_partition_SAPLA_debug(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer);
						//partition_rect_volume(a_node->m_branch[id_branch].m_rect, a_node->rectangle_in_node);
					}

					/*...................*/
#ifdef _DEBUG
					assert_node(*a_node);
#endif
					/*...................*/

				}
				else if (a_node->IsInternalNode()) {//Leaf Node of sup node
					/*...................*/
#ifdef _DEBUG
					assert_rectangle_global_minmax(a_branch->m_child->rectangle_in_node);
#endif
					/*...................*/
					if (!a_node->node_to_parent_branch_pointer) {//Root
						break;
					}
					else {// a_node is internal node, not Root
						break;
					}
					
				}
				//get_convex_hull_by_node(*a_node);// 211225 compute convex hull could put in split function 
				//partition_rect_volume(a_node->m_branch[a_node->m_count - 1].m_rect, a_node->rectangle_in_node);

				/*210711*/
				//a_node->rectangle_in_node.linked_list_partition_vector.emplace_back(*a_node->m_branch[a_node->m_count - 1].m_rect.rectangle_entry_struct.list_partition_pointer);
				break;
			}
			default:
				assert(0);
				break;
			}
			
			break;
		}
		default:
			assert(0);
			break;
		}
		/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

		return false;//0 not split
	}
	else {//  Split. a_node needs split

	     /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     Split      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
		/*................................*/
#ifdef _DEBUG
		ASSERT(a_newNode);
		assert(a_node->m_count == TMAXNODES);
		if (representation_option == 8 && *a_newNode != nullptr) {
			assert_node(**a_newNode);
		}
#endif
		/*................................*/

		/*++++++++++++++++++++++  partition 210706 a_node->IsLeaf() +++++++++++++++++++*/
		// before split node, parition new branches and all sub branches in old node. Thus they all have same partitions
		switch (representation_option) {
		case 2: // PLA
		case 3://APCA
		case 4://PAA
		case 7://PAALM
		case 9://SAX
		case 5://CHEBY
		case 6://ICDE07
		case 8: {// Initial 200706
			if (a_node->IsLeaf()) {

				switch (option_tree_struct.type_tree) {
				case 1: {

					/*......................................*/
#ifdef _DEBUG
					assert(option_tree_struct.type_representation == representation_option);
					assert(a_branch->m_rect.rectangle_entry_struct.list_partition_pointer);
#endif
					/*......................................*/
					for (int id_branch = 0; id_branch < a_node->m_count; id_branch++) {
						APLA::get_same_partition_SAPLA_limit_only_move_1(*a_branch->m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_branch->m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, option_tree_struct);
						//APLA::get_same_partition_SAPLA_210818(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_branch->m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_branch->m_rect.rectangle_entry_struct.list_partition_pointer, number_points);
						//APLA::get_same_partition_SAPLA_debug(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_branch->m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_branch->m_rect.rectangle_entry_struct.list_partition_pointer);
						/*--------------------210817 So far not need reconstructed time series---------*/
						//APLA::getAPLAReconstructSeries(*a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.reconstruct_time_series_vector_pointer);
						//APLA::getAPLAReconstructSeries(*a_branch->m_rect.rectangle_entry_struct.list_partition_pointer, *a_branch->m_rect.reconstruct_time_series_vector_pointer);
						/*-----------------------------------------------------------------------------*/
						//partition_rect_volume(a_node->m_branch[id_branch].m_rect, a_node->rectangle_in_node);
					}

					//get_convex_hull_by_node(*a_node);
				}
				case 2: {// convex tree

					/*......................................*/
#ifdef _DEBUG
					assert(option_tree_struct.type_representation == representation_option);
					assert(TMINNODES < TMAXNODES && TMINNODES >= 0);
					assert(!a_branch->m_child);
#endif
					/*......................................*/

					/*for conevex hull*/
					long double distance_max_node_0 = -1;
					long double distance_max_node_1 = -1;
					int id_min_0;
					int id_max_0;
					int id_min_1;
					int id_max_1;
					/**/

					long double distance_node_new_branch_max = -INF;
					long double distance_node_new_branch_current =  -INF;
					int id_branch_in_node_max = 0;
					int id_data_new_branch = a_branch->m_data;
					int id_node_data_with_branch = 0;

					vector<long double> distance_vector(a_node->m_count);
					for (int id_branch = 0; id_branch < a_node->m_count; id_branch++) {
						distance_vector[id_branch] = distance_node_new_branch_current = APLA::get_same_partition_SAPLA_limit_only_move_1(*a_branch->m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_branch->m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, option_tree_struct);
						if (distance_node_new_branch_current > distance_node_new_branch_max) {
							distance_node_new_branch_max = distance_node_new_branch_current;
							id_branch_in_node_max = id_branch;
							id_node_data_with_branch = a_node->m_branch[id_branch].m_data;
						}

						/*......................................*/
#ifdef _DEBUG
						if (option_tree_struct.type_representation != 3) {
							assert_endpoint_a_b(*a_branch->m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_branch->m_rect.rectangle_entry_struct.list_partition_pointer);
							assert_endpoint_a_b(*a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
						}
#endif
						/*......................................*/
						
						//APLA::get_same_partition_SAPLA_210818(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_branch->m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_branch->m_rect.rectangle_entry_struct.list_partition_pointer, number_points);
						//APLA::get_same_partition_SAPLA_debug(*a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer, *a_branch->m_rect.original_time_series_vector_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_branch->m_rect.rectangle_entry_struct.list_partition_pointer);
						/*--------------------210817 So far not need reconstructed time series---------*/
						//APLA::getAPLAReconstructSeries(*a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.reconstruct_time_series_vector_pointer);
						//APLA::getAPLAReconstructSeries(*a_branch->m_rect.rectangle_entry_struct.list_partition_pointer, *a_branch->m_rect.reconstruct_time_series_vector_pointer);
						/*-----------------------------------------------------------------------------*/
						//partition_rect_volume(a_node->m_branch[id_branch].m_rect, a_node->rectangle_in_node);
					}

					/*##########################    Pick Seeds (convex for sup node)  and choose partition   ###############################*/
					//get_convex_hull_by_node(*a_node);
					//after partition by new branch, compute max distance for seeds
					long double distance_convex_hull_updated = get_convex_hull_volume_by_distance_approximation_global(a_node->m_branch, a_node->m_count, a_node->rectangle_in_node);
					
					vector<int> m_partition_vector(a_node->m_count + 1, -1);
					vector<int> m_count(2, 1);

					/*|||||||||||||||||||||| sup new node ||||||||||||||||||||||||||||||||||*/
					Node* node_sup_new = new Node;
					if (!a_node->node_to_parent_branch_pointer) {
						/*..................................................................................................................................*/
#ifdef _DEBUG
						//assert(!(*a_newNode)->node_to_parent_branch_pointer);

#endif
						/*....................*/
						node_sup_new->m_count = 0;
						node_sup_new->m_level = 1;
						node_sup_new->m_branch = new Branch[TMAXNODES];

						/*++++++++++++++++++   partition 210618   +++++++++++++++++++++*/
						for (int id_branch = 0; id_branch < TMAXNODES; id_branch++) {
							node_sup_new->m_branch[id_branch].branch_to_parent_node_pointer = node_sup_new;
						}
						/*---------210709 rectangle in Node-------------*/
						allocate_rectangle(0, node_sup_new->rectangle_in_node);
						/*----------------------------------------------*/
						/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
					}
					/*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

					/*|||||||||||||||||||||||||    seeds in convex hull |||||||||||||||||||||||||||||||||||||||||*/
					if (distance_node_new_branch_max <= distance_convex_hull_updated) {

						/*!!!!!!!!!! convex hull for sup new node !!!!!!!!!!*/
						if (!a_node->node_to_parent_branch_pointer) {//Root
							/*..................................................................................................................................*/
#ifdef _DEBUG
							//assert(!(*a_newNode)->node_to_parent_branch_pointer);

#endif
							/*....................*/
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_convex_hull_updated;

							node_sup_new->rectangle_in_node.seed_0 = a_node->m_branch[a_node->rectangle_in_node.seed_0].m_data;
							node_sup_new->rectangle_in_node.seed_1 = a_node->m_branch[a_node->rectangle_in_node.seed_1].m_data;

							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);

							/*------------  211225  -------------*/
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer;
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer;
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer;
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer;
							/*-----------------------------------*/
						}
						/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


						int seed_0 = a_node->rectangle_in_node.seed_0;
						int seed_1 = a_node->rectangle_in_node.seed_1;

						m_partition_vector[seed_0] = 0;
						m_partition_vector[seed_1] = 1;

						// For new branch
						if (distance_vector[seed_0] < distance_vector[seed_1]) {
							m_partition_vector[a_node->m_count] = 0;
							++m_count[0];

							/*for conevex hull*/
							if (distance_max_node_0 < distance_vector[seed_0]) {
								distance_max_node_0 = distance_vector[seed_0];
								id_min_0 = seed_0;
								id_max_0 = TMAXNODES;
							}
							/**/
						}
						else {
							m_partition_vector[a_node->m_count] = 1;
							++m_count[1];

							/*for conevex hull*/
							if (distance_max_node_1 < distance_vector[seed_1]) { 
								distance_max_node_1 = distance_vector[seed_1];
								id_min_1 = seed_1;
								id_max_1 = TMAXNODES;
							}
							/**/
						}

						for (int id_branch = 0; id_branch < a_node->m_count; ++id_branch) {
							if (m_partition_vector[id_branch] == -1) {

								/*......................................*/
#ifdef _DEBUG
								assert(seed_0 != id_branch && seed_1 != id_branch);
#endif
								/*......................................*/

								long double distance_0 = APLA::get_sqrt_distance_sapla_same_endpoints(*a_node->m_branch[seed_0].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
								long double distance_1 = APLA::get_sqrt_distance_sapla_same_endpoints(*a_node->m_branch[seed_1].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);

								if (distance_0 <= distance_1) {

									if (m_count[0] + TMINNODES >= TMAXNODES + 1) {
										m_partition_vector[id_branch] = 1;
										++m_count[1];

										/*for conevex hull*/
										if (distance_max_node_1 < distance_1) {
											distance_max_node_1 = distance_1;
											id_min_1 = seed_1;
											id_max_1 = id_branch;
										}
										/**/
									}
									else {
										m_partition_vector[id_branch] = 0;
										++m_count[0];

										/*for conevex hull*/
										if (distance_max_node_0 < distance_0) {
											distance_max_node_0 = distance_0;
											id_min_0 = seed_0;
											id_max_0 = id_branch;
										} 
										/**/
									}
								}
								else if (distance_0 > distance_1) {

									if (m_count[1] + TMINNODES >= TMAXNODES + 1) {
										m_partition_vector[id_branch] = 0;
										++m_count[0];

										/*for conevex hull*/
										if (distance_max_node_0 < distance_0) {
											distance_max_node_0 = distance_0;
											id_min_0 = seed_0;
											id_max_0 = id_branch;
										}
										/**/
									}
									else {
										m_partition_vector[id_branch] = 1;
										++m_count[1];

										/*for conevex hull*/
										if (distance_max_node_1 < distance_1) {
											distance_max_node_1 = distance_1;
											id_min_1 = seed_1;
											id_max_1 = id_branch;
										}
										/**/
									}
								}
							}
						}
					}// seeds in new branch
					else {
						/*||||||||||||||||||  seeds in new branch ||||||||||||||||||||||||| */

						/*!!!!!!!!!! convex hull for sup new node !!!!!!!!!!*/
						if (!a_node->node_to_parent_branch_pointer) {
							/*..................................................................................................................................*/
#ifdef _DEBUG
							//assert(!(*a_newNode)->node_to_parent_branch_pointer);

#endif
							/*....................*/
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_node_new_branch_max;

							node_sup_new->rectangle_in_node.seed_0 = a_node->m_branch[id_branch_in_node_max].m_data;
							node_sup_new->rectangle_in_node.seed_1 = a_branch->m_data;//a_branch->m_rect.rectangle_entry_struct.original_time_series_vector_pointer.m_data;

							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*a_node->m_branch[id_branch_in_node_max].m_rect.rectangle_entry_struct.list_partition_pointer);
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*a_branch->m_rect.rectangle_entry_struct.list_partition_pointer);

							/*------------  211225  -------------*/
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = a_node->m_branch[id_branch_in_node_max].m_rect.rectangle_entry_struct.list_original_pointer;
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = a_branch->m_rect.rectangle_entry_struct.list_original_pointer;
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = a_node->m_branch[id_branch_in_node_max].m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
							node_sup_new->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = a_branch->m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
							/*-----------------------------------*/
						}
						/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

						int seed_0 = id_branch_in_node_max;
						int seed_1 = a_node->m_count;
						m_partition_vector[seed_0] = 0;
						m_partition_vector[seed_1] = 1;

						for (int id_branch = 0; id_branch < a_node->m_count; ++id_branch) {
							if (m_partition_vector[id_branch] == -1) {

								/*......................................*/
#ifdef _DEBUG
								assert(seed_0 != id_branch && seed_1 != id_branch);
#endif
								/*......................................*/

								long double distance_0 = APLA::get_sqrt_distance_sapla_same_endpoints(*a_node->m_branch[seed_0].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
								long double distance_1 = distance_vector[id_branch];// new branch and sub branch in node

								if (distance_0 <= distance_1) {

									if (m_count[0] + TMINNODES >= TMAXNODES + 1) {
										m_partition_vector[id_branch] = 1;
										++m_count[1];

										/*for conevex hull*/
										if (distance_max_node_1 < distance_1) {
											distance_max_node_1 = distance_1;
											id_min_1 = seed_1;
											id_max_1 = id_branch;
										} 
										/**/
									}
									else {
										m_partition_vector[id_branch] = 0;
										++m_count[0];

										/*for conevex hull*/
										if (distance_max_node_0 < distance_0) {
											distance_max_node_0 = distance_0;
											id_min_0 = seed_0;
											id_max_0 = id_branch;
										}
										/**/
									}
								}
								else if (distance_0 > distance_1) {

									if (m_count[1] + TMINNODES >= TMAXNODES + 1) {
										m_partition_vector[id_branch] = 0;
										++m_count[0];

										/*for conevex hull*/
										if (distance_max_node_0 < distance_0) {
											distance_max_node_0 = distance_0;
											id_min_0 = seed_0;
											id_max_0 = id_branch;
										}
										/**/
									}
									else {
										m_partition_vector[id_branch] = 1;
										++m_count[1];
										
										/*for conevex hull*/
										if (distance_max_node_1 < distance_1) {
											distance_max_node_1 = distance_1;
											id_min_1 = seed_1;
											id_max_1 = id_branch;
										}
										/**/
									}
								}
							}
						}

					}

					/*..................................................................................................................................*/
#ifdef _DEBUG
					assert((m_count[0] + m_count[1]) == a_node->m_count + 1);
					assert((m_count[0] >= TMINNODES) && (m_count[1] >= TMINNODES));
					for (auto&& au : m_partition_vector) {
						assert(au != -1 && (au == 0 || au == 1));
					}
					assert(a_node->m_level == 0);
#endif
					/*...................................................................................................................................*/
					/*#################################################################################################################*/
					

					/*#####################################        Load Nodes         #################################################*/

					vector<Branch*> address_branch_vector;
					vector<int> od_in_new_node_vector(TMAXNODES + 1, -1);
					/*================ new splited node==========================*/
					*a_newNode = new Node;
					(*a_newNode)->m_count = (*a_newNode)->m_level = 0;
					(*a_newNode)->m_branch = new Branch[TMAXNODES];
					for (int id_branch = 0; id_branch < TMAXNODES; id_branch++) {
						//(*a_newNode)->m_branch[id_branch].branch_to_parent_node_pointer = (*a_newNode);
						address_branch_vector.emplace_back(&(a_node->m_branch[id_branch]));
					}
					allocate_rectangle(0, (*a_newNode)->rectangle_in_node);
					a_node->m_count = 0;
					/*==========================================*/
					a_node->m_branch = new Branch[TMAXNODES];
					Node* targetNodes[] = { a_node, (*a_newNode) };
					for (int id_branch = 0; id_branch < TMAXNODES; ++id_branch) {
						od_in_new_node_vector[id_branch] = targetNodes[m_partition_vector[id_branch]]->m_count;

						targetNodes[m_partition_vector[id_branch]]->m_branch[targetNodes[m_partition_vector[id_branch]]->m_count] = *address_branch_vector[id_branch];// a_node->m_branch[id_branch];
						targetNodes[m_partition_vector[id_branch]]->m_branch[targetNodes[m_partition_vector[id_branch]]->m_count].branch_to_parent_node_pointer = targetNodes[m_partition_vector[id_branch]];
						++targetNodes[m_partition_vector[id_branch]]->m_count;

					}

					od_in_new_node_vector[TMAXNODES] = targetNodes[m_partition_vector[TMAXNODES]]->m_count;
					targetNodes[m_partition_vector[TMAXNODES]]->m_branch[targetNodes[m_partition_vector[TMAXNODES]]->m_count] = *a_branch;
					targetNodes[m_partition_vector[TMAXNODES]]->m_branch[targetNodes[m_partition_vector[TMAXNODES]]->m_count].branch_to_parent_node_pointer = targetNodes[m_partition_vector[TMAXNODES]];
					++targetNodes[m_partition_vector[TMAXNODES]]->m_count;
					/*#################################################################################################################*/
					address_branch_vector.emplace_back(a_branch);

					/*#####################################     get node cover  min max convex hull For new splited nodes      #################################################*/
					/*..................................................................................................................................*/
#ifdef _DEBUG
					//long double temp_distance_max_node_0 = get_convex_hull_volume_by_distance_approximation_global(a_node->m_branch, a_node->m_count, a_node->rectangle_in_node);
					//long double temp_distance_max_node_1 = get_convex_hull_volume_by_distance_approximation_global((*a_newNode)->m_branch, (*a_newNode)->m_count, (*a_newNode)->rectangle_in_node);

					//assert(distance_max_node_0 == temp_distance_max_node_0 && distance_max_node_1 == temp_distance_max_node_1);
					//assert(a_node->rectangle_in_node.seed_0 == od_in_new_node_vector[id_min_0] && a_node->rectangle_in_node.seed_1 == od_in_new_node_vector[id_max_0]);
					//assert((*a_newNode)->rectangle_in_node.seed_0 == od_in_new_node_vector[id_min_1] && (*a_newNode)->rectangle_in_node.seed_1 == od_in_new_node_vector[id_max_1]);
					//assert(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_0 == address_branch_vector[id_min_0]->m_data);
					//assert(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_1 == address_branch_vector[id_max_0]->m_data);
					//assert((*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_0 == address_branch_vector[id_min_1]->m_data);
					//assert((*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_1 == address_branch_vector[id_max_1]->m_data);
					/*assert_same_vector(*(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer), *(address_branch_vector[id_min_0]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer));
					assert_same_vector(*(a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer), *(address_branch_vector[id_max_0]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer));
					assert_same_vector(*(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer, *(address_branch_vector[id_min_1]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer));
					assert_same_vector(*(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer, *(address_branch_vector[id_max_1]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer));
					
					assert_same_SAPLA(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer, *address_branch_vector[id_min_0]->m_rect.rectangle_entry_struct.list_original_pointer);
					assert_same_SAPLA(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer, *address_branch_vector[id_max_0]->m_rect.rectangle_entry_struct.list_original_pointer);
					assert_same_SAPLA(*(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer, *address_branch_vector[id_min_1]->m_rect.rectangle_entry_struct.list_original_pointer);
					assert_same_SAPLA(*(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer, *address_branch_vector[id_max_1]->m_rect.rectangle_entry_struct.list_original_pointer);

					assert_same_SAPLA(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *address_branch_vector[id_min_0]->m_rect.rectangle_entry_struct.list_partition_pointer);
					assert_same_SAPLA(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, *address_branch_vector[id_max_0]->m_rect.rectangle_entry_struct.list_partition_pointer);
					assert_same_SAPLA(*(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *address_branch_vector[id_min_1]->m_rect.rectangle_entry_struct.list_partition_pointer);
					assert_same_SAPLA(*(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, *address_branch_vector[id_max_1]->m_rect.rectangle_entry_struct.list_partition_pointer);*/
#endif
				    /*...................................................................................................................................*/
					
					

					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_max_node_0;

					a_node->rectangle_in_node.seed_0 = od_in_new_node_vector[id_min_0];
					a_node->rectangle_in_node.seed_1 = od_in_new_node_vector[id_max_0];
					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_0 = address_branch_vector[id_min_0]->m_data;
					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_1 = address_branch_vector[id_max_0]->m_data;

					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*address_branch_vector[id_min_0]->m_rect.rectangle_entry_struct.list_partition_pointer);
					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*address_branch_vector[id_max_0]->m_rect.rectangle_entry_struct.list_partition_pointer);

					/*------------  211225  -------------*/
					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = address_branch_vector[id_min_0]->m_rect.rectangle_entry_struct.list_original_pointer;
					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = address_branch_vector[id_max_0]->m_rect.rectangle_entry_struct.list_original_pointer;
					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = address_branch_vector[id_min_0]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
					a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = address_branch_vector[id_max_0]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
					/*-----------------------------------*/

					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_max_node_1;

					(*a_newNode)->rectangle_in_node.seed_0 = od_in_new_node_vector[id_min_1];
					(*a_newNode)->rectangle_in_node.seed_1 = od_in_new_node_vector[id_max_1];
					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_0 = address_branch_vector[id_min_1]->m_data;
					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.seed_1 = address_branch_vector[id_max_1]->m_data;

					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*address_branch_vector[id_min_1]->m_rect.rectangle_entry_struct.list_partition_pointer);
					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*address_branch_vector[id_max_1]->m_rect.rectangle_entry_struct.list_partition_pointer);

					/*------------  211225  -------------*/
					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = address_branch_vector[id_min_1]->m_rect.rectangle_entry_struct.list_original_pointer;
					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = address_branch_vector[id_max_1]->m_rect.rectangle_entry_struct.list_original_pointer;
					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = address_branch_vector[id_min_1]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
					(*a_newNode)->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = address_branch_vector[id_max_1]->m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
					/*-----------------------------------*/

					long double temp_distance_max_node_0 = get_convex_hull_volume_by_distance_approximation_global(a_node->m_branch, a_node->m_count, a_node->rectangle_in_node);
					long double temp_distance_max_node_1 = get_convex_hull_volume_by_distance_approximation_global((*a_newNode)->m_branch, (*a_newNode)->m_count, (*a_newNode)->rectangle_in_node);

					/*......................................*/
#ifdef _DEBUG
					if (option_tree_struct.type_representation != 3) {
						assert_endpoint_a_b(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer, *a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
						assert_endpoint_a_b(*a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer, *a_node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
					}
#endif
					/*......................................*/
					
					/*#################################################################################################################*/


					/*#####################################        point to sup new node       #################################################*/
					//Branch branch_0;// Not poitner ???????????????????
					//branch_0.m_child = a_node;
					//branch_0.node_rect_pointer = make_shared<Rect>(a_node->rectangle_in_node);
					///*------------210814-----------*/
					//a_node->node_to_parent_branch_pointer = &branch_0;
					///*-----------------------------*/
					//branch_0.branch_to_parent_node_pointer = node_sup_new;

					//node_sup_new->m_branch[node_sup_new->m_count] = branch_0;
					//
					//++node_sup_new->m_count;

					//Branch branch_1;// Not poitner ???????????????????
					//branch_1.m_child = (*a_newNode);
					//branch_1.node_rect_pointer = make_shared<Rect>((*a_newNode)->rectangle_in_node);
					///*------------210814-----------*/
					//(*a_newNode)->node_to_parent_branch_pointer = &branch_1;
					///*-----------------------------*/
					//branch_1.branch_to_parent_node_pointer = node_sup_new;

					//node_sup_new->m_branch[node_sup_new->m_count] = branch_1;
					//
					//++node_sup_new->m_count;
					// 
					// 
					//Branch branch_0;// Not poitner ???????????????????

					if (!a_node->node_to_parent_branch_pointer) {
						/*..................................................................................................................................*/
#ifdef _DEBUG
						assert(!(*a_newNode)->node_to_parent_branch_pointer);

#endif
						/*....................*/
						node_sup_new->m_branch[node_sup_new->m_count].m_child = a_node;
						//node_sup_new->m_branch[node_sup_new->m_count].node_rect_pointer = make_shared<Rect>(a_node->rectangle_in_node);
						/*------------210814-----------*/
						a_node->node_to_parent_branch_pointer = &node_sup_new->m_branch[node_sup_new->m_count];
						/*-----------------------------*/
						node_sup_new->m_branch[node_sup_new->m_count].branch_to_parent_node_pointer = node_sup_new;

						//node_sup_new->m_branch[node_sup_new->m_count] = branch_0;

						++node_sup_new->m_count;

						//Branch branch_1;// Not poitner ???????????????????
						node_sup_new->m_branch[node_sup_new->m_count].m_child = (*a_newNode);
						//node_sup_new->m_branch[node_sup_new->m_count].node_rect_pointer = make_shared<Rect>((*a_newNode)->rectangle_in_node);
						/*------------210814-----------*/
						(*a_newNode)->node_to_parent_branch_pointer = &node_sup_new->m_branch[node_sup_new->m_count];
						/*-----------------------------*/
						node_sup_new->m_branch[node_sup_new->m_count].branch_to_parent_node_pointer = node_sup_new;

						//node_sup_new->m_branch[node_sup_new->m_count] = branch_1;

						++node_sup_new->m_count;

					}
					else {
						delete node_sup_new;
						node_sup_new = nullptr;
					}
					/*#################################################################################################################*/

					/*..................................................................................................................................*/
#ifdef _DEBUG
					assert((a_node->m_count + (*a_newNode)->m_count) == TMAXNODES + 1);
					//vector<int> id_data_vector_0;

					//vector<int> id_data_vector_1;
					//cout << "Original Address: " << endl;
					//for (int id_branch = 0; id_branch < address_branch_vector.size(); id_branch++) {
					//	cout<< address_branch_vector[id_branch] << endl;
					//}
					//cout << "\n Node Address: " << endl;
					//for (int id_branch = 0; id_branch < a_node->m_count; id_branch++) {
					//	cout << &(a_node->m_branch[id_branch]) << endl;
					//	id_data_vector_0.emplace_back(a_node->m_branch[id_branch].m_data);
					//}
					//cout << "\n New Node Address: " << endl;
					//for (int id_branch = 0; id_branch < (*a_newNode)->m_count; id_branch++) {
					//	cout << &((*a_newNode)->m_branch[id_branch]) << endl;
					//	id_data_vector_1.emplace_back((*a_newNode)->m_branch[id_branch].m_data);
					//}

					////not same;

					//for (auto&& au: id_data_vector_0) {
					//	if (find(id_data_vector_1.begin(), id_data_vector_1.end(), au) != id_data_vector_1.end()) {
					//		assert(0);
					//	}
					//}
#endif
					/*...................................................................................................................................*/
					return true;//1 already split a_node into a_newNode.
					break;
				/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
				}
				default:
					assert(0);
					break;
				}
			}
			else {//Internal node
				switch (option_tree_struct.type_tree) {
				case 1: {
					break;
				}
				case 2: {

					vector<Branch*> address_branch_vector;

					for (int id_branch = 0; id_branch < a_node->m_count; id_branch++) {
						address_branch_vector.emplace_back(&(a_node->m_branch[id_branch]));
					}

					address_branch_vector.emplace_back(a_branch);

					/*for conevex hull*/
					long double distance_current = -1;
					long double distance_max = -1;
					int seed_0=0;
					int seed_1=1;
					/**/

					for (int id_branch_a = 0; id_branch_a < a_node->m_count; ++id_branch_a) {
						for (int id_branch_b = id_branch_a + 1; id_branch_b < a_node->m_count + 1; ++id_branch_b) {
							distance_current = combine_compute_partition_node_by_convex(address_branch_vector[id_branch_a]->m_child->rectangle_in_node, address_branch_vector[id_branch_b]->m_child->rectangle_in_node);
							if (distance_current > distance_max) {
								distance_max = distance_current;
								seed_0 = id_branch_a;
								seed_1 = id_branch_b;
							}
						
						}
					}

					vector<int> m_partition_vector(a_node->m_count + 1, -1);
					vector<int> m_count(2, 1);
					m_partition_vector[seed_0] = 0;
					m_partition_vector[seed_1] = 1;

					for (int id_branch = 0; id_branch < a_node->m_count + 1; ++id_branch) {
						if (m_partition_vector[id_branch] == -1) {
							long double distance_0 = combine_compute_partition_node_by_convex(address_branch_vector[seed_0]->m_child->rectangle_in_node, address_branch_vector[id_branch]->m_child->rectangle_in_node);
							long double distance_1 = combine_compute_partition_node_by_convex(address_branch_vector[seed_1]->m_child->rectangle_in_node, address_branch_vector[id_branch]->m_child->rectangle_in_node);
							if (distance_0 <= distance_1) {

								if (m_count[0] + TMINNODES >= TMAXNODES + 1) {
									m_partition_vector[id_branch] = 1;
									++m_count[1];

									
								}
								else {
									m_partition_vector[id_branch] = 0;
									++m_count[0];

									
								}
							}
							else if (distance_0 > distance_1) {

								if (m_count[1] + TMINNODES >= TMAXNODES + 1) {
									m_partition_vector[id_branch] = 0;
									++m_count[0];

									
								}
								else {
									m_partition_vector[id_branch] = 1;
									++m_count[1];
									
								}
							}
						
						}
					}

					/*..................................................................................................................................*/
#ifdef _DEBUG
					assert((m_count[0] + m_count[1]) == a_node->m_count + 1);
					assert((m_count[0] >= TMINNODES) && (m_count[1] >= TMINNODES));
					for (auto&& au : m_partition_vector) {
						assert(au != -1 && (au == 0 || au == 1));
					}
					assert(a_node->m_level > 0);
#endif
					/*...................................................................................................................................*/


					/*#####################################        Load Nodes         #################################################*/
					
					vector<int> od_in_new_node_vector(TMAXNODES + 1, -1);
					/*================ new splited node==========================*/
					*a_newNode = new Node;
					(*a_newNode)->m_count = 0;
					(*a_newNode)->m_branch = new Branch[TMAXNODES];
					(*a_newNode)->m_level = a_node->m_level;
					allocate_rectangle(0, (*a_newNode)->rectangle_in_node);
					a_node->m_count = 0;
					/*==========================================*/
					a_node->m_branch = new Branch[TMAXNODES];
					Node* targetNodes[] = { a_node, (*a_newNode) };
					for (int id_branch = 0; id_branch < TMAXNODES+1; ++id_branch) {
						od_in_new_node_vector[id_branch] = targetNodes[m_partition_vector[id_branch]]->m_count;

						targetNodes[m_partition_vector[id_branch]]->m_branch[targetNodes[m_partition_vector[id_branch]]->m_count] = *address_branch_vector[id_branch];// a_node->m_branch[id_branch];
						targetNodes[m_partition_vector[id_branch]]->m_branch[targetNodes[m_partition_vector[id_branch]]->m_count].branch_to_parent_node_pointer = targetNodes[m_partition_vector[id_branch]];
						++targetNodes[m_partition_vector[id_branch]]->m_count;

					}

					/*#################################################################################################################*/

					/*..................................................................................................................................*/
#ifdef _DEBUG
					assert((a_node->m_count + (*a_newNode)->m_count) == TMAXNODES + 1);
#endif
					/*...................................................................................................................................*/

					return true;
					break;
				}
				default:
					assert(0);
					break;
				}
			}
			/*....................*/
#ifdef _DEBUG
			assert_node(*a_node);
#endif
			/*....................*/
			break;
		}
		default:
			assert(0);
			break;
		}
		/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

		SplitNode(a_node, a_branch, a_newNode);

		return true;//1 already split a_node into a_newNode.
	}
}

// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::DisconnectBranch(Node* a_node, int a_index)
{
	ASSERT(a_node && (a_index >= 0) && (a_index < TMAXNODES));
	ASSERT(a_node->m_count > 0);

	// Remove element by swapping with the last element to prevent gaps in array
	a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];

	--a_node->m_count;
}


// Pick a branch.  Pick the one that will need the smallest increase in area to accommodate the new rectangle.
// This will result in the least total area for the covering rectangles in the current node.
// Used in InsertRectRec() for inernal node
// In case of a tie, pick the one which was smaller before, to get the best resolution when searching.
RTREE_TEMPLATE
int RTREE_PARTITION_QUAL::PickBranch(Rect* a_rect, Node* a_node)
{
	
	/*......................................*/
	
#ifdef _DEBUG
	ASSERT(a_rect && a_node);
	is_rectangle_entry(*a_rect);
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6://ICDE07
	case 8: {// Initial 200706
		assert_rectangle(*a_rect);
		assert_node(*a_node);
		break;
	}
	default:
		assert(0);
		break;
	}
#endif
	/*......................................*/
	int best = 0;
	/*++++++++++++++++++++++  partition 210706  +++++++++++++++++++*/
		// before split node, parition new branches and all sub branches in old node. Thus they all have same partitions
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706

		best = pick_branch(*a_rect, *a_node);
		return best;
		break;
	}
	default:{
		bool firstTime = true;
		ELEMTYPEREAL increase;
		long double bestIncr = -INF;
		long double area;
		ELEMTYPEREAL bestArea;
		

		Rect tempRect;
		allocate_minmax_vector_in_rectangle(tempRect);
		//tempRect.m_ min = new ELEMTYPE[NUMDIMS];
		//tempRect.m_max = new ELEMTYPE[NUMDIMS];

		for (int index = 0; index < a_node->m_count; ++index)
		{
			Rect* curRect = &a_node->m_branch[index].m_rect;// current branch convering

			area = CalcRectVolume(curRect);// origianl branch volume

			//tempRect = CombineRect(a_rect, curRect);
			CombineRect(a_rect, curRect, &tempRect);// get volume of combined covering set.

			increase = CalcRectVolume(&tempRect) - area;

			if ((increase < bestIncr) || firstTime)
			{
				best = index;
				bestArea = area;
				bestIncr = increase;
				firstTime = false;
			}
			else if ((increase == bestIncr) && (area < bestArea))
			{
				best = index;
				bestArea = area;
				bestIncr = increase;
			}
		}

		free_minmax_vector_in_rectangle(tempRect);
		//delete[] tempRect.m_ min;
		//tempRect.m_ min = nullptr;
		//delete[] tempRect.m_max;
		//tempRect.m_max = nullptr;
		break;
	}
	}
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	return best;
}

//************************************
// Method:pick_branch
// Qualifier: 1  a_rect is new entry, has not benn paritioned.
//            2 
// date:210810
// author:
//************************************
RTREE_TEMPLATE
template<typename T, typename Y>//210810 pick branch by distance
int RTREE_PARTITION_QUAL::pick_branch(T& const a_rect, Y& const a_node) {
	/*......................*/
#ifdef _DEBUG
	is_rectangle_entry(a_rect);
	assert_rectangle(a_rect);
	assert_node(a_node);
#endif
	/*......................*/
	const int option_compare = 0;
	bool first_time = true;
	int id_branch_best = 0;
	long double difference_volume_current = 0;
	long double difference_volume_best = INF;
	long double volume_current = INF;
	long double volume_combine = INF;
	long double volume_best = INF;

	/*T rectangle_combine_temp;
	allocate_rectangle(0, rectangle_combine_temp);*/

	for (int id_branch = 0; id_branch < a_node.m_count; ++id_branch) {

		switch (option_tree_struct.type_tree) {
		case 1: {
			/*........................................................................................*/
#ifdef _DEBUG
			assert_rectangle_same(a_node.m_branch[id_branch].m_rect, a_node.m_branch[id_branch].m_child->rectangle_in_node);
#endif
			/*........................................................................................*/
			//update partition by a_rect, get new partiiton, volume of rectangle
			//volume_combine = combine_compute_partition_rectangle_by_entry(a_rect, a_node.m_branch[id_branch].m_rect, rectangle_combine_temp);
			volume_combine = combine_compute_partition_rectangle_by_entry_210905(a_rect, a_node.m_branch[id_branch].m_rect);
			a_node.m_branch[id_branch].m_child->rectangle_in_node = a_node.m_branch[id_branch].m_rect;
			break;
		}
		case 2: {
			//update partition by a_rect, get new partiiton, volume of rectangle
			//volume_combine = combine_compute_partition_rectangle_by_entry(a_rect, a_node.m_branch[id_branch].m_rect, rectangle_combine_temp);
			volume_combine = combine_compute_partition_rectangle_by_convex(a_rect, a_node.m_branch[id_branch].m_child->rectangle_in_node);
			
			break;
		}
		default:
			assert(0);
		}
		/*##################################################################################*/
		switch (option_compare) {
		case 0: {
			  if (volume_combine < volume_best) {
				  id_branch_best = id_branch;
				  volume_best = volume_combine;
			  }
			  break;
		}
		case 1:{// origianl Rtree,
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			volume_current = a_node.m_branch[id_branch].m_rect.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;
			difference_volume_current = volume_combine - volume_current;

			if ((difference_volume_current < difference_volume_best) || first_time) {
				id_branch_best = id_branch;
				volume_best = volume_current;
				difference_volume_best = difference_volume_current;
				first_time = false;
			}
			else if ((difference_volume_current == difference_volume_best) && (volume_current < volume_best)) {
				id_branch_best = id_branch;
				volume_best = volume_current;
				difference_volume_best = difference_volume_current;
			}
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			break;
		}
		default:
			assert(0);
		}
		//clear_vector_in_rectangle(rectangle_combine_temp);// clear vector
	}

	return id_branch_best;
}

//Combine two rectangles into larger one containing both
RTREE_TEMPLATE
typename RTREE_PARTITION_QUAL::Rect& RTREE_PARTITION_QUAL::CombineRect(Rect* a_rectA, Rect* a_rectB)
{
	ASSERT(a_rectA && a_rectB);
	ASSERT(!a_rectA->m_min.empty()  && !a_rectB->m_min.empty());

	Rect newRect;
	//newRect.m _min = new ELEMTYPE[NUMDIMS];
	//newRect.m_max = new ELEMTYPE[NUMDIMS];
	allocate_minmax_vector_in_rectangle(newRect);

	for (int index = 0; index < NUMDIMS; ++index) {

		newRect.m_min[index] = (Min)(a_rectA->m_min[index], a_rectB->m_min[index]);
		newRect.m_max[index] = (Max)(a_rectA->m_max[index], a_rectB->m_max[index]);
	}

	return newRect;
}

//***************************************************************
// Method:CombineRect
// Qualifier: combine two branches
// Input:
// Output: 
// Note:
// date: 210618
// author:
//***************************************************************
RTREE_TEMPLATE
typename RTREE_PARTITION_QUAL::Rect RTREE_PARTITION_QUAL::CombineRect(Rect* a_rectA, Rect* a_rectB, Rect* a_rect) {
	//ASSERT(a_rectA && a_rectB);
	//ASSERT(a_rectA->m_min && a_rectB->m_min && a_rect->m_min);

	for (int index = 0; index < NUMDIMS; ++index) {

		a_rect->m_min[index] = (Min)(a_rectA->m_min[index], a_rectB->m_min[index]);
		a_rect->m_max[index] = (Max)(a_rectA->m_max[index], a_rectB->m_max[index]);
	}

	return *a_rect;
}

//***************************************************************
// Method:combine_compute_partition_rectangle_by_entry
// Qualifier: combine rectangle_entry in rectangle_b by vector
// Input:
// Output: 
// Note: The purpuse of function is get volume of combined rectangle for pick branch
// date: 210728
// author:
//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
long double RTREE_PARTITION_QUAL::combine_compute_partition_rectangle_by_entry(T& const rectangle_entry, Y& const rectangle_node, U& const rectangle_combine) {
	/*......................................*/
#ifdef _DEBUG
	//assert(rectangle_entry->is_entry());
	//assert(rectangle_node.is_node());
	is_rectangle_entry(rectangle_entry);
	assert_rectangle(rectangle_entry);
	assert_rectangle(rectangle_node);
	assert_rectangle(rectangle_combine);
#endif
	/*......................................*/

	copy_rect_original_series_approximation_vector_pointer(rectangle_node, rectangle_combine);
	copy_rect_original_series_approximation_vector_pointer(rectangle_entry, rectangle_combine);

	/*........................................................................................*/
#ifdef _DEBUG
	assert(rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.size_entry() == 1 + rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry());
#endif
	/*........................................................................................*/

	/**/
	option_tree_struct.type_representation = representation_option;
	
	/**/

	for (int id_entry = 0; id_entry < rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry(); id_entry++) {
		APLA::get_same_partition_SAPLA_limit_only_move_1(*rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_entry], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_entry], option_tree_struct);
		//APLA::get_same_partition_SAPLA_210818(*rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_entry], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_entry], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()], number_points);
		//APLA::get_same_partition_SAPLA_debug(*rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_entry], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_entry], *rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()]);
		/*--------------------210817 So far not need reconstructed time series---------*/
		//APLA::getAPLAReconstructSeries(*rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_entry], *rectangle_combine.reconstruct_time_series_pointer_vector[id_entry]);
		//APLA::getAPLAReconstructSeries(*rectangle_combine.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()], *rectangle_combine.reconstruct_time_series_pointer_vector[rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry()]);
	    /*----------------------------------------------------------------------------*/
	}
	get_convex_hull_by_rectangle(rectangle_node);
	//get_convex_hull_by_rectangle(rectangle_combine);
	/*......................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_entry);
	assert_rectangle(rectangle_node);
	//assert_rectangle(rectangle_combine);
	/*if (rectangle_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size() > 1) {
		long double distance_1 = get_sqrt_distance_sapla_same_endpoints(*rectangle_entry.rectangle_entry_struct.list_partition_pointer, *rectangle_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
		long double distance_2 = get_sqrt_distance_sapla_same_endpoints(*rectangle_entry.rectangle_entry_struct.list_partition_pointer, *rectangle_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
		long double distance_max_0 = max(distance_1, distance_2);
		long double distance_max_1 = max(rectangle_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt, distance_max_0);
		assert(rectangle_combine.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt == distance_max_1);
	}*/
#endif
	/*......................................*/

	return get_convex_hull_by_rectangle(rectangle_combine);
}

//The purpuse of function is get volume of combined rectangle for pick branch
RTREE_TEMPLATE
template<typename T, typename Y>//210905 return difference of volume. rectangle_entry is entry. No need temp combine rectangle, speed up
long double RTREE_PARTITION_QUAL::combine_compute_partition_rectangle_by_entry_210905(T& const rectangle_entry, Y& const rectangle_node) {
	/*......................................*/
#ifdef _DEBUG
	//assert(rectangle_entry->is_entry());
	//assert(rectangle_node.is_node());
	is_rectangle_entry(rectangle_entry);
	assert_rectangle(rectangle_entry);
	assert_rectangle(rectangle_node);
#endif
	/*......................................*/

	option_tree_struct.type_representation = representation_option;

	long double distance_sqrt_SAPLA_with_new_entry_max = -INF;
	long double distance_sqrt_SAPLA_with_new_entry_current = -INF;
	int number_points = INF;
	const int size_entry = rectangle_node.rectangle_inter_node_struct.vector_vector_struct.size_entry();

	for (int id_entry = 0; id_entry < size_entry; id_entry++) {
		distance_sqrt_SAPLA_with_new_entry_current = APLA::get_same_partition_SAPLA_limit_only_move_1(*rectangle_entry.rectangle_entry_struct.original_time_series_vector_pointer, *rectangle_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_entry], *rectangle_entry.rectangle_entry_struct.list_partition_pointer, *rectangle_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_entry], option_tree_struct);
		if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_SAPLA_with_new_entry_current) {
			distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_SAPLA_with_new_entry_current;
		}
	}
	distance_sqrt_SAPLA_with_new_entry_current = get_convex_hull_by_rectangle(rectangle_node);

	if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_SAPLA_with_new_entry_current) {
		distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_SAPLA_with_new_entry_current;
	}

	/*......................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_entry);
	assert_rectangle(rectangle_node);

	T rectangle_combine_temp;
	allocate_rectangle(0, rectangle_combine_temp);
	copy_rect_original_series_approximation_vector_pointer(rectangle_node, rectangle_combine_temp);
	copy_rect_original_series_approximation_vector_pointer(rectangle_entry, rectangle_combine_temp);
	const long double distance_sqrt_SAPLA_with_new_entry_max_test = get_convex_hull_by_rectangle(rectangle_combine_temp);
	assert(float(distance_sqrt_SAPLA_with_new_entry_max) == float(distance_sqrt_SAPLA_with_new_entry_max_test));
	clear_vector_in_rectangle(rectangle_combine_temp);// clear vector
#endif
	/*......................................*/

	return distance_sqrt_SAPLA_with_new_entry_max;
}

//211225 For the convex hull. The purpuse of function is get volume of combined rectangle for pick branch
RTREE_TEMPLATE
template<typename T, typename Y>// return difference of volume. rectangle_entry is entry. No need temp combine rectangle, speed up
long double RTREE_PARTITION_QUAL::combine_compute_partition_rectangle_by_convex(T& const rectangle_entry, Y& const rectangle_node) {
	/*......................................*/
#ifdef _DEBUG
	//assert(rectangle_entry->is_entry());
	//assert(rectangle_node.is_node());
	is_rectangle_entry(rectangle_entry);
	assert_rectangle(rectangle_entry);
	assert_rectangle(rectangle_node);
	assert(option_tree_struct.type_representation == representation_option);
#endif
	/*......................................*/

	long double distance_sqrt_SAPLA_with_new_entry_max = -INF;
	long double distance_sqrt_SAPLA_with_new_entry_current = -INF;


	long double distance_sqrt_with_min = distance_sqrt_SAPLA_with_new_entry_max = APLA::get_same_partition_SAPLA_limit_only_move_1(*rectangle_entry.rectangle_entry_struct.original_time_series_vector_pointer, *rectangle_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer, *rectangle_entry.rectangle_entry_struct.list_partition_pointer, *rectangle_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, option_tree_struct);
	long double distance_sqrt_with_max = APLA::get_same_partition_SAPLA_limit_only_move_1(*rectangle_entry.rectangle_entry_struct.original_time_series_vector_pointer, *rectangle_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer, *rectangle_entry.rectangle_entry_struct.list_partition_pointer, *rectangle_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, option_tree_struct);

	if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_with_max) {
		distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_with_max;
	}

	rectangle_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = APLA::get_sqrt_distance_sapla_same_endpoints(*rectangle_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);

	if (distance_sqrt_SAPLA_with_new_entry_max < rectangle_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt) {
		distance_sqrt_SAPLA_with_new_entry_max = rectangle_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;
	}

	/*......................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_entry);
	assert_rectangle(rectangle_node);
#endif
	/*......................................*/

	return distance_sqrt_SAPLA_with_new_entry_max;
}


//211227 The purpuse of function is get volume of internal nodt for pick branch
RTREE_TEMPLATE
template<typename T, typename Y>// return difference of volume. rectangle_node_0 is node rectanle. No need temp combine rectangle, speed up
long double RTREE_PARTITION_QUAL::combine_compute_partition_node_by_convex(T& const rectangle_node_0, Y& const rectangle_node_1) {
	/*......................................*/
#ifdef _DEBUG
	//assert(rectangle_node_0->is_entry());
	//assert(rectangle_node.is_node());
	assert_rectangle(rectangle_node_0);
	assert_rectangle(rectangle_node_1);
	assert(option_tree_struct.type_representation == representation_option);

	APLA::assert_has_same_endpoints(*rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	APLA::assert_has_same_endpoints(*rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
#endif
	/*......................................*/

	long double distance_sqrt_SAPLA_with_new_entry_max = -INF;
	long double distance_sqrt_SAPLA_with_new_entry_current = -INF;


	long double distance_sqrt_min0_min1 = distance_sqrt_SAPLA_with_new_entry_max = APLA::get_same_partition_SAPLA_limit_only_move_1(*rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer, *rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer, *rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, option_tree_struct);
	long double distance_sqrt_max0_max1 = APLA::get_same_partition_SAPLA_limit_only_move_1(*rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer, *rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer, *rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, *rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, option_tree_struct);
	if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_max0_max1) {
		distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_max0_max1;
	}
	
	long double distance_sqrt_min0_max0 = rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = APLA::get_sqrt_distance_sapla_same_endpoints(*rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min0_max0) {
		distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min0_max0;
	}
	
	long double distance_sqrt_min1_max1 = rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = APLA::get_sqrt_distance_sapla_same_endpoints(*rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min1_max1) {
		distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min1_max1;
	}
	
	long double distance_sqrt_min0_max1 = APLA::get_sqrt_distance_sapla_same_endpoints(*rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min0_max1) {
		distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min0_max1;
	}
	
	long double distance_sqrt_min1_max0 = APLA::get_sqrt_distance_sapla_same_endpoints(*rectangle_node_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_node_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);

	if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min1_max0) {
		distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min1_max0;
	}

	

	/*......................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_node_0);
	assert_rectangle(rectangle_node_1);
#endif
	/*......................................*/

	return distance_sqrt_SAPLA_with_new_entry_max;
}

//RTREE_TEMPLATE
//typename RTREE_PARTITION_QUAL::Rect& RTREE_PARTITION_QUAL::CombineRect(Rect* a_rectA, Rect* a_rectB, Rect& a_rect) {
//	ASSERT(a_rectA && a_rectB);
//
//	for (int index = 0; index < NUMDIMS; ++index) {
//
//		a_rect.m_ min[index] = (Min)(a_rectA->m_ min[index], a_rectB->m_ min[index]);
//		a_rect.m_max[index] = (Max)(a_rectA->m_max[index], a_rectB->m_max[index]);
//	}
//
//	return a_rect;
//}


RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::assertRect(Rect& a_rect) {
	for (int index = 0; index < NUMDIMS; ++index) {
		ASSERT(a_rect.m_min[index] < a_rect.m_max[index]);
	}
}

// Split a node.(internal or leaf ?)
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
// Used in AddBranch()
// a_branch is new coming branch
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::SplitNode(Node* a_node, const Branch* a_branch, Node** a_newNode)
{
	ASSERT(a_node);
	ASSERT(a_branch);

	/*..............................................................................*/
#ifdef _DEBUG
	//update_branch_pointer_vector(*a_node);
	//update_branch_pointer_vector(**a_newNode);
#endif
	/*..............................................................................*/

	// Could just use local here, but member or external is faster since it is reused
	PartitionVars localVars;

	/*++++++++++++    Initial PartitionVars localVars    +++++++++++++++++++++++++++++++++++*/

	/*~~~~~~~~~~~~~~~210706 vector instead pointer ~~~~~~~~~~~~~~~~~~~~~*/
	//localVars.m_partition = new int[TMAXNODES + 1];
	localVars.m_partition.assign(TMAXNODES + 1, -1);
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	localVars.m_branchBuf = new Branch[TMAXNODES + 1];// is the branch structure for the new branch
	allocate_minmax_vector_in_rectangle(localVars.m_cover[0]);// vector min max 
	allocate_minmax_vector_in_rectangle(localVars.m_cover[1]);
	/*~~~~~~~~ 210726 allocate rectangle: new pointer  ~~~~~~~~~~~~~~~~*/
	allocate_rectangle(0, localVars.m_coverSplit);// new rectanle
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//localVars.m_coverSplit.m_ min = new ELEMTYPE[NUMDIMS];
	//localVars.m_coverSplit.m_max = new ELEMTYPE[NUMDIMS];
	//localVars.m_cover[0].m_min = new ELEMTYPE[NUMDIMS];
	//localVars.m_cover[0].m_max = new ELEMTYPE[NUMDIMS];
	//localVars.m_cover[1].m_ min = new ELEMTYPE[NUMDIMS];
	//localVars.m_cover[1].m_max = new ELEMTYPE[NUMDIMS];
	//memory_account[7] += sizeof(localVars.m_partition)*(TMAXNODES + 1);
	//memory_account[6] += sizeof(localVars.m_branchBuf)*(TMAXNODES + 1);

	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	PartitionVars* parVars = &localVars;

	// Load all the branches into a buffer, initialize old node
	// Leaf or Internal Node?
	GetBranches(a_node, a_branch, parVars);

	// The core is to get a_parVars->m_partition[index]; 0/1 group id of each branch.
	ChoosePartition(parVars, TMINNODES);

	// Create a new node to hold (about) half of the branches
	*a_newNode = AllocNode();
	//pointer[3] = *a_newNode;
	(*a_newNode)->m_level = a_node->m_level;

	// Put branches from buffer into 2 nodes according to the chosen partition
	a_node->m_count = 0;
	
	/*==================  partition 210618 =========================*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		clear_vector_in_rectangle(a_node->rectangle_in_node);
		break;
	}
	default: {
		assert(0);
		break;
	}
	}
	/*==================================================*/

	/*.....      210729    ...........*/
#ifdef _DEBUG
	assert((*a_newNode)->m_count == 0);
#endif
	/*................................*/

	//Input: a_parVars->m_partition[index]; 0/1 group id of each branch.
	LoadNodes(a_node, *a_newNode, parVars);

	ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);

	/*--------------   Clear Vector   ---------------------*/
	clear_rectangle_PartitionVars(localVars.m_coverSplit);
	parVars = nullptr;
	/*-----------------------------------------------------*/

	/*.....      210729    ...........*/
#ifdef _DEBUG
	//update_branch_pointer_vector(*a_node);
	//update_branch_pointer_vector(**a_newNode);
#endif
	/*................................*/
}


// Calculate the n-dimensional volume of a rectangle      : can get the area of rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::RectVolume(Rect* a_rect)
{
	ASSERT(a_rect);

	ELEMTYPEREAL volume = (ELEMTYPEREAL)1;

	for (int index = 0; index < NUMDIMS; ++index)
	{
		volume *= a_rect->m_max[index] - a_rect->m_min[index];
		//assert(volume != INF);
		//cout << "volume: " << volume << endl;
	}

	if (volume == INF || !isfinite(volume)) {
		volume = 0;
	}

	assert(volume != INF);
	ASSERT(!isinf(volume));
	//ASSERT(volume >= DBL_TRUE_MIN);
	ASSERT(isfinite(volume));
	//ASSERT(isnormal(volume));
	//ASSERT(volume > (ELEMTYPEREAL)0);

	return volume;
}

// Calculate the n-dimensional APCA Point volume of a rectangle      : can get the area of rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::APCARectVolume(Rect* a_rect)
{
	ASSERT(a_rect);

	ELEMTYPEREAL f_temp_difference = NULL;
	/*--------------------------------------Old version-----------------------------------------*/
	ELEMTYPEREAL volume = (ELEMTYPEREAL)a_rect->m_max[0] + 1;
	/*-------------------------------------------------------------------------------------------*/
	/*--------------------------------------New version-----------------------------------------*/
	//ELEMTYPEREAL volume = 0;
	/*-------------------------------------------------------------------------------------------*/
	//cout << "max[" << 0 << "] = " << a_rect->m_max[0]+1<< " - min[" << -1 << "] = " << 0 << endl;
	/*--------------------------------------Old version-----------------------------------------*/
	for (int index = 1; index < NUMDIMS; index++)
		/*-------------------------------------------------------------------------------------------*/
		//for (int index = 0; index < NUMDIMS; index++)
	{

		ASSERT(a_rect->m_max[index] >= a_rect->m_min[index]);

		/*--------------------------------------Old version-----------------------------------------*/
		if (index & 1) {//odd max value - min value
			//test
			//if (a_rect->m_max[index] < a_rect->m_ min[index]) {
				//cout << "max: " << a_rect->m_max[index] << ", min: " << a_rect->m_ min[index] << endl;
				//ASSERT(0);
			//}
			f_temp_difference = a_rect->m_max[index] - a_rect->m_min[index];
			//#ifdef _DEBUG
			//			cout<<" temp difference: "<< f_temp_difference <<endl;
			//#endif
						//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_ min[index] << endl;


						//if (f_temp_difference < DBL_TRUE_MIN) continue;//191128 will influence result
			volume *= f_temp_difference;

		}
		else {//even
			//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_ min[index] << endl;
			//ASSERT(a_rect->m_max[index] == a_rect->m_ min[index]);
			volume *= (a_rect->m_max[index] - a_rect->m_min[index - 2]);

		}
		/*-------------------------------------------------------------------------------------------*/
		/*--------------------------------------New version-----------------------------------------*/
		//volume += (a_rect->m_max[index] - a_rect->m_ min[index]);
		/*-------------------------------------------------------------------------------------------*/
		//cout << volume << endl;
	}


#ifdef _DEBUG
	//if (volume < DBL_TRUE_MIN) {
	//	for (int index = 1; index < NUMDIMS; index++) {
	//		if (index & 1) {//odd
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_ min[index] << endl;
	//		}
	//		else {//even
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index - 2 << "] = " << a_rect->m_ min[index-2] << endl;
	//		}
	//	}
	//}
	ASSERT(!isinf(volume));
	ASSERT(volume >= DBL_TRUE_MIN);
	ASSERT(isfinite(volume));
	ASSERT(isnormal(volume));
	ASSERT(volume > (ELEMTYPEREAL)0);
#endif
	//cout << volume << endl;
	return fabs(volume);
}

//191128 accumulation, min - min - 2
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::apca_rect_accumulation_volume(Rect* a_rect) {
	ASSERT(a_rect);

	ELEMTYPEREAL f_temp_difference = NULL;
	/*--------------------------------------Old version-----------------------------------------*/
	ELEMTYPEREAL volume = (ELEMTYPEREAL)a_rect->m_max[0] + 1;
	/*-------------------------------------------------------------------------------------------*/
	/*--------------------------------------New version-----------------------------------------*/
	//ELEMTYPEREAL volume = 0;
	/*-------------------------------------------------------------------------------------------*/
	//cout << "max[" << 0 << "] = " << a_rect->m_max[0]+1<< " - min[" << -1 << "] = " << 0 << endl;
	/*--------------------------------------Old version-----------------------------------------*/
	for (int index = 1; index < NUMDIMS; index++)
		/*-------------------------------------------------------------------------------------------*/
		//for (int index = 0; index < NUMDIMS; index++)
	{

		ASSERT(a_rect->m_max[index] >= a_rect->m_min[index]);

		/*--------------------------------------Old version-----------------------------------------*/
		if (index & 1) {//odd max value - min value
			//test
			//if (a_rect->m_max[index] < a_rect->m_ min[index]) {
				//cout << "max: " << a_rect->m_max[index] << ", min: " << a_rect->m_ min[index] << endl;
				//ASSERT(0);
			//}
			f_temp_difference = a_rect->m_max[index] - a_rect->m_min[index];
			//#ifdef _DEBUG
			//			cout<<" temp difference: "<< f_temp_difference <<endl;
			//#endif
						//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m _min[index] << endl;


						//if (f_temp_difference < DBL_TRUE_MIN) continue;//191128 will influence result
			volume += f_temp_difference;

		}
		else {//even min id
			//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m _min[index] << endl;
			//ASSERT(a_rect->m_max[index] == a_rect->m_ min[index]);
			volume += (a_rect->m_max[index] - a_rect->m_min[index - 2]);

		}
		/*-------------------------------------------------------------------------------------------*/
		/*--------------------------------------New version-----------------------------------------*/
		//volume += (a_rect->m_max[index] - a_rect->m _min[index]);
		/*-------------------------------------------------------------------------------------------*/
		//cout << volume << endl;
	}


#ifdef _DEBUG
	//if (volume < DBL_TRUE_MIN) {
	//	for (int index = 1; index < NUMDIMS; index++) {
	//		if (index & 1) {//odd
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m _min[index] << endl;
	//		}
	//		else {//even
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index - 2 << "] = " << a_rect->m _min[index-2] << endl;
	//		}
	//	}
	//}
	ASSERT(!isinf(volume));
	ASSERT(volume >= DBL_TRUE_MIN);
	ASSERT(isfinite(volume));
	ASSERT(isnormal(volume));
	ASSERT(volume > (ELEMTYPEREAL)0);
#endif
	//cout << volume << endl;
	return fabs(volume);
}

// Calculate the n-dimensional volume of a rectangle      : can get the area of rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::PLARectVolume(Rect* a_rect)
{
	//PrintRect(*a_rect);
#ifdef _DEBUG
	assert(0);
	ASSERT(a_rect);
#endif
	ELEMTYPEREAL  volume = (ELEMTYPEREAL)1;
	volume = 0;
	int count = 0;

	for (int index = 0; index < NUMDIMS; ++index)
	{
		ASSERT(a_rect->m_max[index] >= a_rect->m_min[index]);
		if (a_rect->m_max[index] == a_rect->m_min[index]) {
			count++;
			continue;
		}
		cout << "m_max: " << a_rect->m_max[index] << " m_min: " << a_rect->m_min[index] << endl;
		volume *= (a_rect->m_max[index] - a_rect->m_min[index]);// *10.0;
		cout << "   volume: " << volume << endl;
#ifdef _DEBUG
		if (volume == 0) {
			cout << a_rect->m_min[index] << endl;
			cout << a_rect->m_max[index] << endl;
			cout << "   volume: " << volume << endl;
		}
#endif
	}
#ifdef _DEBUG
	ASSERT(!isinf(volume));
	ASSERT(volume >= DBL_TRUE_MIN);
	ASSERT(isfinite(volume));
	ASSERT(isnormal(volume));
	//cout << "volume: " << volume << endl;
	ASSERT(volume > (ELEMTYPEREAL)0);
#endif
	if (count == NUMDIMS) volume = 0;
	return volume;
}

// Calculate the n-dimensional volume of a rectangle For PLA & CHEBY      : can get the area of rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::ChebyshevRectVolume(Rect* a_rect)
{
	//PrintRect(*a_rect);
	ASSERT(a_rect);
	/*--------------------------------------Old version-----------------------------------------*/
	ELEMTYPEREAL  volume = (ELEMTYPEREAL)1;
	ELEMTYPEREAL  volume0 = (ELEMTYPEREAL)1;
	/*------------------------------------------------------------------------------------------*/
	/*--------------------------------------new version-----------------------------------------*/
	//volume = 0;
	/*------------------------------------------------------------------------------------------*/
	int count = 0;

	for (int index = 0; index < NUMDIMS; ++index)
	{
		//cout << "m_max: " << a_rect->m_max[index] << " m _min: " << a_rect->m _min[index] << " max-min:" << a_rect->m_max[index] - a_rect->m _min[index] << endl;
		//ASSERT(a_rect->m_max[index] != 0 && a_rect->m _min[index]!=0);
#ifdef _DEBUG
		ASSERT(a_rect->m_max[index] >= a_rect->m_min[index]);
		if (a_rect->m_max[index] == a_rect->m_min[index]) {
			//cout << "m_max: " << a_rect->m_max[index] << " m _min: " << a_rect->m_ min[index] << " max-min:" << a_rect->m_max[index] - a_rect->m _min[index] << endl;
			if (a_rect->m_max[index] != 0) {
				volume0 *= fabs(a_rect->m_max[index]);
			}

			//cout << "   volume0: " << volume0 << endl;
			count++;
			continue;
		}
#endif
		//cout << "m_max: " << a_rect->m_max[index] << " m _min: " << a_rect->m _min[index]<<" max-min:"<< a_rect->m_max[index]- a_rect->m _min[index] << endl;

		/*--------------------------------------Old version-----------------------------------------*/
		//volume *= (a_rect->m_max[index] - a_rect->m _min[index]) *10.0;
		/*------------------------------------------------------------------------------------------*/
		/*--------------------------------------new version-----------------------------------------*/
		volume += (a_rect->m_max[index] - a_rect->m_min[index]);//191118 change volume calculation of PLA & Cheby from multiple to accummulation
		/*------------------------------------------------------------------------------------------*/
		//cout << "   volume: " << volume << endl;
#ifdef _DEBUG
		ASSERT(volume != 0);
		ASSERT(!isinf(volume));
		ASSERT(volume >= DBL_TRUE_MIN);
		ASSERT(isfinite(volume));
		ASSERT(isnormal(volume));
		/*ASSERT(isnormal(volume0));
		ASSERT(isnormal(volume));*/
		ASSERT(volume > (ELEMTYPEREAL)0);
#endif
	}

	//cout << "volume: " << volume << endl;
	if (count == NUMDIMS) {
#ifdef _DEBUG
		//ASSERT(!isinf(volume0));
		//ASSERT(volume0 >= DBL_TRUE_MIN);
		//ASSERT(isfinite(volume0));
		//ASSERT(isnormal(volume0));
		//volume = 0;
#endif
		volume = 0;
	}

	return volume;
}

//***************************************************************
// Method:apla_rect_volume
// Qualifier:   apla segment volume calculation
// Input:
// Output: compute volume of APLA segment
// date:191117
// author:
//***************************************************************
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::apla_rect_volume(Rect* a_rect) {
#ifdef _DEBUG
	ASSERT(a_rect);
#endif

	ELEMTYPEREAL volume = 0;

	for (int index = 0; index < NUMDIMS; ++index) {
#ifdef _DEBUG
		assert(a_rect->m_max[index] >= a_rect->m_min[index]);
#endif
		volume += a_rect->m_max[index] - a_rect->m_min[index];
	}

#ifdef _DEBUG
	ASSERT(!isinf(volume));
	ASSERT(volume >= DBL_TRUE_MIN);
	ASSERT(isfinite(volume));
	ASSERT(isnormal(volume));
	ASSERT(volume > (ELEMTYPEREAL)0);
	//cout <<"volume: "<< volume << endl;
#endif

	return volume;
}

//***************************************************************
// Method: partition_rect_volume
// Qualifier:   
// Input: 1 rectangle_new: sub rectangle. 2 node: leaf node.
// Output: compute volume of APLA apprximation.
// date:210719
// author:
//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y>
long double RTREE_PARTITION_QUAL::partition_rect_volume(const T& const rectangle_new, Y& const node) {
	assert(0);
	long double volume = 0;
//	/*...........................................................................................................................................................*/
//#ifdef _DEBUG
//	//ASSERT(a_rect);
//	assert(node->m_count > 1 && node->IsLeaf());
//	assert(rectangle_new.original_time_series_vector_pointer->size() == rectangle_new.reconstruct_time_series_vector_pointer->size() && rectangle_new.list_partition_pointer->size() >= rectangle_new.list_original_pointer->size());
//	assert(rectangle_new.original_time_series_vector_pointer->size() == node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size() && node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size() == node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.size());
//	assert(rectangle_new.original_time_series_vector_pointer->size() == rectangle_new.list_partition_pointer->back().right_endpoint);
//#endif
//	/*...........................................................................................................................................................*/
//
//	
//
//	for (int id_point = 0; id_point < rectangle_new.original_time_series_vector_pointer->size(); id_point++) {
//		node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point] = max(node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point], (*rectangle_new.reconstruct_time_series_vector_pointer)[id_point]);
//		node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point] = min(node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point], (*rectangle_new.reconstruct_time_series_vector_pointer)[id_point]);
//	}
//
//	if (node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == rectangle_new.list_partition_pointer->size()) {
//
//		/*...........................................................................................................................................................*/
//#ifdef _DEBUG
//		assert(node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->size() && node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() == node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->size());
//		assert(node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size());
//		assert(node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == rectangle_new.list_partition_pointer->size());
//#endif
//		/*...........................................................................................................................................................*/
//	}
//	else {
//		node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->clear();
//		node->rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->clear();
//		node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
//		node->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
//	}
//
//	/*...................................................................*/
//#ifdef _DEBUG
//	ASSERT(!isinf(volume));
//	ASSERT(volume >= DBL_TRUE_MIN);
//	ASSERT(isfinite(volume));
//	ASSERT(isnormal(volume));
//	ASSERT(volume > (ELEMTYPEREAL)0);
//	//cout <<"volume: "<< volume << endl;
//#endif
//	/*...................................................................*/

	return volume;
}

//***************************************************************
	// Method: scan_each_node_for_variant
	// Qualifier:   
	// Input: Node
	// Output:  Scan each sub node to get variant.
	// date:210806
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y>//210806 Scan each sub node to get variant.
Y& RTREE_PARTITION_QUAL::scan_each_node_for_variant(T& const node, Y& const rectangle) {
	/*.............................*/
#ifdef _DEBUG
	assert_node(node);
#endif
	/*.............................*/

	if (node.IsInternalNode()) {  // Internal node
		for (int id_branch = 0; id_branch < node.m_count; ++id_branch) {
			scan_each_node_for_variant(*node.m_branch[id_branch].m_child, rectangle);
		}
	}
	else { // Leaf node

		/*.............................*/
#ifdef _DEBUG
		assert(node.m_count == node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size());
		assert(node.IsLeaf());
#endif
		/*.............................*/

		for (int id_branch = 0; id_branch < node.m_count; id_branch++) {
			rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.emplace(node.m_branch[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer);
		}
	}

	return rectangle;
}

//***************************************************************
	// Method: get_convex_hull_volume_by_branch_pointer
	// Qualifier:   
	// Input: Pointer of Branch.
	// Output: Combind convex hull of these branches and volume of this convex hull.
	// date:210724
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
long double& RTREE_PARTITION_QUAL::get_convex_hull_volume_by_minmax_reconstruct_point(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined) {
	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert(size_branch > 0);
	assert(branch_pointer);
	assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*...........................................................................................................................................................*/
	assert(0);
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  Min Max Point  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*--------------------------- max(at+b)-min(at+b) ---------------------------------*/
	rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, -INF);
	/*----------------------------------------------------------------------------------*/

	rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point = 0;// max reconstruct point difference
	rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = 0;// max reconstruct point difference square root

	for (int id_point = 0; id_point < branch_pointer[0].m_rect.rectangle_entry_struct.original_time_series_vector_pointer->size(); id_point++) {

		for (int id_branch = 0; id_branch < size_branch; id_branch++) {
			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			assert(branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer->size() >= branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_original_pointer->size());
			//assert(branch_pointer[id_branch].m_rect.original_time_series_vector_pointer->size() == branch_pointer[id_branch].m_rect.reconstruct_time_series_vector_pointer->size());
			assert(rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.size() == rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size());
			//assert(rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.size() == branch_pointer[id_branch].m_rect.reconstruct_time_series_vector_pointer->size());
			assert(branch_pointer[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer->size() - 1 == branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer->back().right_endpoint);
			APLA::assert_has_same_endpoints(*branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
#endif
			/*...........................................................................................................................................................*/

			//rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point] = max(rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point], (*branch_pointer[id_branch].m_rect.reconstruct_time_series_vector_pointer)[id_point]);
			//rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point] = min(rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point], (*branch_pointer[id_branch].m_rect.reconstruct_time_series_vector_pointer)[id_point]);
		}

		/*...........................................................................................................................................................*/
#ifdef _DEBUG
		assert(rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point] >= rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point]);
#endif
		/*...........................................................................................................................................................*/

		const long double difference_point = rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point] - rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point];
		rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point += difference_point;
		rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt += powl(difference_point, 2);
	}

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*...........................................................................................................................................................*/

	return rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = sqrtl(rectangle_combined.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt);
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
}

RTREE_TEMPLATE
template<typename T>//210806 vector in node
long double& RTREE_PARTITION_QUAL::get_convex_hull_by_minmax_reconstruct_point(T& const rectangle_in_node) {
	/*.............................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_in_node);
#endif
	/*.............................*/

	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  Min Max Point  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	/*--------------------------- max(at+b)-min(at+b) ---------------------------------*/
	rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.assign(1024, INF);
	rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.assign(1024, -INF);
	/*----------------------------------------------------------------------------------*/

	rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point = 0;// max reconstruct point difference
	rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = 0;// max reconstruct point difference square root

	for (int id_point = 0; id_point < rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size(); id_point++) {

		for (int id_branch = 0; id_branch < rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(); id_branch++) {
			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			assert(rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.size() == rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size());
			assert(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[0]->size() - 1 == rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->back().right_endpoint);
			APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0], *rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch]);
#endif
			/*...........................................................................................................................................................*/

			//rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point] = max(rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point], (*rectangle_in_node.reconstruct_time_series_pointer_vector[id_branch])[id_point]);
			//rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point] = min(rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point], (*rectangle_in_node.reconstruct_time_series_pointer_vector[id_branch])[id_point]);
		}

		/*...........................................................................................................................................................*/
#ifdef _DEBUG
		assert(rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point] >= rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point]);
#endif
		/*...........................................................................................................................................................*/

		const long double difference_point = rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector[id_point] - rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector[id_point];
		rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point += difference_point;
		rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt += powl(difference_point, 2);
	}

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert_rectangle(rectangle_in_node);
#endif
	/*...........................................................................................................................................................*/

	return rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = sqrtl(rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt);
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

}

//***************************************************************
	// Method: get_convex_hull_volume_by_distance_approximation_global
	// Qualifier: Entry in branch. In Leaf Node.       max approximation distance and these two approximations
	// Input: Pointer of Branch.
	// Output: Combind convex hull of these branches and volume of this convex hull.
	// date:210724
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
long double& RTREE_PARTITION_QUAL::get_convex_hull_volume_by_distance_approximation_global(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined) {
	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert(size_branch > 0);
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	assert_branch_rect_in_node(rectangle_combined, branch_pointer);
	is_rectangle_node(rectangle_combined);
	is_branch_entry(*branch_pointer);
#endif
	/*...........................................................................................................................................................*/

	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;
	long double distance_current = -INF;
	rectangle_combined.seed_0 = rectangle_combined.seed_1 = -1;

	/*====================================   For global approximation. Bubble sort: compare every two of branches  =================================================*/
	for (int id_branch_a = 0; id_branch_a < size_branch - 1; ++id_branch_a) {
		for (int id_branch_b = id_branch_a + 1; id_branch_b < size_branch; ++id_branch_b) {

			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
#endif
			/*...........................................................................................................................................................*/

			distance_current = APLA::get_sqrt_distance_sapla_same_endpoints(*branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);

			//CombineRect(&a_parVars->m_branchBuf[id_branch_a].m_rect, &a_parVars->m_branchBuf[id_branch_b].m_rect, &oneRect);
			//waste = CalcRectVolume(&oneRect) - area[id_branch_a] - area[id_branch_b];

			//cout << "*******waste: " << waste << endl;
			if (distance_current > rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt) {

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_current;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt, 2);

				rectangle_combined.seed_0 = id_branch_a;
				rectangle_combined.seed_1 = id_branch_b;

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_0 = branch_pointer[id_branch_a].m_data;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_1 = branch_pointer[id_branch_b].m_data;
			}
			//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
		}
		//cout << endl;
	}

	/*------------  get two approximation with max SAPLA distance  -------------*/
	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*branch_pointer[rectangle_combined.seed_0].m_rect.rectangle_entry_struct.list_partition_pointer);
	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*branch_pointer[rectangle_combined.seed_1].m_rect.rectangle_entry_struct.list_partition_pointer);
	/*---------------------------------------------------------------------------*/

	/*------------  211225  -------------*/
	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = branch_pointer[rectangle_combined.seed_0].m_rect.rectangle_entry_struct.list_original_pointer;
	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = branch_pointer[rectangle_combined.seed_1].m_rect.rectangle_entry_struct.list_original_pointer;
	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = branch_pointer[rectangle_combined.seed_0].m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = branch_pointer[rectangle_combined.seed_1].m_rect.rectangle_entry_struct.original_time_series_vector_pointer;
	/*-----------------------------------*/
	/*==========================================================================================================================================================*/

	/*........................................................................*/
#ifdef _DEBUG
	assert(rectangle_combined.seed_0 != rectangle_combined.seed_1);
	//assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*........................................................................*/

	return rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;
}


//211227
RTREE_TEMPLATE
template<typename T, typename Y, typename U>// In intranl Node. 
long double& RTREE_PARTITION_QUAL::get_convex_hull_volume_by_distance_internal_node_global(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined) {
	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert(size_branch > 0);
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	assert_branch_rect_in_node(rectangle_combined, branch_pointer);
	is_rectangle_node(rectangle_combined);
#endif
	/*...........................................................................................................................................................*/

	rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;
	long double distance_current = -INF;
	rectangle_combined.seed_0 = rectangle_combined.seed_1 = -1;
	long double distance_sqrt_SAPLA_with_new_entry_max = -INF;
	/*====================================   For global approximation. Bubble sort: compare every two of branches  =================================================*/
	for (int id_branch_a = 0; id_branch_a < size_branch - 1; ++id_branch_a) {
		for (int id_branch_b = id_branch_a + 1; id_branch_b < size_branch; ++id_branch_b) {

			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			//APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
			APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_child->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *branch_pointer[id_branch_b].m_child->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
			APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_child->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, *branch_pointer[id_branch_b].m_child->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
#endif
			/*...........................................................................................................................................................*/

			//combine_compute_partition_node_by_convex(branch_pointer[id_branch_a].m_child->m_child->rectangle_in_node, branch_pointer[id_branch_b].m_child->m_child->rectangle_in_node);
			auto& const minmax_struct_a = branch_pointer[id_branch_a].m_child->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct;
			auto& const minmax_struct_b = branch_pointer[id_branch_b].m_child->rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct;
			//1
			long double distance_sqrt_min0_min1 = APLA::get_same_partition_SAPLA_limit_only_move_1(*minmax_struct_a.original_time_series_vector_min_pointer, *minmax_struct_b.original_time_series_vector_min_pointer, *minmax_struct_a.list_partition_global_min_pointer, *minmax_struct_b.list_partition_global_min_pointer, option_tree_struct);
			if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min0_min1) {
				distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min0_min1;

				/*~~~~~~~~~~~~~~~~~~~~~~ max convex~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_sqrt_SAPLA_with_new_entry_max;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(distance_sqrt_SAPLA_with_new_entry_max, 2);

				rectangle_combined.seed_0 = id_branch_a;
				rectangle_combined.seed_1 = id_branch_b;

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_0 = minmax_struct_a.seed_0;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_1 = minmax_struct_b.seed_0;

				/*------------  get two approximation with max SAPLA distance  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*minmax_struct_a.list_partition_global_min_pointer);
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*minmax_struct_b.list_partition_global_min_pointer);
				/*---------------------------------------------------------------------------*/

				/*------------  211225  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = minmax_struct_a.list_original_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = minmax_struct_b.list_original_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = minmax_struct_a.original_time_series_vector_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = minmax_struct_b.original_time_series_vector_min_pointer;
				/*-----------------------------------*/
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			}
			
			//2
			long double distance_sqrt_max0_max1 = APLA::get_same_partition_SAPLA_limit_only_move_1(*minmax_struct_a.original_time_series_vector_max_pointer, *minmax_struct_b.original_time_series_vector_max_pointer, *minmax_struct_a.list_partition_global_max_pointer, *minmax_struct_b.list_partition_global_max_pointer, option_tree_struct);
			if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_max0_max1) {
				distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_max0_max1;

				/*~~~~~~~~~~~~~~~~~~~~~~ max convex~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_sqrt_SAPLA_with_new_entry_max;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(distance_sqrt_SAPLA_with_new_entry_max, 2);

				rectangle_combined.seed_0 = id_branch_a;
				rectangle_combined.seed_1 = id_branch_b;

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_0 = minmax_struct_a.seed_1;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_1 = minmax_struct_b.seed_1;

				/*------------  get two approximation with max SAPLA distance  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*minmax_struct_a.list_partition_global_max_pointer);
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*minmax_struct_b.list_partition_global_max_pointer);
				/*---------------------------------------------------------------------------*/

				/*------------  211225  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = minmax_struct_a.list_original_max_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = minmax_struct_b.list_original_max_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = minmax_struct_a.original_time_series_vector_max_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = minmax_struct_b.original_time_series_vector_max_pointer;
				/*-----------------------------------*/
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			}
			//3
			long double distance_sqrt_min0_max0 = minmax_struct_a.volume_minmax_whole_sqrt = APLA::get_sqrt_distance_sapla_same_endpoints(*minmax_struct_a.list_partition_global_min_pointer, *minmax_struct_a.list_partition_global_max_pointer);
			if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min0_max0) {
				distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min0_max0;
				/*~~~~~~~~~~~~~~~~~~~~~~ max convex~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_sqrt_SAPLA_with_new_entry_max;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(distance_sqrt_SAPLA_with_new_entry_max, 2);

				rectangle_combined.seed_0 = id_branch_a;
				rectangle_combined.seed_1 = id_branch_b;

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_0 = minmax_struct_a.seed_0;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_1 = minmax_struct_a.seed_1;

				/*------------  get two approximation with max SAPLA distance  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*minmax_struct_a.list_partition_global_min_pointer);
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*minmax_struct_a.list_partition_global_max_pointer);
				/*---------------------------------------------------------------------------*/

				/*------------  211225  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = minmax_struct_a.list_original_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = minmax_struct_a.list_original_max_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = minmax_struct_a.original_time_series_vector_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = minmax_struct_a.original_time_series_vector_max_pointer;
				/*-----------------------------------*/
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			}

			long double distance_sqrt_min1_max1 = minmax_struct_b.volume_minmax_whole_sqrt = APLA::get_sqrt_distance_sapla_same_endpoints(*minmax_struct_b.list_partition_global_min_pointer, *minmax_struct_b.list_partition_global_max_pointer);
			if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min1_max1) {
				distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min1_max1;
				/*~~~~~~~~~~~~~~~~~~~~~~ max convex~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_sqrt_SAPLA_with_new_entry_max;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(distance_sqrt_SAPLA_with_new_entry_max, 2);

				rectangle_combined.seed_0 = id_branch_a;
				rectangle_combined.seed_1 = id_branch_b;

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_0 = minmax_struct_b.seed_0;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_1 = minmax_struct_b.seed_1;

				/*------------  get two approximation with max SAPLA distance  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*minmax_struct_b.list_partition_global_min_pointer);
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*minmax_struct_b.list_partition_global_max_pointer);
				/*---------------------------------------------------------------------------*/

				/*------------  211225  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = minmax_struct_b.list_original_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = minmax_struct_b.list_original_max_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = minmax_struct_b.original_time_series_vector_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = minmax_struct_b.original_time_series_vector_max_pointer;
				/*-----------------------------------*/
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			}
			//4
			long double distance_sqrt_min0_max1 = APLA::get_sqrt_distance_sapla_same_endpoints(*minmax_struct_a.list_partition_global_min_pointer, *minmax_struct_b.list_partition_global_max_pointer);
			if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min0_max1) {
				distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min0_max1;
				/*~~~~~~~~~~~~~~~~~~~~~~ max convex~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_sqrt_SAPLA_with_new_entry_max;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(distance_sqrt_SAPLA_with_new_entry_max, 2);

				rectangle_combined.seed_0 = id_branch_a;
				rectangle_combined.seed_1 = id_branch_b;

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_0 = minmax_struct_a.seed_0;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_1 = minmax_struct_b.seed_1;

				/*------------  get two approximation with max SAPLA distance  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*minmax_struct_a.list_partition_global_min_pointer);
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*minmax_struct_b.list_partition_global_max_pointer);
				/*---------------------------------------------------------------------------*/

				/*------------  211225  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = minmax_struct_a.list_original_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = minmax_struct_b.list_original_max_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = minmax_struct_a.original_time_series_vector_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = minmax_struct_b.original_time_series_vector_max_pointer;
				/*-----------------------------------*/
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			}
			//5
			long double distance_sqrt_min1_max0 = APLA::get_sqrt_distance_sapla_same_endpoints(*minmax_struct_b.list_partition_global_min_pointer, *minmax_struct_a.list_partition_global_max_pointer);

			if (distance_sqrt_SAPLA_with_new_entry_max < distance_sqrt_min1_max0) {
				distance_sqrt_SAPLA_with_new_entry_max = distance_sqrt_min1_max0;
				/*~~~~~~~~~~~~~~~~~~~~~~ max convex~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_sqrt_SAPLA_with_new_entry_max;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(distance_sqrt_SAPLA_with_new_entry_max, 2);

				rectangle_combined.seed_0 = id_branch_a;
				rectangle_combined.seed_1 = id_branch_b;

				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_0 = minmax_struct_b.seed_0;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.seed_1 = minmax_struct_a.seed_1;

				/*------------  get two approximation with max SAPLA distance  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*minmax_struct_b.list_partition_global_min_pointer);
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*minmax_struct_a.list_partition_global_max_pointer);
				/*---------------------------------------------------------------------------*/

				/*------------  211225  -------------*/
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer = minmax_struct_b.list_original_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer = minmax_struct_a.list_original_max_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer = minmax_struct_b.original_time_series_vector_min_pointer;
				rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer = minmax_struct_a.original_time_series_vector_max_pointer;
				/*-----------------------------------*/
				/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			}

		}
		//cout << endl;
	}

	
	/*==========================================================================================================================================================*/

	/*........................................................................*/
#ifdef _DEBUG
	assert(rectangle_combined.seed_0 != rectangle_combined.seed_1);
	//assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*........................................................................*/

	return rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;
}

//***************************************************************
	// Method: get_convex_hull_by_distance_approximation_global
	// Qualifier: max approximation distance and these two approximations
	// Input: Pointer of Branch.
	// Output: Combind convex hull of these branches and volume of this convex hull.
	// date:210724
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T>
long double& RTREE_PARTITION_QUAL::get_convex_hull_by_distance_approximation_global(T& const rectangle_in_node) {
	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	is_rectangle_node(rectangle_in_node);
	assert_rectangle(rectangle_in_node);
#endif
	/*...........................................................................................................................................................*/

	rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;
	long double distance_current = -INF;
	rectangle_in_node.seed_0 = rectangle_in_node.seed_1 = -1;

	/*====================================   For global approximation. Bubble sort: compare every two of branches  =================================================*/
	for (int id_branch_a = 0; id_branch_a < rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size() - 1; ++id_branch_a) {
		for (int id_branch_b = id_branch_a + 1; id_branch_b < rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(); ++id_branch_b) {
			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a], *rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]);
#endif
			/*...........................................................................................................................................................*/

			distance_current = APLA::get_sqrt_distance_sapla_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a], *rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]);

			//CombineRect(&a_parVars->m_branchBuf[id_branch_a].m_rect, &a_parVars->m_branchBuf[id_branch_b].m_rect, &oneRect);
			//waste = CalcRectVolume(&oneRect) - area[id_branch_a] - area[id_branch_b];

			//cout << "*******waste: " << waste << endl;
			if (distance_current > rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt) {

				rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = distance_current;
				rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = powl(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt, 2);

				rectangle_in_node.seed_0 = id_branch_a;
				rectangle_in_node.seed_1 = id_branch_b;
			}
			//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
		}
		//cout << endl;
	}

	/*------------  get two approximation with max SAPLA distance  -------------*/
	rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_in_node.seed_0]);
	rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_in_node.seed_1]);
	/*---------------------------------------------------------------------------*/

	/*==========================================================================================================================================================*/

	/*.............................................................*/
#ifdef _DEBUG
	assert(rectangle_in_node.seed_0 != rectangle_in_node.seed_1);
	//assert_branch_rect_in_node(rectangle_in_node, branch_pointer);
#endif
	/*..............................................................*/

	return rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;
}


//***************************************************************
	// Method: get_convex_hull_volume_by_distance_approximation_local
	// Qualifier: max approximation distance and these two approximations
	// Input: Pointer of Branch.
	// Output: Combind convex hull of these branches and volume of this convex hull.
	// date:210724
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
long double& RTREE_PARTITION_QUAL::get_convex_hull_volume_by_distance_approximation_local(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined) {
	/*=====================================   For local approximation. Bubble sort: compare every two of branches   ================================================*/

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert(size_branch > 0);
	if (rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size()) {
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer);
	}
	else {
		assert(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->empty() && !branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->empty());
	}
	//assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*...........................................................................................................................................................*/
	long double distance_power_segment_current = -INF;
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.assign(branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size(), -INF);
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.assign(branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size(), -INF);
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.assign(branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size(), -1);
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.assign(branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size(), -1);

	for (int id_segment = 0; id_segment < branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size(); id_segment++) {

		for (int id_branch_a = 0; id_branch_a < size_branch - 1; ++id_branch_a) {
			for (int id_branch_b = id_branch_a + 1; id_branch_b < size_branch; ++id_branch_b) {

				/*...........................................................................................................................................................*/
#ifdef _DEBUG
				assert(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.size() == rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.size());
				APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
				APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
#endif
				/*...........................................................................................................................................................*/

				distance_power_segment_current = APLA::get_segment_distance_LB(branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_segment), branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_segment));

				if (distance_power_segment_current > rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_segment]) {

					rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_segment] = distance_power_segment_current;
					rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector[id_segment] = sqrtl(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_segment]);

					rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_segment] = id_branch_a;
					rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_segment] = id_branch_b;

				}
				//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
			}
			//cout << endl;
		}

		/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
		if (rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size()) {
			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			APLA::assert_has_same_endpoints(*branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
#endif
			/*...........................................................................................................................................................*/

			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->get(id_segment) = branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_segment]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_segment);
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->get(id_segment) = branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_segment]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_segment);
		}
		else {
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->emplace_back(branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_segment]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_segment));
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->emplace_back(branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_segment]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_segment));
		}
		/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
	}
	/*==========================================================================================================================================================*/

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer);
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	//assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*...........................................................................................................................................................*/

	//node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = accumulate(node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.begin(), node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.end(), 0);
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment = accumulate(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.begin(), rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.end(), 0);
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = sqrtl(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment);

	return rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt;
}

//***************************************************************
	// Method: get_convex_hull_by_distance_approximation_local
	// Qualifier: max approximation distance and these two approximations. vector in node
	// Input: Pointer of Branch.
	// Output: Combind convex hull of these branches and volume of this convex hull.
	// date:210806
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T>
long double& RTREE_PARTITION_QUAL::get_convex_hull_by_distance_approximation_local(T& const rectangle_in_node) {
	/*=====================================   For local approximation. Bubble sort: compare every two of branches   ================================================*/

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	if (rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size()) {
		APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]);
	}
	else {
		assert(rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->empty() && !rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->empty());
	}
	assert_rectangle(rectangle_in_node);
#endif
	/*...........................................................................................................................................................*/
	long double distance_power_segment_current = -INF;
	rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.assign(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size(), -INF);
	rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.assign(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size(), -INF);
	rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.assign(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size(), -1);
	rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.assign(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size(), -1);

	for (int id_segment = 0; id_segment < rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size(); id_segment++) {

		for (int id_branch_a = 0; id_branch_a < rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size() - 1; ++id_branch_a) {
			for (int id_branch_b = id_branch_a + 1; id_branch_b < rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size(); ++id_branch_b) {

				/*...........................................................................................................................................................*/
#ifdef _DEBUG
				assert(rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.size() == rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.size());
				APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
#endif
				/*...........................................................................................................................................................*/

				distance_power_segment_current = APLA::get_segment_distance_LB(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a]->get(id_segment), rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]->get(id_segment));

				if (distance_power_segment_current > rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_segment]) {

					rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_segment] = distance_power_segment_current;
					rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector[id_segment] = sqrtl(rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_segment]);

					rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_segment] = id_branch_a;
					rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_segment] = id_branch_b;

				}
				//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
			}
			//cout << endl;
		}

		/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
		if (rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size()) {
			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0], *rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
#endif
			/*...........................................................................................................................................................*/

			rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->get(id_segment) = rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_segment]]->get(id_segment);
			rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->get(id_segment) = rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_segment]]->get(id_segment);
		}
		else {
			rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->emplace_back(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_segment]]->get(id_segment));
			rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->emplace_back(rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_segment]]->get(id_segment));
		}
		/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
	}
	/*==========================================================================================================================================================*/

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]);
	APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
	APLA::assert_has_same_endpoints(*rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	//assert_branch_rect_in_node(rectangle_in_node, branch_pointer);
#endif
	/*...........................................................................................................................................................*/

	//node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = accumulate(node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.begin(), node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.end(), 0);
	rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment = accumulate(rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.begin(), rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.end(), 0);
	rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = sqrtl(rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment);

	return rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt;
}


// return sqrt distance
RTREE_TEMPLATE
template<typename T, typename Y>//210811 vector in node
long double RTREE_PARTITION_QUAL::compute_distance_max_between_rectangle_node(T& const rectangle_in_node_0, Y& const rectangle_in_node_1) {

	/*.................................................*/
#ifdef _DEBUG
	assert(is_rectangle_node(rectangle_in_node_0));
	assert(is_rectangle_node(rectangle_in_node_1));
	assert_rectangle(rectangle_in_node_0);
	assert_rectangle(rectangle_in_node_1);
#endif
	/*..................................................*/

	long double distance_max = -INF;
	long double distance_max_sqrt = -INF;
	long double distance_current = -INF;

	for (int id_branch_a = 0; id_branch_a < rectangle_in_node_0.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(); ++id_branch_a) {
		for (int id_branch_b = 0; id_branch_b < rectangle_in_node_1.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(); ++id_branch_b) {
			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			APLA::assert_has_same_endpoints(*rectangle_in_node_0.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a], *rectangle_in_node_1.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]);
#endif
			/*...........................................................................................................................................................*/

			distance_current = APLA::get_sqrt_distance_sapla_same_endpoints(*rectangle_in_node_0.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a], *rectangle_in_node_1.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]);

			//CombineRect(&a_parVars->m_branchBuf[id_branch_a].m_rect, &a_parVars->m_branchBuf[id_branch_b].m_rect, &oneRect);
			//waste = CalcRectVolume(&oneRect) - area[id_branch_a] - area[id_branch_b];

			//cout << "*******waste: " << waste << endl;
			if (distance_current > distance_max_sqrt) {

				distance_max_sqrt = distance_current;
				distance_max = powl(distance_max_sqrt, 2);

				//rectangle_in_node.seed_0 = id_branch_a;
				//rectangle_in_node.seed_1 = id_branch_b;
			}
			//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
		}
		//cout << endl;
	}

	return distance_max_sqrt;
}

//***************************************************************
	// Method: get_convex_hull_volume_by_MBR_local
	// Qualifier: max approximation distance and these two approximations
	// Input: Pointer of Branch.
	// Output: Combind MBR convex hull of these branches and volume of MBR.
	// date:211011
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y, typename U>//211011 entry in branch. In Leaf Node. 
long double& RTREE_PARTITION_QUAL::get_convex_hull_volume_by_MBR_local(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined) {
	/*...............................................................*/
#ifdef _DEBUG
	assert(size_branch > 0 && size_branch <= this->TMAXNODESs);
	assert_branch_rect_in_node(rectangle_combined, branch_pointer);
	is_rectangle_node(rectangle_combined);
	is_branch_entry(*branch_pointer);
#endif
	/*................................................................*/

	const size_t& const size_list_partition = branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size();
	long double distance_power_segment_current = -INF;

	const vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>& const max_vector = rectangle_combined.rectangle_inter_node_struct.minmax_MBR_struct.MBR_list_max_vector;
	const vector<DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>>& const min_vector = rectangle_combined.rectangle_inter_node_struct.minmax_MBR_struct.MBR_list_min_vector;
	max_vector.assign(size_list_partition, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>());
	min_vector.assign(size_list_partition, DoublyLinkedList<APLA::AREA_COEFFICIENT_SPEED_NO_MINMAX>());

	for (int id_seg = 0; id_seg < size_list_partition; id_seg++) {

		/*..................................................................*/
#ifdef _DEBUG
		assert(max_vector[id_seg].empty() && min_vector[id_seg].empty());
#endif
		/*..................................................................*/

		max_vector[id_seg].emplace_back(branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg));
		min_vector[id_seg].emplace_back(branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg));

		for (int id_branch_a = 1; id_branch_a < size_branch - 1; ++id_branch_a) {
			//for (int id_branch_b = id_branch_a + 1; id_branch_b < size_branch; ++id_branch_b) {
				/*...........................................................................................................................................................*/
#ifdef _DEBUG
				//APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
#endif
				/*...........................................................................................................................................................*/

				//distance_power_segment_current = APLA::get_segment_distance_LB(branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg), branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg));

				if (distance_power_segment_current > rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_seg]) {

					rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_seg] = distance_power_segment_current;
					rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector[id_seg] = sqrtl(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector[id_seg]);

					rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_seg] = id_branch_a;
					//rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_seg] = id_branch_b;

				}
				//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
			//}
			//cout << endl;
		}

		/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
		if (rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() == size_list_partition) {
			/*...........................................................................................................................................................*/
#ifdef _DEBUG
			APLA::assert_has_same_endpoints(*branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
#endif
			/*...........................................................................................................................................................*/

			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->get(id_seg) = branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_seg]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg);
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->get(id_seg) = branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_seg]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg);
		}
		else {
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->emplace_back(branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector[id_seg]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg));
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->emplace_back(branch_pointer[rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector[id_seg]].m_rect.rectangle_entry_struct.list_partition_pointer->get(id_seg));
		}
		/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
	}
	/*==========================================================================================================================================================*/

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer);
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
	APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	//assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*...........................................................................................................................................................*/

	//node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = accumulate(node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.begin(), node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.end(), 0);
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment = accumulate(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.begin(), rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.end(), 0);
	rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = sqrtl(rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment);

	return rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt;

}

//***************************************************************
	// Method: get_convex_hull_by_entry
	// Qualifier: Leaf Node
	// Input:1 Branch pointer: ; 2 size_branch: size of branch; 3 rectangle_combined: Rectangle;
	// Output: compute convex hull of leaf node and Volume of convex hull.
	// date:210726
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T, typename Y, typename U>
long double RTREE_PARTITION_QUAL::get_convex_hull_by_entry(const T* const branch_pointer, const Y& const size_branch, U& const rectangle_combined) {
	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	//ASSERT(a_rect);
	assert(branch_pointer[0].m_rect.rectangle_entry_struct.list_original_pointer);
	assert(size_branch > 0);
	assert(branch_pointer);
	assert_branch_rect_in_node(rectangle_combined, branch_pointer);
	is_rectangle_node(rectangle_combined);
	is_branch_entry(*branch_pointer);
	//assert(rectangle_new.original_time_series_vector_pointer.size() == rectangle_new.reconstruct_time_series_vector_pointer.size() && rectangle_new.list_partition_pointer->size() >= rectangle_new.list_original_pointer->size());
	//assert(rectangle_new.original_time_series_vector_pointer.size() == node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size() && node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size() == node->rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.size());
	//assert(rectangle_new.original_time_series_vector_pointer.size() == rectangle_new.list_partition_pointer->back().right_endpoint);
#endif
	/*...........................................................................................................................................................*/

	switch (option_tree_struct.type_volume) {//Leaf Node-> Entry
	case 0: {// global distance
		/*...........................................................................................................................................................*/
#ifdef _DEBUG
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
#endif
		/*...........................................................................................................................................................*/
		
		rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;
		if (rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() != branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size()) {
			rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
			rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
		}

		/*=========   Leaf Node .For global approximation. Bubble sort: compare every two of branches  ============*/
		get_convex_hull_volume_by_distance_approximation_global(branch_pointer, size_branch, rectangle_combined);
		/*=========================================================================================================*/
		break;
	}
	default: {//Leaf Node
		assert(0);
		/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  Min Max Point  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
		get_convex_hull_volume_by_minmax_reconstruct_point(branch_pointer, size_branch, rectangle_combined);
		/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

		/*...........................................................................................................................................................*/
#ifdef _DEBUG
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
#endif
		/*...........................................................................................................................................................*/

		rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;

		if (rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() != branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer->size()) {
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->clear();
			rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->clear();
			rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
			rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
		}

		/*==============   For global approximation. Bubble sort: compare every two of branches  ======================*/
		get_convex_hull_volume_by_distance_approximation_global(branch_pointer, size_branch, rectangle_combined);
		/*=============================================================================================================*/

		/*=================   For local approximation. Bubble sort: compare every two of branches   ====================*/
		get_convex_hull_volume_by_distance_approximation_local(branch_pointer, size_branch, rectangle_combined);
		/*==============================================================================================================*/

		/*...........................................................................................................................................................*/
#ifdef _DEBUG
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *branch_pointer[0].m_rect.rectangle_entry_struct.list_partition_pointer);
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
		APLA::assert_has_same_endpoints(*rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
#endif
		/*...........................................................................................................................................................*/
		break;
	}
	}

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*...........................................................................................................................................................*/

	return rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;
}

//***************************************************************
// Method: get_convex_hull_by_node
// Qualifier:   
// Input:1 node: leaf node.
// Output: compute convex hull of leaf node.
// date:210720
// author:
//***************************************************************
RTREE_TEMPLATE
template<typename T>
long double RTREE_PARTITION_QUAL::get_convex_hull_by_node(T& const node) {
	/*............................................................*/
#ifdef _DEBUG
	assert_node(node);
#endif
	/*............................................................*/
	if (node.m_count == 1) return 0;

	if (node.IsLeaf()) {// Leaf Node from Entry
		/*............................................................*/
#ifdef _DEBUG
		assert(node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.size_entry() == node.m_count && node.m_count > 1);
#endif
		/*............................................................*/
		//For leaf node
		return get_convex_hull_by_entry(node.m_branch, node.m_count, node.rectangle_in_node);
	}
	else if (node.IsInternalNode()) {// Internal Node for sub nodes
		/*............................................................*/
#ifdef _DEBUG
		assert(node.m_count <= node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.size_entry());
#endif
		/*............................................................*/
		//For internal node
		get_convex_hull_by_node_rectangle(node);
	}
	else {
		return 0;
	}

}

//***************************************************************
	// Method: get_convex_hull_by_node_rectangle
	// Qualifier:   
	// Input:1 node: internal node.
	// Output: compute convex hull of internal node.
	// date:210804
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T>
long double RTREE_PARTITION_QUAL::get_convex_hull_by_node_rectangle(T& const node) {

	/*...   210804   ....*/
#ifdef _DEBUG
	assert_node(node);
#endif
	/*...................*/

	get_convex_hull_by_rectangle(node.rectangle_in_node);

	//switch (option_tree_struct.type_volume) {
	//case 0: {// global distance
	//	node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;
	//	if (node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() != node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size()) {
	//		node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
	//		node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
	//	}

	//	/*====================================   For global approximation. Bubble sort: compare every two of branches  =============================================*/
	//	get_convex_hull_by_distance_approximation_global(node.rectangle_in_node);
	//	/*==========================================================================================================================================================*/
	//	break;
	//}
	//default: {
	//	assert(0);
	//	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  Min Max Point  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	//	get_convex_hull_by_minmax_reconstruct_point(node.rectangle_in_node);
	//	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	//	node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;

	//	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	//	if (node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() != node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size()) {
	//		node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->clear();
	//		node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->clear();
	//		node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
	//		node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
	//	}

	//	/*====================================   For global approximation. Bubble sort: compare every two of branches  =============================================*/
	//	get_convex_hull_by_distance_approximation_global(node.rectangle_in_node);
	//	/*==========================================================================================================================================================*/

	//	/*=====================================   For local approximation. Bubble sort: compare every two of branches   ============================================*/
	//	get_convex_hull_by_distance_approximation_local(node.rectangle_in_node);
	//	/*==========================================================================================================================================================*/

	//	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	//	break;
	//}
	//}

	/*....................*/
#ifdef _DEBUG
	assert_node(node);
#endif
	/*.....................*/

	return node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;

}

//***************************************************************
	// Method: get_convex_hull_by_node_rectangle
	// Qualifier:   
	// Input:1 node: internal node.
	// Output: compute convex hull of internal node.
	// date:210810
	// author:
	//***************************************************************
RTREE_TEMPLATE
template<typename T>
long double RTREE_PARTITION_QUAL::get_convex_hull_by_rectangle(T& const rectangle) {
	/*........   210804   ......*/
#ifdef _DEBUG
	assert_rectangle(rectangle);
#endif
	/*..........................*/

	switch (option_tree_struct.type_volume) {
	case 0: {// global distance

		if (rectangle.rectangle_inter_node_struct.vector_vector_struct.size_entry() == 0) {
			return 0;
		}
		else if (rectangle.rectangle_inter_node_struct.vector_vector_struct.size_entry() == 1) {
			return rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = 0;
		}

		rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;

		if (rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() != rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size()) {
			rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
			rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
		}

		/*======   For global approximation. Bubble sort: compare every two of branches  =======*/
		get_convex_hull_by_distance_approximation_global(rectangle);
		/*=======================================================================================*/
		
		break;
	}
	default: {
		assert(0);

		if (rectangle.rectangle_inter_node_struct.vector_vector_struct.size_entry() == 0) {
			return 0;
		}
		else if (rectangle.rectangle_inter_node_struct.vector_vector_struct.size_entry() == 1) {
			return rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = rectangle.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt = rectangle.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment = rectangle.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt = rectangle.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point = 0;
		}

		/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  Min Max Point  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
		get_convex_hull_by_minmax_reconstruct_point(rectangle);
		/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

		rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt = rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole = -INF;

		if (rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() != rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[0]->size()) {
			rectangle.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->clear();
			rectangle.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer->clear();
			rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->clear();
			rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->clear();
		}

		/*==   For global approximation. Bubble sort: compare every two of branches  ====*/
		get_convex_hull_by_distance_approximation_global(rectangle);
		/*================================================================================*/

		/*======= For local approximation. Bubble sort: compare every two of branches ====*/
		get_convex_hull_by_distance_approximation_local(rectangle);
		/*================================================================================*/
		break;
	}
	}


	/*....................*/
#ifdef _DEBUG
	assert_rectangle(rectangle);
#endif
	/*.....................*/
	/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	return rectangle.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt;
}

// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::RectSphericalVolume(Rect* a_rect)
{

	ASSERT(a_rect);

	ELEMTYPEREAL sumOfSquares = (ELEMTYPEREAL)0;
	ELEMTYPEREAL radius;

	for (int index = 0; index < NUMDIMS; ++index)
	{
		ELEMTYPEREAL halfExtent = ((ELEMTYPEREAL)a_rect->m_max[index] - (ELEMTYPEREAL)a_rect->m_min[index]) * 0.5f;
		sumOfSquares += halfExtent * halfExtent;
	}

	radius = (ELEMTYPEREAL)sqrt(sumOfSquares);

	// Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
	if (NUMDIMS == 3)
	{
		return (radius * radius * radius * m_unitSphereVolume);
	}
	else if (NUMDIMS == 2)
	{
		return (radius * radius * m_unitSphereVolume);
	}
	else
	{
		//cout << (ELEMTYPEREAL)(pow(radius, NUMDIMS) * m_unitSphereVolume) << endl;
		return (ELEMTYPEREAL)(pow(radius, NUMDIMS) * m_unitSphereVolume);
	}
}


// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::APCARectSphericalVolume(Rect* a_rect)
{

	ASSERT(a_rect);

	ELEMTYPEREAL sumOfSquares = (ELEMTYPEREAL)0;
	ELEMTYPEREAL radius;
	ELEMTYPEREAL halfExtent = ((ELEMTYPEREAL)(a_rect->m_max[0] + 1) * 0.5f);
	//cout <<endl<< a_rect->m_max[0]+1 << " - " << 0 << endl;
	sumOfSquares += halfExtent * halfExtent;

	for (int index = 1; index < NUMDIMS; ++index)
	{
		if (index & 1) {
			halfExtent = ((ELEMTYPEREAL)a_rect->m_max[index] - (ELEMTYPEREAL)a_rect->m_min[index]) * 0.5;
			//cout << a_rect->m_max[index] << " - " << a_rect->m _min[index] << endl;
		}
		else {
			halfExtent = ((ELEMTYPEREAL)a_rect->m_max[index] - (ELEMTYPEREAL)a_rect->m_min[index - 2]) * 0.5;
			//cout << a_rect->m_max[index] << " - " << a_rect->m _min[index-2] << endl;
		}

		sumOfSquares += halfExtent * halfExtent;
	}
	//cout << "sumOfSquares: " << sumOfSquares << endl;
	radius = (ELEMTYPEREAL)sqrt(sumOfSquares);
	//cout << "raius: " << radius << endl;

	// Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
	if (NUMDIMS == 3)
	{
		return (radius * radius * radius * m_unitSphereVolume);
	}
	else if (NUMDIMS == 2)
	{
		return (radius * radius * m_unitSphereVolume);
	}
	else
	{
		//cout << "Volume: " << pow(radius, NUMDIMS) * m_unitSphereVolume << endl;
		return (ELEMTYPEREAL)(pow(radius, NUMDIMS) * m_unitSphereVolume);
	}
}


// Use one of the methods to calculate retangle volume
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_PARTITION_QUAL::CalcRectVolume(Rect* a_rect)
{
#ifdef RTREE_USE_APCA_SPHERICAL_VOLUME
	return APCARectSphericalVolume(a_rect); // APCA data point. NUMDIMS < 51 ! Slower but helps certain merge cases 
#else // RTREE_USE_SPHERICAL_VOLUME

	/*..............................*/
#ifdef _DEBUG
	assert(representation_option);
#endif
	/*..............................*/

	//return RectVolume(a_rect);
	//return PLARectVolume(a_rect);
	if (representation_option == 1 || representation_option == 2) {//representation_option=1 for 1 PAA, 2 APCA 180916
		//return APCARectVolume(a_rect);          // APCA data point. NUMDIMS < 51   ! Faster, worse, NUMDIMS can be very big.
		return RectVolume(a_rect);//191120 :  Rtree, though insert MBR, the min id == max id. For Root and internal node, the min id may != max id. While for leaf node, min id == max id. So for Rtree volume calculation, do not need to worry min id == max id. So use volume *= max[i] - min[i] to calculate rectangle volume
		//return apla_rect_volume(a_rect);//191128 accumulation instead multiple
		//return APCARectVolume(a_rect);//191128 multiple min id - (min id -2)
		//return apca_rect_accumulation_volume(a_rect);//191128
	}
	else if (representation_option == 3 || representation_option == 4) {//representation_option=2 for PLA, CHEBY 180916
		//return ChebyshevRectVolume(a_rect);
		return RectVolume(a_rect);//191121 : Rtree.
		//return apla_rect_volume(a_rect);//191128 accumulation instead multiple
		//return APCARectVolume(a_rect);//191128 multiple min id - (min id -2)
	}
	else if (representation_option == 5 || representation_option == 6) {//5 APLA, 6 ICDE07 191117
		return RectVolume(a_rect);//191128 original Rtree volume calculation, *, max - max ,min -min
		//return apla_rect_volume(a_rect);//191128 accumulation instead multiple
		//return APCARectVolume(a_rect);//191128 multiple min id - (min id -2)
		//return apca_rect_accumulation_volume(a_rect);//191128
	}
	else if (representation_option == 8) {//210618. 8 SAPLA
		return RectVolume(a_rect);
	}
	else {
		return RectVolume(a_rect);
		//assert(0);
	}

#endif // RTREE_USE_SPHERICAL_VOLUME  
}


// Load branch buffer with branches from full node plus the extra new branch.
//************************************
// Method:GetBranches
// Qualifier: 1 Load branch buffer with branches from full node plus the extra branch for Split Node() function.
//            2 Compute combined rectangle and volume
//            3 This function GetBranches() only used in Split Node() function
//            4 Used in SplitNode()
// Input: a_node: original node with full branches. a_branch: new branch. a_parVars: structure to store variants of old node and new branch.
// Output: 1 Put each branch of node and new branch in a_parVars.
//         2 Compute combined voume and get combined MBR
// date: 210714
// author:
//************************************
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::GetBranches(Node* a_node, const Branch* a_branch, PartitionVars* a_parVars)
{
	ASSERT(a_node);
	ASSERT(a_branch);
	ASSERT(a_node->m_count == TMAXNODES);

	/*+++++++++++++++++++++++++++++++      210726 Copy each branchs of node and new branch in a_parVars     ++++++++++++++++++++++++++++++++++*/
	// Load the branch buffer
	for (int index = 0; index < TMAXNODES; ++index)
	{
		/*......................   210706   .....................*/
#ifdef _DEBUG
		const Branch& test_branch0 = a_parVars->m_branchBuf[index];
		const Branch& test_branch1 = a_node->m_branch[index];
#endif
		/*.......................................................*/

		a_parVars->m_branchBuf[index] = a_node->m_branch[index];
		//assertRect(a_parVars->m_branchBuf[index].m_rect);

		/*==============  partition 210618 Copy copy all rects in sub branches =========================*/
		switch (representation_option) {
		case 2://PLA
		case 3://APCA
		case 4://PAA
		case 7://PAALM
		case 9://SAX
		case 5://CHEBY
		case 6: // ICDE07
		case 8: {// Initial 200706
			// a_parVars->m_coverSplit copy all rects in sub branches
			copy_rect_original_series_approximation_vector_pointer(a_node->m_branch[index].m_rect, a_parVars->m_coverSplit);
			a_parVars->m_branchBuf[index].branch_to_parent_node_pointer = nullptr;
			/*..........................      210729    ..........................*/
#ifdef _DEBUG
			assert_rectangle(a_node->m_branch[index].m_rect);
			assert_rectangle(a_parVars->m_coverSplit);
#endif
			/*....................................................................*/
			break;
		}
		default:
			break;
		}
		/*================================================================================*/
	}

	/*..........................      210706    ..........................*/
#ifdef _DEBUG
	const Branch& test_branch2 = a_parVars->m_branchBuf[TMAXNODES];
	const Branch& test_branch3 = *a_branch;
#endif
	/*....................................................................*/

	a_parVars->m_branchBuf[TMAXNODES] = *a_branch;

	a_parVars->m_branchCount = TMAXNODES + 1;

	/*==============  partition 210618 =========================*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		// copy coefficients
		copy_rect_original_series_approximation_vector_pointer(a_branch->m_rect, a_parVars->m_coverSplit);
		// get Voulme of combined convex hull
		//a_parVars->m_coverSplitArea = get_convex_hull_by_rectangle(a_parVars->m_coverSplit);
		/*..........................      210729    ..........................*/
#ifdef _DEBUG
		assert_rectangle(a_branch->m_rect);
		assert_rectangle(a_parVars->m_coverSplit);
#endif
		/*....................................................................*/
		break;
	}
	default:
		break;
	}
	/*==========================================================*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	// Calculate rect containing all in the set
	//o a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;

	/*+++++++++++++++++++++++            210715 This part for combined volume calculation.          +++++++++++++++++++++++++++++++++++++++++++*/
	
	/*======================================           partition 210618           ===============================*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		// get Voulme of combined convex hull

//		a_parVars->m_coverSplitArea = get_convex_hull_by_rectangle(a_parVars->m_coverSplit);;
//		/*..........................      210729    ..........................*/
//#ifdef _DEBUG
//		assert_rectangle(a_branch->m_rect);
//		assert_rectangle(a_parVars->m_coverSplit);
//#endif
//		/*....................................................................*/
		break;
	}
	default:{// Original Rtree
		assert(0);
		// Combined Rectangle (a_parVars->m_coverSplit) for combined Volume(a_parVars->m_coverSplitArea) from (original node + new branch)
	// a_parVars->m_coverSplit only for combined volume calculation.

		copy(a_parVars->m_branchBuf[0].m_rect.m_min.begin(), a_parVars->m_branchBuf[0].m_rect.m_min.begin() + NUMDIMS, a_parVars->m_coverSplit.m_min.begin());
		copy(a_parVars->m_branchBuf[0].m_rect.m_max.begin(), a_parVars->m_branchBuf[0].m_rect.m_max.begin() + NUMDIMS, a_parVars->m_coverSplit.m_max.begin());

		for (int index = 1; index < TMAXNODES + 1; ++index) {
			// a_parVars->m_coverSplit = CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
			CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect, &a_parVars->m_coverSplit);
		}
		/*cout<<"Get Branches()"<<endl;
		for (int index = 0; index < NUMDIMS; ++index) {
			cout << "(min:" << a_parVars->m_coverSplit.m _min[index] << ", max:" << a_parVars->m_coverSplit.m_max[index] << ") ";
		}
		cout << endl;*/
		a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);
		//cout << "a_parVars->m_coverSplitArea: "<< a_parVars->m_coverSplitArea << endl;
		//cout << "coverSplitArea(): "<< a_parVars->m_coverSplitArea << endl;
		break;
	}
	}
	/*===========================================================================================================*/
	/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
}


// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// 
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::ChoosePartition(PartitionVars* a_parVars, int a_minFill)
{
	//cout << "***Choos ePartition()" << endl;
	ASSERT(a_parVars);

	ELEMTYPEREAL biggestDiff;//210702 biggest volume difference after cover

	// 1 group: group id 0/1; 2 chosen: id of branch; 3 betterGroup: 0/1, branch belong to wihch group;
	int group, chosen = 0, betterGroup = 0;//210702 group is which to choose 0 /1. better Group is the group to be put branch 0 /1. choose is branch id that is put into better group. 

	/*------  1 m_count[0] = 0, 2 m_area[0] =0, 3 m_partition[a_maxRects] = 0  --------*/
	InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
	/*---------------------------------------------------------------------------------*/

	PickSeeds(a_parVars);

	/*:::::::::::::::::::::::::::::::::::::::::::::::: Binary partition :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		/*=============================================================================================================*/
		choose_partition_by_distance(a_minFill, *a_parVars);
		/*=============================================================================================================*/
		break;
	}
	default: {
		assert(0);
		/*=============================================================================================================*/
	//m_count[0] is the number of branches in group0.  a_parVars->m_total is the total number of branches.
		while (((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
			&& (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill))// branch number is group 0 cannot to many that, branch number in gorup 1 < minimal branch number
			&& (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill)))
		{
			biggestDiff = (ELEMTYPEREAL)-1;
			for (int index = 0; index < a_parVars->m_total; ++index)// except seeds(two branches has biggesd differece area), other branches try in each group to compare are for group selection
			{
				// a_parVars->m_partition[index] is assigned by Classify(chosen, betterGroup, a_parVars)
				if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index])//Except two seeds branches. m_partition: assign the group id of each branch. seed0 and seed1 alreay assigend in PickSeed s(a_parVars)->Classify()
				{
					Rect* curRect = &a_parVars->m_branchBuf[index].m_rect;// this branch not in seeds and not in groups
					Rect rect0;//Combine seed0 and candidate branch(Put candidate branch in group 0)
					Rect rect1;//Combine seed1 and candidate branch(Put candidate branch in group 1)
					allocate_minmax_vector_in_rectangle(rect0);
					allocate_minmax_vector_in_rectangle(rect1);
					/*rect0.m_min = new ELEMTYPE[NUMDIMS];
					rect0.m_max = new ELEMTYPE[NUMDIMS];
					rect1.m_min = new ELEMTYPE[NUMDIMS];
					rect1.m_max = new ELEMTYPE[NUMDIMS];*/
					//Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
					//Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);

					CombineRect(curRect, &a_parVars->m_cover[0], &rect0);
					CombineRect(curRect, &a_parVars->m_cover[1], &rect1);

					ELEMTYPEREAL growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];//Area difference between combined group and original group 
					ELEMTYPEREAL growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];

					free_minmax_vector_in_rectangle(rect0);
					free_minmax_vector_in_rectangle(rect1);
					/*delete[] rect0.m _min;
					rect0.m _min = nullptr;
					delete[] rect0.m_max;
					rect0.m_max = nullptr;
					delete[] rect1.m _min;
					rect1.m _min = nullptr;
					delete[] rect1.m_max;
					rect1.m_max = nullptr;*/

					ELEMTYPEREAL diff = growth1 - growth0;//difference of their difference. the smaller, the better to combine
					if (diff >= 0)// group 0 has the biggesd area difference between original and combined
					{
						group = 0;// put this branch in group 0
					}
					else {// now diff < 0
						group = 1;// put this branch in group 1
						diff = -diff;//make difference of area differeence between two groups positive.
					}

					if (diff > biggestDiff)// First, find branch cause biggest diffrence in two groups, then put the rest branches in un fulled group
					{
						biggestDiff = diff;
						chosen = index;// the id of branches has biggest differnce
						betterGroup = group;//id of group
					}
					else if ((diff == biggestDiff) && (a_parVars->m_count[group] < a_parVars->m_count[betterGroup]))// if difference same, put branches in group that has fewer branches
					{
						chosen = index;
						betterGroup = group;
					}

				}
			}
			//         1 Assign group id(0/1) of each branch in a_parVars->m_partition[a_index].
			//         2 a_parVars->m_cover[0]: Combined rectangle for group[a_group]
			//         3 Volume of  a_parVars->m_cover[0]
			Classify(chosen, betterGroup, a_parVars);
		}

		// If one group too full, put remaining rects in the other
		if ((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total) {
			if (a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill) {
				group = 1;
			}
			else {
				group = 0;
			}

			for (int index = 0; index < a_parVars->m_total; ++index) {
				if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index]) {
					//         1 Assign group id(0/1) of each branch in a_parVars->m_partition[a_index].
					//         2 a_parVars->m_cover[0]: Combined rectangle for group[a_group]
					//         3 Volume of  a_parVars->m_cover[0]
					Classify(index, group, a_parVars);
				}
			}
		}
		/*=============================================================================================================*/
		break;
	}
	}
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	ASSERT((a_parVars->m_count[0] + a_parVars->m_count[1]) == a_parVars->m_total);
	ASSERT((a_parVars->m_count[0] >= a_parVars->m_minFill) && (a_parVars->m_count[1] >= a_parVars->m_minFill));

	/*..........................      210729    ..........................*/
#ifdef _DEBUG
	for (auto&& au : a_parVars->m_partition) {
		assert(au != -1);
	}
#endif
	/*....................................................................*/
}

//210729 choose_partition_by_ distance, already get two seeds.
//210904 Because Rtree is balance tree. all branches are same entry or sub nodes
RTREE_TEMPLATE
template<typename T, typename Y>
void RTREE_PARTITION_QUAL::choose_partition_by_distance(const T& const a_minFill, Y& const a_parVars) {
	/*..................................................................................................................................*/
#ifdef _DEBUG
	assert(TMAXNODES + 1 == a_parVars.m_total && a_minFill < TMAXNODES&& a_minFill >= 0);
	//assert_branch_rect_in_node(a_parVars.m_coverSplit, a_parVars.m_branchBuf);
	//assert_same_SAPLA(*a_parVars.m_coverSplit.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_0]);
	//assert_same_SAPLA(*a_parVars.m_coverSplit.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, *a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_1]);
#endif
	/*...................................................................................................................................*/

	
	if(!a_parVars.m_branchBuf[0].m_child){//branches are all entry
		//if (is_branch_pointer_all_entry(a_parVars.m_total, a_parVars.m_branchBuf)) {
		/*..................................................................................................................................*/
#ifdef _DEBUG
		assert(is_branch_pointer_all_entry(a_parVars.m_total, a_parVars.m_branchBuf));
		//if (option_tree_struct.type_tree == 1) {
			assert(a_parVars.m_coverSplit.seed_0 == a_parVars.seed_0 && a_parVars.m_coverSplit.seed_1 == a_parVars.seed_1);
			assert(a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_0] == a_parVars.m_branchBuf[a_parVars.m_coverSplit.seed_0].m_rect.rectangle_entry_struct.list_partition_pointer);
			assert(a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_1] == a_parVars.m_branchBuf[a_parVars.m_coverSplit.seed_1].m_rect.rectangle_entry_struct.list_partition_pointer);
			long double max_distance = APLA::get_sqrt_distance_sapla_same_endpoints(*a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_0], *a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_1]);
		//}
#endif
		/*...................................................................................................................................*/
		for (int id_branch = 0; id_branch < a_parVars.m_total; ++id_branch) {

			if (PartitionVars::NOT_TAKEN == a_parVars.m_partition[id_branch]) {
				long double distance_0 = APLA::get_sqrt_distance_sapla_same_endpoints(*a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_0], *a_parVars.m_branchBuf[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
				long double distance_1 = APLA::get_sqrt_distance_sapla_same_endpoints(*a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_1], *a_parVars.m_branchBuf[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
				/*..................................................................................................................................*/
#ifdef _DEBUG
				assert(max_distance >= distance_0 && max_distance >= distance_1);
				assert(a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch] == a_parVars.m_branchBuf[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
#endif
				/*...................................................................................................................................*/

				if (distance_0 <= distance_1) {

					if(a_parVars.m_count[0] + a_parVars.m_minFill >= a_parVars.m_total){
						a_parVars.m_partition[id_branch] = 1;
						++a_parVars.m_count[1];
					}
					else {
						a_parVars.m_partition[id_branch] = 0;
						++a_parVars.m_count[0];
					}
				}
				else if (distance_0 > distance_1) {

					if (a_parVars.m_count[1] + a_parVars.m_minFill >= a_parVars.m_total) {
						a_parVars.m_partition[id_branch] = 0;
						++a_parVars.m_count[0];
					}
					else {
						a_parVars.m_partition[id_branch] = 1;
						++a_parVars.m_count[1];
					}

					
				}
				else {
					assert(0);
				}
			}
		}
	}
	else {// branch all node
		/*..................................................................................................................................*/
#ifdef _DEBUG
		assert(is_branch_pointer_all_node(a_parVars.m_total, a_parVars.m_branchBuf));
		const long double max_distance = compute_distance_max_between_rectangle_node(a_parVars.m_branchBuf[a_parVars.seed_0].m_rect, a_parVars.m_branchBuf[a_parVars.seed_1].m_rect);
#endif
		/*...................................................................................................................................*/
		long double distance_0 = -INF;
		long double distance_1 = -INF;

		for (int id_branch = 0; id_branch < a_parVars.m_total; ++id_branch) {

			if (PartitionVars::NOT_TAKEN == a_parVars.m_partition[id_branch]) {

				//if (is_branch_entry(a_parVars.m_branchBuf[id_branch]) && is_branch_node(a_parVars.m_branchBuf[id_branch])) {// entry node
				//	assert(0);
				//}
				//else if (is_branch_entry(a_parVars.m_branchBuf[id_branch]) && is_branch_node(a_parVars.m_branchBuf[id_branch])) {//node entry
				//	assert(0);
				//}
				//else if (is_branch_node(a_parVars.m_branchBuf[id_branch]) && is_branch_node(a_parVars.m_branchBuf[id_branch])) {//node node
					distance_0 = compute_distance_max_between_rectangle_node(a_parVars.m_branchBuf[a_parVars.seed_0].m_rect, a_parVars.m_branchBuf[id_branch].m_rect);
					distance_1 = compute_distance_max_between_rectangle_node(a_parVars.m_branchBuf[a_parVars.seed_1].m_rect, a_parVars.m_branchBuf[id_branch].m_rect);
//				}
//				else if (is_branch_entry(a_parVars.m_branchBuf[id_branch]) && is_branch_entry(a_parVars.m_branchBuf[id_branch])) {//entry entry
//					assert(0);
//					/*..................................................................................................................................*/
//#ifdef _DEBUG
//					assert(a_parVars.m_coverSplit.seed_0 == a_parVars.seed_0 && a_parVars.m_coverSplit.seed_1 == a_parVars.seed_1);
//#endif
//					/*...................................................................................................................................*/
//					distance_0 = APLA::get_sqrt_distance_sapla_same_endpoints(*a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_0], *a_parVars.m_branchBuf[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
//					distance_1 = APLA::get_sqrt_distance_sapla_same_endpoints(*a_parVars.m_coverSplit.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[a_parVars.m_coverSplit.seed_1], *a_parVars.m_branchBuf[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
//				}
//				else {
//					assert(0);
//				}

				/*..................................................................................................................................*/
#ifdef _DEBUG
				assert(max_distance >= distance_0 && max_distance >= distance_1);
#endif
				/*...................................................................................................................................*/

				if (distance_0 <= distance_1) {

					if (a_parVars.m_count[0] + a_parVars.m_minFill >= a_parVars.m_total) {
						a_parVars.m_partition[id_branch] = 1;
						++a_parVars.m_count[1];
					}
					else {
						a_parVars.m_partition[id_branch] = 0;
						++a_parVars.m_count[0];
					}

				}
				else if (distance_0 > distance_1) {

					if (a_parVars.m_count[1] + a_parVars.m_minFill >= a_parVars.m_total) {
						a_parVars.m_partition[id_branch] = 0;
						++a_parVars.m_count[0];
					}
					else {
						a_parVars.m_partition[id_branch] = 1;
						++a_parVars.m_count[1];
					}
				}
				else {
					assert(0);
				}
			}
		}
	}

	/*..................................................................................................................................*/
#ifdef _DEBUG
	assert((a_parVars.m_count[0] + a_parVars.m_count[1]) == a_parVars.m_total);
	assert((a_parVars.m_count[0] >= a_parVars.m_minFill) && (a_parVars.m_count[1] >= a_parVars.m_minFill));
	for (auto&& au : a_parVars.m_partition) {
		assert(au != -1);
	}
#endif
	/*...................................................................................................................................*/
}


// Copy branches from the buffer into two nodes according to the partition.
// Used in SplitNode()
// Notice: a_parVars->m_partition[index];//0/1 group id of each branch.
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars)
{
	ASSERT(a_nodeA);
	ASSERT(a_nodeB);
	ASSERT(a_parVars);

	/*..........................      210729    ..........................*/
#ifdef _DEBUG
	for (auto&& au : a_parVars->m_partition) {
		assert(au != -1);
	}
#endif
	/*....................................................................*/

	for (int index = 0; index < a_parVars->m_total; ++index)
	{
		ASSERT(a_parVars->m_partition[index] == 0 || a_parVars->m_partition[index] == 1);

		const int targetNodeIndex = a_parVars->m_partition[index];//0/1 group id of each branch.
		Node* targetNodes[] = { a_nodeA, a_nodeB };

		// It is assured that AddBranch here will not cause a node split. 
		bool nodeWasSplit = AddBranch(&a_parVars->m_branchBuf[index], targetNodes[targetNodeIndex], NULL);
		ASSERT(!nodeWasSplit);
	}
}


// Initialize a PartitionVars structure.
// used in ChoosePartition()
RTREE_TEMPLATE
inline void RTREE_PARTITION_QUAL::InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill)
{
	/*..........................      210729    ..........................*/
#ifdef _DEBUG
	ASSERT(a_parVars);
	for (auto&& au : a_parVars->m_partition) {
		assert(au == PartitionVars::NOT_TAKEN);//-1
	}
#endif
	/*....................................................................*/

	a_parVars->m_count[0] = a_parVars->m_count[1] = 0; // 0/1 for 2 groups, count how many branches in one group. the number of notes in every group
	a_parVars->m_area[0] = a_parVars->m_area[1] = (ELEMTYPEREAL)0; //  int the are of each two group
	a_parVars->m_total = a_maxRects; // the number of the rectangle/branches. the number of total branches in m_ branchBuf == max node number + 1.
	a_parVars->m_minFill = a_minFill; // is minimum nuber of branch in one node
	//a_parVars->m _branchBuf = new Branch[TMAXNODES+1];
	//a_parVars->m_partition = new int[TMAXNODES+1];

	// PartitionVars::NOT_TAKEN is -1
	//for (int index = 0; index < a_maxRects; ++index) {
		//a_parVars->m_partition[index] = PartitionVars::NOT_TAKEN; // to store the node belong to which group? 0/1. assign the group id of each branch
	//}
}

//Select two entries to be the first elements of the groups.
//[Calculate inefficiency of grouping entries together.] For each pair of entries E1 and E2, compose a rectangle J including E1.I and E2.I. Calculate d = area(J) - area(E1.I) - area(E2.I).
// used in ChoosePartition()
//[Choose the most wasteful pair.] Choose the pair with the largest d.
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::PickSeeds(PartitionVars* a_parVars) {
	//cout << "**** PickSeed s()" << endl;
	int seed0 = 0, seed1 = 0;// id of branchs in group0 and group 1;

	/*==============  get seed 0 and seed 1 for node split =========================*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		// get seed 0 and seed 1 for node split
		/*.......................................*/
#ifdef _DEBUG
		assert(a_parVars->m_total == TMAXNODES + 1);
#endif
		/*.......................................*/
		/*~~~~~~~~~Get seeds~~~~~~~~~~~~~*/
		pick_seeds(*a_parVars);
		seed0 = a_parVars->seed_0;
		seed1 = a_parVars->seed_1;
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		break;
	}
	default:{
		assert(0);
		ELEMTYPEREAL worst, waste;
		vector<ELEMTYPEREAL> area(TMAXNODES + 1, 0);// area of each branch

		/*~~~~~~~~~~~~~~~ Compute volume of each branch ~~~~~~~~~~~~~~~~~~~~~*/
		for (int index = 0; index < a_parVars->m_total; ++index) {
			area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
			//cout << "  area[index]: " << area[index] << endl;
		}
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		worst = -a_parVars->m_coverSplitArea - 1;// covered volume of original node and new branch

		//cout << "****** worst: " << worst << endl;
		for (int indexA = 0; indexA < a_parVars->m_total - 1; ++indexA) // Bubble sort????????
		{
			for (int indexB = indexA + 1; indexB < a_parVars->m_total; ++indexB) {
				Rect oneRect;
				allocate_minmax_vector_in_rectangle(oneRect);
				//oneRect.m _min = new ELEMTYPE[NUMDIMS];
				//oneRect.m_max = new ELEMTYPE[NUMDIMS];

				//Rect oneRect = CombineRect(&a_parVars->m _branchBuf[indexA].m_rect, &a_parVars->m _branchBuf[indexB].m_rect);
				CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect, &oneRect);
				waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
				/*.......................................*/
#ifdef _DEBUG
			//cout << "combined rectangle: \n";
			//for (int segment_id = 0; segment_id < NUMDIMS; segment_id++) {
				//cout << oneRect.m_min[segment_id]<<" , "<< oneRect.m_max[segment_id] <<endl;
			//}
#endif
			/*.......................................*/
				free_minmax_vector_in_rectangle(oneRect);
				/*delete[] oneRect.m _min;
				oneRect.m _min = nullptr;
				delete[] oneRect.m_max;
				oneRect.m_max = nullptr;*/
				//cout << "*******waste: " << waste << endl;
				if (waste > worst)
				{
					worst = waste;
					seed0 = indexA;
					seed1 = indexB;
				}
				//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
			}
			//cout << endl;
		}

		area.clear();
		area.shrink_to_fit();
		break;
	}
	}
	/*==============================================================================*/

	/*=============== Put seeds in each group =======================================*/
	//         1 Assign group id(0/1) of each branch in a_parVars->m_partition[a_index].
	//         2 a_parVars->m_cover[0]: Combined rectangle for group[a_group]
	//         3 Volume of  a_parVars->m_cover[0]
	//cout << "******** seed0: "   << seed0 << endl;
	Classify(seed0, 0, a_parVars);//Combine branch[seed0] in group 0, update group0, Calculate Volume of group0. put seed 0 in group 0, update group0
	//cout << "******** seed1:  " << seed1 << endl;
	Classify(seed1, 1, a_parVars);//Combine branch[seed1] in group 1, update group1, Calculate Volume of group1. put seed 1 in group 1, update group1
	/*==============================================================================*/
}

//***************************************************************
	// Method:Classify
	// Qualifier: Put a branch in one of the groups.
	// Input: 1 a_index is id of branch; 2 a_group is id of group 0/1; 3 a_parVars: candidate branches
	// Output: Updated rectangle a_parVars->m_cover[a_index] and volume
	//         1 Assign group id(0/1) of each branch in a_parVars->m_partition[a_index].
	//         2 a_parVars->m_cover[0]: Combined rectangle for group[a_group]
	//         3 Volume of  a_parVars->m_cover[0]
	// notice: Brute force to get. Change original approximation point
	// date:210704
	// author:
	//***************************************************************
// Put a id of branch(seed) in one of the groups.
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::Classify(int a_index, int a_group, PartitionVars* a_parVars) {
	//cout << "*********Classify()" << endl;
	//cout << "seed: " << a_index << " group:"<< a_group<<"a_parVars->m_partition[a_index]: "<<a_parVars->m_partition[a_index] << endl;
	ASSERT(a_parVars);
	ASSERT(PartitionVars::NOT_TAKEN == a_parVars->m_partition[a_index]);

	/*~ a_parVars->m_partition store group id of each branch. this branch id: a_index assign to group: a_group(0/1), initial is -1~*/
	a_parVars->m_partition[a_index] = a_group;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*~~~~     Put branch[a_index] in group[a_group], update combined rect a_parVars->m_cover[a_group], the number of branches in each group      ~~~~*/
	if (a_parVars->m_count[a_group] == 0) {// No branch in group[a_group], put branch[a_index] in it
		//a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
		copy_n(a_parVars->m_branchBuf[a_index].m_rect.m_min.begin(), NUMDIMS, a_parVars->m_cover[a_group].m_min.begin());
		copy_n(a_parVars->m_branchBuf[a_index].m_rect.m_max.begin(), NUMDIMS, a_parVars->m_cover[a_group].m_max.begin());
	}
	else {// group[a_group] already has other branches, so combine branch[a_index] now.

		//a_parVars->m_cover[a_group] = CombineRect(&a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group]);
		CombineRect(&a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group], &a_parVars->m_cover[a_group]);
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*~~~~~~~  1 Update volume(a_parVars->m_cover[a_group]) of group[a_group]; 2 size of group[a_group]  ~~~~~~~~~~~~~~~~~~~~*/
	a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);// Update volume of group[a_group]
	++a_parVars->m_count[a_group]; //The number of branches of group[a_group] plus one.
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
}

RTREE_TEMPLATE
template<typename T>//Because Rtree is balance tree, branch all entry or nodes
long double RTREE_PARTITION_QUAL::pick_seeds(T& const a_parVars) {
	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert(a_parVars.m_total == TMAXNODES + 1);
	assert_branch_array(a_parVars.m_total, a_parVars.m_branchBuf);
#endif
	/*...........................................................................................................................................................*/

	long double volume_minmax_whole_sqrt = -INF;
	long double volume_minmax_whole = -INF;
	long double distance_current = -INF;
	a_parVars.seed_0 = -1;
	a_parVars.seed_1 = -1;

	if (!a_parVars.m_branchBuf[0].m_child) {// all branch are entry, is leaf node. !branch.m_child && branch.m_data != INF
		/*...........................................................................................................................................................*/
#ifdef _DEBUG
		assert(is_branch_pointer_all_entry(a_parVars.m_total, a_parVars.m_branchBuf));
#endif
		/*...........................................................................................................................................................*/

		a_parVars.m_coverSplitArea = get_convex_hull_by_entry(a_parVars.m_branchBuf, a_parVars.m_total, a_parVars.m_coverSplit);
		a_parVars.seed_0 = a_parVars.m_coverSplit.seed_0;
		a_parVars.seed_1 = a_parVars.m_coverSplit.seed_1;
	}
	/*====================================   For global approximation. Bubble sort: compare every two of branches  =================================================*/
	else {// all branch are nodes. balance tree
		/*...........................................................................................................................................................*/
#ifdef _DEBUG
		assert(is_branch_pointer_all_node(a_parVars.m_total, a_parVars.m_branchBuf));
#endif
		/*...........................................................................................................................................................*/

		for (int id_branch_a = 0; id_branch_a < a_parVars.m_total - 1; ++id_branch_a) {
			for (int id_branch_b = id_branch_a + 1; id_branch_b < a_parVars.m_total; ++id_branch_b) {
				/*...........................................................................................................................................................*/
#ifdef _DEBUG
			//APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
#endif
			/*...........................................................................................................................................................*/

				//if (is_branch_entry(a_parVars.m_branchBuf[id_branch_a]) && is_branch_node(a_parVars.m_branchBuf[id_branch_a])) {// entry node
				//	assert(0);
				//}
				//else if (is_branch_entry(a_parVars.m_branchBuf[id_branch_b]) && is_branch_node(a_parVars.m_branchBuf[id_branch_a])) {//node entry
				//	assert(0);
				//}
				//else if (is_branch_node(a_parVars.m_branchBuf[id_branch_a]) && is_branch_node(a_parVars.m_branchBuf[id_branch_b])) {//node node
					distance_current = compute_distance_max_between_rectangle_node(a_parVars.m_branchBuf[id_branch_a].m_rect, a_parVars.m_branchBuf[id_branch_b].m_rect);
				//}
				//else if (is_branch_entry(a_parVars.m_branchBuf[id_branch_a]) && is_branch_entry(a_parVars.m_branchBuf[id_branch_b])) {//entry entry
				//	assert(0);
				//	distance_current = APLA::get_sqrt_distance_sapla_same_endpoints(*a_parVars.m_branchBuf[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *a_parVars.m_branchBuf[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
				//}
				//else {
				//	assert(0);
				//}

				//CombineRect(&a_parVars->m_branchBuf[id_branch_a].m_rect, &a_parVars->m_branchBuf[id_branch_b].m_rect, &oneRect);
				//waste = CalcRectVolume(&oneRect) - area[id_branch_a] - area[id_branch_b];

				//cout << "*******waste: " << waste << endl;
				if (distance_current > volume_minmax_whole_sqrt) {

					volume_minmax_whole_sqrt = distance_current;
					volume_minmax_whole = powl(volume_minmax_whole_sqrt, 2);

					a_parVars.seed_0 = id_branch_a;
					a_parVars.seed_1 = id_branch_b;
				}
				//cout <<"seed0: " << seed0 << " seed1: " << seed1 << " waste: " << waste << " worst: " << worst << endl;
			}
			//cout << endl;
		}
		a_parVars.m_coverSplitArea = volume_minmax_whole_sqrt;
	}
	/*------------  get two approximation with max SAPLA distance  -------------*/
	//rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->copy(*branch_pointer[a_parVars.seed_0].m_rect.rectangle_entry_struct.list_partition_pointer);
	//rectangle_combined.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer->copy(*branch_pointer[a_parVars.seed_1].m_rect.rectangle_entry_struct.list_partition_pointer);
	/*---------------------------------------------------------------------------*/

	/*==========================================================================================================================================================*/

	/*...........................................................................................................................................................*/
#ifdef _DEBUG
	assert(a_parVars.seed_0 != a_parVars.seed_1);
	//assert_branch_rect_in_node(rectangle_combined, branch_pointer);
#endif
	/*...........................................................................................................................................................*/

	return a_parVars.m_coverSplitArea;

}

// Delete a data rectangle from an index structure.   // delete the data, not the Node,internalNode??????
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns 1 if record not found, 0 if success.
// RemoveRect provides for eliminating the root.
RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::RemoveRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root)
{
	ASSERT(a_rect && a_root);
	ASSERT(*a_root);

	ListNode* reInsertList = NULL;

	if (!RemoveRectRec(a_rect, a_id, *a_root, &reInsertList))  // Returns true if record not found, false if success.
	{
		// 1 Found and deleted a data item. 2 Reinsert any branches from eliminated nodes
		// reInsertList links from leaf node to top
		while (reInsertList)
		{
			Node* tempNode = reInsertList->m_node;

			// For the 
			for (int index = 0; index < tempNode->m_count; ++index) {
				// TODO go over this code. should I use (tempNode->m_level - 1)?
				InsertRect(tempNode->m_branch[index], a_root, tempNode->m_level);
			}

			ListNode* remLNode = reInsertList;
			reInsertList = reInsertList->m_next;

			FreeNode(remLNode->m_node);
			FreeListNode(remLNode);
		}

		// Check for redundant root (not leaf, 1 child) and eliminate TODO replace
		// if with while? In case there is a whole branch of redundant roots...
		if ((*a_root)->m_count == 1 && (*a_root)->IsInternalNode())
		{
			Node* tempNode = (*a_root)->m_branch[0].m_child;

			ASSERT(tempNode);
			FreeNode(*a_root);
			*a_root = tempNode;
		}
		return false;
	}
	else
	{
		return true;
	}
}


// Delete a rectangle from non-root part of an index structure.
// Called by RemoveRect.  Descends tree recursively,
// merges branches on the way back up.
// Returns 1 if record not found, 0 if success.
RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::RemoveRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, ListNode** a_listNode) {
	ASSERT(a_rect && a_node && a_listNode);
	ASSERT(a_node->m_level >= 0);

	if (a_node->IsInternalNode()) {  // not a leaf node

		for (int index = 0; index < a_node->m_count; ++index) {

			if (Overlap(a_rect, &(a_node->m_branch[index].m_rect))) {

				//20803 if node a_node->m_branch[index].m_child is leaf node, will disconnect branch below(else part), and retrun false.
				if (!RemoveRectRec(a_rect, a_id, a_node->m_branch[index].m_child, a_listNode)) {

					if (a_node->m_branch[index].m_child->m_count >= TMINNODES)//a_node->m_branch[index].m_child leaf node?
					{
						// child removed, just resize parent rect
						//a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
						NodeCover(a_node->m_branch[index].m_child, a_node->m_branch[index].m_rect);
					}
					else {
						// child removed, not enough entries in node a_node->m_branch[index].m_child, eliminate node
						// a_listNode point from leaf to top
						ReInsert(a_node->m_branch[index].m_child, a_listNode);// Is a_node->m_branch[index].m_child is a leaf node?? a_listNode only store leaf node?
						DisconnectBranch(a_node, index); // Must return after this call as count has changed
					}
					return false;
				}
			}
		}
		return true;
	}
	else { // A leaf node

		for (int index = 0; index < a_node->m_count; ++index)
		{
			if (a_node->m_branch[index].m_data == a_id)
			{
				DisconnectBranch(a_node, index); // Must return after this call as count has changed
				return false;
			}
		}
		return true;
	}
}


// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::Overlap(Rect* a_rectA, Rect* a_rectB) const
{
	ASSERT(a_rectA && a_rectB);

	for (int index = 0; index < NUMDIMS; ++index)
	{
		if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
			a_rectB->m_min[index] > a_rectA->m_max[index])
		{
			return false;
		}
	}
	return true;
}


// Add a node to the reinsertion list.  All its branches will later be reinserted into the index structure.
RTREE_TEMPLATE
void RTREE_PARTITION_QUAL::ReInsert(Node* a_node, ListNode** a_listNode)
{
	ListNode* newListNode;

	newListNode = AllocListNode();
	newListNode->m_node = a_node;
	newListNode->m_next = *a_listNode;
	*a_listNode = newListNode;
}


// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
RTREE_TEMPLATE
bool RTREE_PARTITION_QUAL::Search(Node* a_node, Rect* a_rect, int& a_foundCount, std::function<bool(const DATATYPE&)> callback) const
{
	ASSERT(a_node);
	ASSERT(a_node->m_level >= 0);
	ASSERT(a_rect);

	if (a_node->IsInternalNode())
	{
		// This is an internal node in the tree
		for (int index = 0; index < a_node->m_count; ++index)
		{
			if (Overlap(a_rect, &a_node->m_branch[index].m_rect))
			{
				if (!Search(a_node->m_branch[index].m_child, a_rect, a_foundCount, callback))
				{
					// The callback indicated to stop searching
					return false;
				}
			}
		}
	}
	else
	{
		// This is a leaf node
		for (int index = 0; index < a_node->m_count; ++index)
		{
			if (Overlap(a_rect, &a_node->m_branch[index].m_rect))
			{
				DATATYPE& id = a_node->m_branch[index].m_data;// id is the index of data, not Node
				++a_foundCount;

				if (callback && !callback(id))
				{
					return false; // Don't continue searching
				}
			}
		}
	}

	return true; // Continue searching
}

RTREE_TEMPLATE
template<typename T>//210813
void RTREE_PARTITION_QUAL::update_branch_pointer_vector(T& const node) {
	return;
	/*==================  partition 210618 =========================*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		if (&node) {
			assert_node(node);
			node.branch_ponter_vector.clear();
			node.branch_ponter_vector.shrink_to_fit();

			for (int id_branch = 0; id_branch < node.m_count; id_branch++) {
				node.branch_ponter_vector.emplace_back(&node.m_branch[id_branch]);
			}
		}
		break;
	}
	default: {
		break;
	}
	}
	/*==================================================*/
}

/*=======================================================        Is Function     =====================================================================================*/
RTREE_TEMPLATE
template<typename T>
inline bool RTREE_PARTITION_QUAL::is_rectangle_entry(const T& const rectangle) {
	if (rectangle.rectangle_entry_struct.list_original_pointer) return true;
	else {
		assert(0);
		return false;
	}
}

RTREE_TEMPLATE
template<typename T>
inline bool RTREE_PARTITION_QUAL::is_rectangle_node(const T& const rectangle) {
	if (!rectangle.rectangle_entry_struct.list_original_pointer && rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer) return true;
	else return false;
}

RTREE_TEMPLATE
template<typename T>
inline bool RTREE_PARTITION_QUAL::is_branch_entry(const T& const branch) {
	if (!branch.m_child && branch.m_data != INF) {
		assert(branch.m_rect.rectangle_entry_struct.list_partition_pointer && branch.m_rect.rectangle_entry_struct.original_time_series_vector_pointer);
		return true;
	}
	assert(0);
	return false;
}

RTREE_TEMPLATE
template<typename T>
inline bool RTREE_PARTITION_QUAL::is_branch_node(const T& const branch) {
	if (branch.m_child && branch.m_data == INF) {
		assert(is_rectangle_node(branch.m_rect));
		assert_rectangle_same(branch.m_rect, branch.m_child->rectangle_in_node);
		assert_node(*branch.m_child);

		return true;
	}
	return false;
}

RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::is_branch_pointer_all_entry(const T& const size_array, const Y& const branch_pointer) {

	for (int id_branch = 0; id_branch < size_array; id_branch++) {
		if (!is_branch_entry(branch_pointer[id_branch])) {
			assert(0);
			return false;
		}
	}

	return true;
}

RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::is_branch_pointer_all_node(const T& const size_array, const Y& const branch_pointer) {
	for (int id_branch = 0; id_branch < size_array; id_branch++) {
		if (!is_branch_node(branch_pointer[id_branch])) {
			return false;
		}
	}
	return true;
}
/*====================================================================================================================================================================*/

/*============================================================              210618 Assert            =================================================================*/

//211008 assert option tree
RTREE_TEMPLATE
template<typename T>
inline bool RTREE_PARTITION_QUAL::assert_option_tree(const T& const option_tree_struct) {
	if (option_tree_struct.type_tree != -1 && option_tree_struct.type_representation > 1 && option_tree_struct.type_distance != -1 && option_tree_struct.type_volume != -1) {
		return true;
	}
	assert(0);
	return false;
}


//210725 two vectors have same value and size
RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::assert_same_vector(const T& const vector_0, const Y& const vector_1) {
	assert(vector_0.size() == vector_1.size());
	for (int id_point = 0; id_point < vector_0.size(); id_point++) {
		assert(vector_0[id_point] == vector_1[id_point]);
	}
	return true;
}

// 210729 two SAPLAs have same a&b, right endpoints and rectangle width.
RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::assert_same_SAPLA(const DoublyLinkedList<T>& const doubly_linked_list_0, const DoublyLinkedList<Y>& const doubly_linked_list_1) {
	assert(doubly_linked_list_0.size() == doubly_linked_list_1.size());

	for (int id_segment = 0; id_segment < doubly_linked_list_0.size(); id_segment++) {
		const T& const  segment_0 = doubly_linked_list_0[id_segment];
		const T& const  segment_1 = doubly_linked_list_1[id_segment];
		assert(segment_0.right_endpoint == segment_1.right_endpoint && segment_0.rectangle_width == segment_1.rectangle_width);
		assert(segment_0.apla.a == segment_1.apla.a && segment_0.apla.b == segment_1.apla.b);
	}

	return true;
}

//211227
RTREE_TEMPLATE
template<typename T>
inline bool RTREE_PARTITION_QUAL::assert_rectangle_global_minmax(const T& const rectangle_in_node) {
	assert(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt != -INF);
	assert(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
	assert(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	assert(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_min_pointer);
	assert(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_original_max_pointer);
	assert(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_min_pointer);
	assert(rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.original_time_series_vector_max_pointer);

	return true;
}

//210728
RTREE_TEMPLATE
template<typename T>
bool RTREE_PARTITION_QUAL::assert_rectangle(const T& const rectangle) {
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if (rectangle.rectangle_entry_struct.list_original_pointer != nullptr) {//entry
		//assert(rectangle.original_time_series_vector_pointer->size() == rectangle.reconstruct_time_series_vector_pointer->size());
		assert(rectangle.rectangle_entry_struct.list_original_pointer->size() <= rectangle.rectangle_entry_struct.list_partition_pointer->size());
		if (option_tree_struct.type_tree == 1) {
			assert(rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.empty() && rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.empty() && rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.empty());
		}
		}
	else {// leaf node or internal node
		//assert(!rectangle.original_time_series_vector_pointer && !rectangle.reconstruct_time_series_vector_pointer);
		assert(!rectangle.rectangle_entry_struct.list_original_pointer && !rectangle.rectangle_entry_struct.list_original_pointer);
		assert(rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer && rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
		//assert(!rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.empty() && !rectangle.reconstruct_time_series_pointer_vector.empty() && !rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.empty() && !rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.empty());
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if (option_tree_struct.type_tree == 1) {
		//assert(rectangle.reconstruct_time_series_pointer_vector.size() == rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size());
		assert(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size() == rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size());
		assert(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size() == rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size());
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//assert(rectangle.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.size() == rectangle.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size());
	assert(rectangle.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.size() == rectangle.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.size());
	assert(rectangle.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.size() == rectangle.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.size());

	if (rectangle.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer) {
		APLA::assert_has_same_endpoints(*rectangle.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
		APLA::assert_has_same_endpoints(*rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);

		if (rectangle.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer->size() > 0 && rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer->size() > 0) {
			APLA::assert_has_same_endpoints(*rectangle.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
		}
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	/*for (int id_branch = 0; id_branch < rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size(); id_branch++) {
		assert_same_vector(*rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch], rectangle.original_time_series_vector_vector[id_branch]);
	}*/
	if (option_tree_struct.type_tree == 1) {
		if (!rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.empty()) {
			for (int id_branch_a = 0; id_branch_a < rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size() - 1; ++id_branch_a) {
				for (int id_branch_b = id_branch_a + 1; id_branch_b < rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size(); ++id_branch_b) {
					APLA::assert_has_same_endpoints(*rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a], *rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]);
					//APLA::assert_has_same_endpoints(*rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_a], *rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_b]);
					if (representation_option != 3)
						assert(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_a]->size() == rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_b]->size());
				}
			}
		}
	}

	return true;
}

//210804 Two rectangle has same values. branch.m_rect == node.rectangle_in_node
RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::assert_rectangle_same(const T& const rectangle_0, const Y& const rectangle_1) {
	assert_rectangle(rectangle_0);
	assert_rectangle(rectangle_1);

	//assert(rectangle_0.m_min == rectangle_1.m_min && rectangle_0.m_max == rectangle_1.m_max);
	//assert_same_vector(rectangle_0.m_min, rectangle_1.m_min);
	//assert_same_vector(rectangle_0.m_max, rectangle_1.m_max);

	if (rectangle_0.rectangle_entry_struct.original_time_series_vector_pointer) {// Entry rectangle
		assert(rectangle_0.rectangle_entry_struct.original_time_series_vector_pointer == rectangle_1.rectangle_entry_struct.original_time_series_vector_pointer);
		//assert(rectangle_0.reconstruct_time_series_vector_pointer == rectangle_1.reconstruct_time_series_vector_pointer);
		assert(rectangle_0.rectangle_entry_struct.list_original_pointer == rectangle_1.rectangle_entry_struct.list_original_pointer);
		assert(rectangle_0.rectangle_entry_struct.list_partition_pointer == rectangle_1.rectangle_entry_struct.list_partition_pointer);
	}
	else {// Internal rectangle


		if (rectangle_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer) {// has convex hull volume
			assert(rectangle_0.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment == rectangle_1.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment && rectangle_0.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt == rectangle_1.rectangle_inter_node_struct.minmax_local_struct.volume_minmax_segment_sqrt);
			assert(rectangle_0.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole == rectangle_1.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole && rectangle_0.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt == rectangle_1.rectangle_inter_node_struct.minmax_global_struct.volume_minmax_whole_sqrt);
			assert(rectangle_0.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point == rectangle_1.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point && rectangle_0.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt == rectangle_1.rectangle_inter_node_struct.minmax_point_struct.volume_minmax_point_sqrt);
			assert(rectangle_0.seed_0 == rectangle_1.seed_0 && rectangle_0.seed_1 == rectangle_1.seed_1);

			assert_same_vector(rectangle_0.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector, rectangle_1.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector);
			//(rectangle_0.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector, rectangle_1.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector);
			assert_same_vector(rectangle_0.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector, rectangle_1.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector);
			assert_same_vector(rectangle_0.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector, rectangle_1.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector);
			assert_same_vector(rectangle_0.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector, rectangle_1.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector);
			assert_same_vector(rectangle_0.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector, rectangle_1.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector);

			assert_same_SAPLA(*rectangle_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *rectangle_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);
			assert_same_SAPLA(*rectangle_0.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer, *rectangle_1.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
			if (rectangle_0.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer) {
				assert_same_SAPLA(*rectangle_0.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *rectangle_1.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer);
				assert_same_SAPLA(*rectangle_0.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer, *rectangle_1.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
			}
		}
		if (option_tree_struct.type_tree == 1) {
			// vector of pointer
			assert(rectangle_0.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size() == rectangle_1.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size());
			//assert(rectangle_0.reconstruct_time_series_pointer_vector.size() == rectangle_1.reconstruct_time_series_pointer_vector.size());
			assert(rectangle_0.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size() == rectangle_1.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size());
			assert(rectangle_0.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size() == rectangle_1.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size());
			for (int id_branch = 0; id_branch < rectangle_0.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(); id_branch++) {
				assert(rectangle_0.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch] == rectangle_1.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch]);
				//assert(rectangle_0.reconstruct_time_series_pointer_vector[id_branch] == rectangle_1.reconstruct_time_series_pointer_vector[id_branch]);
				assert(rectangle_0.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch] == rectangle_1.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch]);
				assert(rectangle_0.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch] == rectangle_1.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch]);
			}
		}
	}


	return true;
}

//210728
RTREE_TEMPLATE
template<typename T>
bool RTREE_PARTITION_QUAL::assert_branch(T& const branch) {
	if(option_tree_struct.type_tree == 1){
		assert_rectangle(branch.m_rect);
	}

	//if (branch.node_rect_pointer)
		//assert_rectangle(*branch.node_rect_pointer);

	/*if (branch.branch_to_parent_node_pointer) {
		assert(branch.branch_to_parent_node_pointer->m_branch == &branch);
	}*/

	if (branch.m_child) {// internal node
		assert(branch.m_data == INF);
		assert_node(*branch.m_child);
	}
	else {// entry
		assert(branch.m_data != INF);
	}

	return true;
}

//210729
RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::assert_branch_array(const T& const size_array, Y* branch_pointer) {

	for (int id_branch = 0; id_branch < size_array; id_branch++) {
		assert_branch(branch_pointer[id_branch]);
		if (&branch_pointer[id_branch] != nullptr) {
			//assert(&branch_pointer[id_branch] == &branch_pointer[id_branch].branch_to_parent_node_pointer->m_branch[id_branch]);
		}
	}

	if (size_array > 0) {
		for (int id_branch_a = 0; id_branch_a < size_array - 1; ++id_branch_a) {
			for (int id_branch_b = id_branch_a + 1; id_branch_b < size_array; ++id_branch_b) {
				assert(branch_pointer[id_branch_a].branch_to_parent_node_pointer == branch_pointer[id_branch_b].branch_to_parent_node_pointer);

				if (branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.original_time_series_vector_pointer && branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.original_time_series_vector_pointer) {
					assert(branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.original_time_series_vector_pointer->size() == branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.original_time_series_vector_pointer->size());
					//assert(branch_pointer[id_branch_a].m_rect.reconstruct_time_series_vector_pointer->size() == branch_pointer[id_branch_b].m_rect.reconstruct_time_series_vector_pointer->size());
					if(representation_option != 3)
					assert(branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_original_pointer->size() == branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_original_pointer->size());
					APLA::assert_has_same_endpoints(*branch_pointer[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *branch_pointer[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
				}
				//APLA::assert_has_same_endpoints(*a_node->m_branch[id_branch_a].m_rect.rectangle_entry_struct.list_original_pointer, *a_node->m_branch[id_branch_b].m_rect.rectangle_entry_struct.list_original_pointer);
				//APLA::assert_has_same_endpoints(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_a], *a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_b]);
			}
		}
	}
	return true;
}

//210728
RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::assert_node_branch_pointer(const T& const node, Y* branch_pointer) {
	assert_branch_array(node.m_count, branch_pointer);

	for (int id_branch = 0; id_branch < node.m_count; id_branch++) {
		//assert_branch(branch_pointer[id_branch]);
		if (&branch_pointer[id_branch] != nullptr) {
			assert(&node == branch_pointer[id_branch].branch_to_parent_node_pointer);
			assert(&node.m_branch[id_branch] == &branch_pointer[id_branch]);
		}
	}

	return true;
}

//210730 
RTREE_TEMPLATE
template<typename T, typename Y>
bool RTREE_PARTITION_QUAL::assert_branch_rect_in_node(const T& const rectangle, Y* branch_pointer) {

	assert_rectangle(rectangle);

	if (!branch_pointer || branch_pointer->m_child) {//For leaf node
		assert(branch_pointer->m_data == INF);
		return true;
	}

	if (option_tree_struct.type_tree == 1) {

		assert_branch_array(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(), branch_pointer);

		for (int id_branch = 0; id_branch < rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size(); id_branch++) {

			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			//assert_same_vector(rectangle.original_time_series_vector_vector[id_branch], *branch_pointer[id_branch].m_rect.original_time_series_vector_pointer);
			//assert_same_vector(*rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch], *branch_pointer[id_branch].m_rect.original_time_series_vector_pointer);
			//assert(a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch] == &(branch_pointer[id_branch].m_rect.original_time_series_vector_pointer));

			// pointer same address
			assert(rectangle.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch] == branch_pointer[id_branch].m_rect.rectangle_entry_struct.original_time_series_vector_pointer);
			//assert(rectangle.reconstruct_time_series_pointer_vector[id_branch] == branch_pointer[id_branch].m_rect.reconstruct_time_series_vector_pointer);
			assert(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch] == branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
			assert(rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch] == branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_original_pointer);

			assert_same_SAPLA(*rectangle.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch], *branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
			assert_same_SAPLA(*rectangle.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch], *branch_pointer[id_branch].m_rect.rectangle_entry_struct.list_original_pointer);
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		}
	}

	return true;
}

//210725 Assert node same size and value
RTREE_TEMPLATE
template<typename T>
bool RTREE_PARTITION_QUAL::assert_node_leaf(T& const a_node) {
	if (option_tree_struct.type_tree == 1) {
		assert(a_node.IsLeaf() && a_node.m_branch && a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size() == a_node.m_count);
		//assert(a_node.rectangle_in_node.reconstruct_time_series_pointer_vector.size() == a_node.m_count);

		//assert(a_node.rectangle_in_node.original_time_series_vector_vector.size() == a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size());
		//assert(a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size() == a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector.size());
		//assert(a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector.size() == a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size());
	}
	/*assert(a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_max_vector.size() == a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_point_struct.value_point_min_vector.size());
	assert(a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_min_vector.size() == a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.id_branch_segment_max_vector.size());
	assert(a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_each_segment_max_vector.size() == a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.difference_power_each_segment_max_vector.size());*/

	/*APLA::assert_has_same_endpoints(*a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_max_pointer);
	APLA::assert_has_same_endpoints(*a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer, *a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_max_pointer);
	APLA::assert_has_same_endpoints(*a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_local_struct.list_partition_local_min_pointer, *a_node.rectangle_in_node.rectangle_inter_node_struct.minmax_global_struct.list_partition_global_min_pointer);*/

	assert_node_branch_pointer(a_node, a_node.m_branch);
	assert_branch_rect_in_node(a_node.rectangle_in_node, a_node.m_branch);

	//for (int id_branch = 0; id_branch < a_node->m_count; id_branch++) {
	//	APLA::assert_has_same_endpoints(*a_node->m_branch[0].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
	//	//APLA::assert_has_same_endpoints(*a_node->m_branch[0].m_rect.rectangle_entry_struct.list_original_pointer, *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_original_pointer);
	//	assert(a_node->m_branch[0].branch_to_parent_node_pointer == a_node->m_branch[id_branch].branch_to_parent_node_pointer);

	//	if (&a_node->m_branch[id_branch] != nullptr) {
	//		assert(a_node->m_branch[id_branch].branch_to_parent_node_pointer == a_node);
	//	}

	//	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//	assert_same_vector(a_node->rectangle_in_node.original_time_series_vector_vector[id_branch], a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer);
	//	assert_same_vector(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch], a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer);
	//	assert_same_vector(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch], a_node->rectangle_in_node.original_time_series_vector_vector[id_branch]);
	//	//assert(a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.original_time_series_pointer_vector[id_branch] == &(a_node->m_branch[id_branch].m_rect.original_time_series_vector_pointer));

	//	// pointer same address
	//	assert(a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch] == a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
	//	assert(a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch] == a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_original_pointer);

	//	assert_same_SAPLA(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch], *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_partition_pointer);
	//	assert_same_SAPLA(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch], *a_node->m_branch[id_branch].m_rect.rectangle_entry_struct.list_original_pointer);
	//	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	//
	//	assert_rectangle(a_node->m_branch[id_branch].m_rect);
	//}

	//for (int id_branch_a = 0; id_branch_a < a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size() - 1; ++id_branch_a) {
	//	for (int id_branch_b = id_branch_a + 1; id_branch_b < a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector.size(); ++id_branch_b) {
	//		assert(a_node->m_branch[id_branch_a].m_rect.original_time_series_vector_pointer.size() == a_node->m_branch[id_branch_b].m_rect.original_time_series_vector_pointer.size());
	//		assert(a_node->m_branch[id_branch_a].m_rect.reconstruct_time_series_vector_pointer.size() == a_node->m_branch[id_branch_b].m_rect.reconstruct_time_series_vector_pointer.size());
	//		APLA::assert_has_same_endpoints(*a_node->m_branch[id_branch_a].m_rect.rectangle_entry_struct.list_partition_pointer, *a_node->m_branch[id_branch_b].m_rect.rectangle_entry_struct.list_partition_pointer);
	//		assert(a_node->m_branch[id_branch_a].m_rect.rectangle_entry_struct.list_original_pointer->size() == a_node->m_branch[id_branch_b].m_rect.rectangle_entry_struct.list_original_pointer->size());
	//		assert(a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_a]->size() == a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_b]->size());
	//		//APLA::assert_has_same_endpoints(*a_node->m_branch[id_branch_a].m_rect.rectangle_entry_struct.list_original_pointer, *a_node->m_branch[id_branch_b].m_rect.rectangle_entry_struct.list_original_pointer);
	//		APLA::assert_has_same_endpoints(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a], *a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]);
	//		//APLA::assert_has_same_endpoints(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_a], *a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_b]);
	//	}
	//}

	return true;
}

//************************************
// Method:assert_node()
// Qualifier: from each internal node to each leaf node
// date: 210803
// author:
//************************************
RTREE_TEMPLATE
template<typename T>
bool RTREE_PARTITION_QUAL::assert_node(T& const a_node) {
	/*==================  partition 210618 =========================*/
	switch (representation_option) {
	case 2: // PLA
	case 3://APCA
	case 4://PAA
	case 7://PAALM
	case 9://SAX
	case 5://CHEBY
	case 6: // ICDE07
	case 8: {// Initial 200706
		if (&a_node) {
			assert_node_branch_pointer(a_node, a_node.m_branch);

			assert_branch_rect_in_node(a_node.rectangle_in_node, a_node.m_branch);

			if (option_tree_struct.type_tree == 1) {
				if (a_node.IsInternalNode()) {  // Internal node

					int size_entry_sub = 0;

					for (int id_branch = 0; id_branch < a_node.m_count; ++id_branch) {
						size_entry_sub += a_node.m_branch[id_branch].m_rect.rectangle_inter_node_struct.vector_vector_struct.size_entry();
						assert(a_node.m_branch[id_branch].branch_to_parent_node_pointer == &a_node);
						assert(&a_node.m_branch[id_branch] == a_node.m_branch[id_branch].m_child->node_to_parent_branch_pointer);
						assert_rectangle_same(a_node.m_branch[id_branch].m_child->rectangle_in_node, a_node.m_branch[id_branch].m_rect);
						assert_node(*a_node.m_branch[id_branch].m_child);
					}
					assert(a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.size_entry() == size_entry_sub);
				}
				else { // Leaf node
					assert(a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.size_entry() == a_node.m_count);
					assert_node_leaf(a_node);
				}

				if (a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.size_entry() > 0) {
					for (int id_branch_a = 0; id_branch_a < a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.size_entry() - 1; ++id_branch_a) {
						for (int id_branch_b = id_branch_a + 1; id_branch_b < a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.size_entry(); ++id_branch_b) {

							APLA::assert_has_same_endpoints(*a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_a], *a_node.rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_partition_pointer_vector[id_branch_b]);

							//APLA::assert_has_same_endpoints(*a_node->m_branch[id_branch_a].m_rect.rectangle_entry_struct.list_original_pointer, *a_node->m_branch[id_branch_b].m_rect.rectangle_entry_struct.list_original_pointer);
							//APLA::assert_has_same_endpoints(*a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_a], *a_node->rectangle_in_node.rectangle_inter_node_struct.vector_vector_struct.list_original_pointer_vector[id_branch_b]);
						}
					}
				}
			}
		}
		break;
	}
	default: {
		break;
	}
	}
	/*==================================================*/
	
	return true;
}

/*==========================================================================================================================================================================================================*/

#undef RTREE_TEMPLATE
#undef RTREE_PARTITION_QUAL


#endif //RTREE_H
