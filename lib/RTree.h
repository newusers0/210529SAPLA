﻿#ifndef RTREE_H
#define RTREE_H

// NOTE This file compiles under MSVC 6 SP5 and MSVC .Net 2003 it may not work on other compilers without modification.

// NOTE These next few lines may be win32 specific, you may need to modify them to compile on other platform
//#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <assert.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include<vector>
#include <limits>

//#include "APCA.h"
/*-----210618----------*/
//#include "pch.h"
//#include "lib/doublyLinkedList.h"
//#include "CAPLA.h"
/*---------------*/

using namespace std;
using std::vector;

#define ASSERT assert // RTree uses ASSERT( condition )
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
// RTree.h
//

//#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
//#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, class ELEMTYPEREAL>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, ELEMTYPEREAL>

#define RTREE_DONT_USE_MEMPOOLS // This version does not contain a fixed memory allocator, fill in lines with EXAMPLE to implement one.

//#define RTREE_USE_APCA_SPHERICAL_VOLUME // Better split classification, may be slower on some systems


// Fwd decl
class RTFileStream;  // File I/O helper class, look below for implementation and notes.
//******class RTree;

/// class RTree
/// Implementation of RTree, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use RTree<Object*, float, 3> myTree;
///
/// This modified, templated C++ version by Greg Douglas at Auran (http://www.auran.com)
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
class RTree
{
public:
	struct Node;  // Fwd decl.  Used by other internal structs and iterator
	int NUMDIMS = NULL;
	int TMAXNODES = NULL;
	int TMINNODES = NULL;
	int representation_option = NULL;// 180916 for representation option: PAA, APCA,PLA,CHEBY
	RTree(int NUMDIMS, int &TMAXNODESs);
	RTree(int NUMDIMS, int &TMAXNODESs, const int& representation_option);// 180916 for representation option: PAA, APCA,PLA,CHEBY

public:

	// These constant must be declared after Branch and before Node struct
	// Stuck up here for MSVC 6 compiler.  NSVC .NET 2003 is much happier.
	//enum
	//{
	//	TMAXNODES = TMAXNODES,                         ///< Max elements in node
	//	TMINNODES = TMINNODES,                         ///< Min elements in node
	//};

	//210618 Rtree cover convex hull
	//struct AREA_COEFFICIENT_RTREE_NODE {
	//public:
	//	long double right_endpoint = INF; //The right end point of segment;
	//	long double rectangle_width = INF; //width of rectangle

	//	APLA::APLA_COEFFICIENT apla;//a & b

	//	AREA_COEFFICIENT_RTREE_NODE() {
	//		right_endpoint = INF;
	//		rectangle_width = INF;
	//	}
	//	~AREA_COEFFICIENT_RTREE_NODE() {
	//		//right_endpoint = INF; //The right end point of segment;
	//		//apla.~APLA_COEFFICIENT();
	//	}
	//};

	/// Minimal bounding rectangle (n-dimensional)
	struct Rect
	{
		//ELEMTYPE m_min[NUMDIMS];                      ///< Min dimensions of bounding box 
		ELEMTYPE *m_min = nullptr;
		//ELEMTYPE m_max[NUMDIMS];                      ///< Max dimensions of bounding box 
		ELEMTYPE *m_max = nullptr;

		/*~Rect() {
			delete[] m_min;
			m_min = nullptr;
			delete[] m_max;
			m_max = nullptr;
		}*/

	};

	/// May be data or may be another subtree
	/// The parents level determines this.
	/// If the parents level is 0, then this is data
	// 191111 store the MBR of sub nodes
	struct Branch
	{
		Rect m_rect;                                  ///< Bounds
		Node* m_child = nullptr;                                ///< Child node
		DATATYPE m_data = INF;                              ///< Data Id     :Only for the leaf node??
		/*~Branch() {
			delete[] m_rect.m_min;
			m_rect.m_min = nullptr;
		}*/
	};

	/// Node for each branch level : 191111 node is root, or leaf node, internal node.
	struct Node
	{
		bool IsInternalNode() { return (m_level > 0); } // Not a leaf, but a internal node
		bool IsLeaf() { return (m_level == 0); } // A leaf, contains several entries
		int m_count = NULL;                                  ///< Count  : 191111 the number of sub nodes / branches <= MAX NODES Limits: TMAXNODESs
		int m_level = NULL;                                  ///< Leaf is zero, others positive
		Branch* m_branch = nullptr;                    ///< Branch.size == m_count
	};


public:

	Node * m_root = nullptr;    ///< Root of tree. Initial level is 0
	Node * pointer[10];
	RTree();
	RTree(const RTree& other);
	virtual ~RTree();

	void initialRTree(int NUMDIMS, int &TMAXNODESs, const int& representation_option);// 180917 for other class to use
	/// Insert entry
	/// \param a_min Min of bounding rect
	/// \param a_max Max of bounding rect
	/// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
	void Insert(const ELEMTYPE a_min[], const ELEMTYPE a_max[], const DATATYPE& a_dataId);
	//void Insert(ELEMTYPE a_min[], ELEMTYPE a_max[], const DATATYPE& a_dataId);
	// 191129
	void Insert(const vector<ELEMTYPE>& const a_min, const vector<ELEMTYPE>& const a_max, const DATATYPE& a_dataId);

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
	bool Load(RTFileStream& a_stream);
	/// Save tree contents to file
	bool Save(const char* a_fileName);
	/// Save tree contents to stream
	bool Save(RTFileStream& a_stream);

	void PrintBranchRect(const Branch& m_branch);

	void PrintRect(const Rect& m_rect);

	void printRTree();// 180921

	
	template<typename T>
	void print_tree(T& const node);

	/// Iterator is not remove safe.
	class Iterator
	{
	private:

		enum { MAX_STACK = 1024 }; //201031: 64  Max stack size. Default is 32. Allows almost n^32 where n is number of branches in node

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

			delete[] curBranch.m_rect.m_min;
			curBranch.m_rect.m_min = nullptr;
			delete[] curBranch.m_rect.m_max;
			curBranch.m_rect.m_max = nullptr;

		}

		/// Get the bounds for this node
		void GetBounds(ELEMTYPE a_min[], ELEMTYPE a_max[], const int &NUMDIMS)
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

		friend class RTree; // Allow hiding of non-public functions while allowing manipulation by logical owner
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

		int *m_partition;// = new int;     //  to store the node belong to which group? 0/1????
		int m_total;                       // the number of the rectangle/branches
		int m_minFill;
		int m_count[2];                    //  the number of notes in every group
		Rect m_cover[2];                    // 2 MBRS
		ELEMTYPEREAL m_area[2];            //  the are of each two group

		Branch *m_branchBuf;// = new Branch;
		int m_branchCount;
		Rect m_coverSplit;
		ELEMTYPEREAL m_coverSplitArea;

		~PartitionVars() {
			delete[] m_partition;
			m_partition = nullptr;

			if (m_cover[0].m_min != nullptr) {
				delete[] m_cover[0].m_min;
				m_cover[0].m_min = nullptr;
				delete[] m_cover[0].m_max;
				m_cover[0].m_max = nullptr;
				delete[] m_cover[1].m_min;
				m_cover[1].m_min = nullptr;
				delete[] m_cover[1].m_max;
				m_cover[1].m_max = nullptr;
			}

			if (m_coverSplit.m_min!=nullptr) {
				delete[] m_coverSplit.m_min;
				m_coverSplit.m_min = nullptr;
			}
			if (m_coverSplit.m_max != nullptr) {
				delete[] m_coverSplit.m_max;
				m_coverSplit.m_max = nullptr;
			}

			/*for (int i = 0; i < m_total; i++) {
				if (m_branchBuf[i].m_rect.m_min != nullptr) {
					delete[] m_branchBuf[i].m_rect.m_min;
					m_branchBuf[i].m_rect.m_min = nullptr;
				}
				if (m_branchBuf[i].m_rect.m_max !=nullptr) {
					delete[] m_branchBuf[i].m_rect.m_max;
					m_branchBuf[i].m_rect.m_max = nullptr;
				}

			}*/
			delete[] m_branchBuf;
			m_branchBuf = nullptr;

		}

	};

	Node* AllocNode();
	void AllocBranch(Branch& a_branch);

	void FreeNode(Node* a_node);
	void InitNode(Node* a_node);
	void InitRect(Rect* a_rect);
	bool InsertRectRec(const Branch& a_branch, Node* a_node, Node** a_newNode, int a_level);
	bool InsertRect(const Branch& a_branch, Node** a_root, int a_level);
	//Rect& NodeCover(Node* a_node);
	Rect& NodeCover(Node* a_node, Rect& a_rect);
	bool AddBranch(const Branch* a_branch, Node* a_node, Node** a_newNode);
	void DisconnectBranch(Node* a_node, int a_index);
	int PickBranch(const Rect* a_rect, Node* a_node);
	Rect& CombineRect(const Rect* a_rectA, const Rect* a_rectB);
	Rect CombineRect(const Rect* a_rectA, const Rect* a_rectB, const Rect* a_rect);
	//Rect& CombineRect( Rect* a_rectA,  Rect* a_rectB, Rect& a_rect);
	void assertRect(Rect& a_rect);
	void SplitNode(Node* a_node, const Branch* a_branch, Node** a_newNode);
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
	ELEMTYPEREAL CalcRectVolume(Rect* a_rect);
	void GetBranches(Node* a_node, const Branch* a_branch, PartitionVars* a_parVars);
	void ChoosePartition(PartitionVars* a_parVars, int a_minFill);
	void LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars);
	void InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill);
	void PickSeeds(PartitionVars* a_parVars);
	void Classify(int a_index, int a_group, PartitionVars* a_parVars);
	bool RemoveRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root);
	bool RemoveRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, ListNode** a_listNode);
	ListNode* AllocListNode();
	void FreeListNode(ListNode* a_listNode);
	bool Overlap(Rect* a_rectA, Rect* a_rectB) const;
	void ReInsert(Node* a_node, ListNode** a_listNode);
	bool Search(Node* a_node, Rect* a_rect, int& a_foundCount, std::function<bool(const DATATYPE&)> callback) const;
	void RemoveAllRec(Node* a_node);
	void Reset();
	void CountRec(Node* a_node, int& a_count);

	bool SaveRec(Node* a_node, RTFileStream& a_stream);
	bool LoadRec(Node* a_node, RTFileStream& a_stream);
	void CopyRec(Node* current, Node* other);

	//Node* m_root;                                    ///< Root of tree
	ELEMTYPEREAL m_unitSphereVolume;                 ///< Unit sphere constant for required number of dimensions
};

// Because there is not stream support, this is a quick and dirty file I/O helper.
// Users will likely replace its usage with a Stream implementation from their favorite API.
class RTFileStream
{
	FILE* m_file;

public:


	RTFileStream()
	{
		m_file = NULL;
	}

	~RTFileStream()
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


RTREE_TEMPLATE
RTREE_QUAL::RTree()
{
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

	m_root = AllocNode();
	pointer[0] = m_root;
	m_root->m_level = 0;
	m_unitSphereVolume = (ELEMTYPEREAL)UNIT_SPHERE_VOLUMES[NUMDIMS];
}

RTREE_TEMPLATE
RTREE_QUAL::RTree(const RTree& other) : RTree()
{
	CopyRec(m_root, other.m_root);
}


RTREE_TEMPLATE
RTREE_QUAL::~RTree()
{
#ifdef _DEBUG
	//cout << "RTree clear memory \n";
#endif
	RemoveAll();// 180917 auto free memeory
	Reset(); // Free, or reset node memory
}


RTREE_TEMPLATE
RTREE_QUAL::RTree(int NUMDIMS, int &TMAXNODESs) {//Defualt TMAXNODES=8 TMINNODES = TMAXNODES / 2
	this->NUMDIMS = NUMDIMS;
	TMAXNODES = TMAXNODESs;
	TMINNODES = 2 ;//TMAXNODESs / 2;
	ASSERT(TMAXNODES > TMINNODES);
	ASSERT(TMINNODES > 0);
	//cout << "RTREE: NUMDIMS = " << NUMDIMS << ", TMAXNODES = " << TMAXNODES << ", TMINNODES =" << TMINNODES << endl;
	// Precomputed volumes of the unit spheres for the first few dimensions
	m_root = AllocNode();
	pointer[1] = m_root;
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
// Method:RTree
// Qualifier:
// date:180916
// author: 
//************************************
RTREE_TEMPLATE
RTREE_QUAL::RTree(int NUMDIMS, int &TMAXNODESs,const int& representation_option) {//Defualt TMAXNODES=8 TMINNODES = TMAXNODES / 2
	this->NUMDIMS = NUMDIMS;
	TMAXNODES = TMAXNODESs;
	TMINNODES = 2;// Default TMAXNODESs / 2;
	this->representation_option = representation_option;
	ASSERT(TMAXNODES > TMINNODES);
	ASSERT(TMINNODES > 0);
	//cout << "RTREE: NUMDIMS = " << NUMDIMS << ", TMAXNODES = " << TMAXNODES << ", TMINNODES =" << TMINNODES << endl;
	// Precomputed volumes of the unit spheres for the first few dimensions
	m_root = AllocNode();
	pointer[1] = m_root;
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


RTREE_TEMPLATE
void RTREE_QUAL::initialRTree(int NUMDIMS, int &TMAXNODESs, const int& representation_option) {// 180917 for other class to use
	this->NUMDIMS = NUMDIMS;
	TMAXNODES = TMAXNODESs;
	TMINNODES = 2; //TMAXNODESs / 2;
	this->representation_option = representation_option;
	ASSERT(TMAXNODES > TMINNODES);
	ASSERT(TMINNODES > 0);
	//cout << "RTREE: NUMDIMS = " << NUMDIMS << ", TMAXNODES = " << TMAXNODES << ", TMINNODES =" << TMINNODES << endl;
	// Precomputed volumes of the unit spheres for the first few dimensions
	m_root = AllocNode();
	pointer[1] = m_root;
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
}

RTREE_TEMPLATE
void RTREE_QUAL::Insert(const ELEMTYPE a_min[], const ELEMTYPE a_max[], const DATATYPE& a_dataId)
{
#ifdef _DEBUG
	//cout <<"RTree::Insert() = "<< NUMDIMS << endl;
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
		case 2://APCA
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
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
	branch.m_rect.m_min = new ELEMTYPE[NUMDIMS];
	branch.m_rect.m_max = new ELEMTYPE[NUMDIMS];

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
void RTREE_QUAL::Insert(const vector<ELEMTYPE>& const a_min, const vector<ELEMTYPE>& const a_max, const DATATYPE& const a_dataId) {
#ifdef _DEBUG
	//cout <<"RTree::Insert() = "<< NUMDIMS << endl;
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
		case 2://APCA
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

	/*~~~~~~ Inserted Approximation point ~~~~~~~~*/
	Branch branch;
	branch.m_rect.m_min = new ELEMTYPE[NUMDIMS];
	branch.m_rect.m_max = new ELEMTYPE[NUMDIMS];

	branch.m_data = a_dataId;
	branch.m_child = NULL;

	for (int axis = 0; axis < NUMDIMS; ++axis)
	{
		//cout <<"(a_min : "<< a_min[axis] << ", a_max : " << a_max[axis] <<")  ";
		assert(a_min[axis] != INF && a_max[axis] != INF);
		branch.m_rect.m_min[axis] = a_min[axis];
		branch.m_rect.m_max[axis] = a_max[axis];
	}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	//cout << endl;
	InsertRect(branch, &m_root, 0);
}


RTREE_TEMPLATE
void RTREE_QUAL::Remove(const ELEMTYPE a_min[], const  ELEMTYPE a_max[], const DATATYPE& a_dataId)
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
int RTREE_QUAL::Search(const ELEMTYPE a_min[], const ELEMTYPE a_max[], std::function<bool(const DATATYPE&)> callback) const
{
#ifdef _DEBUG
	for (int index = 0; index < NUMDIMS; ++index)
	{
		ASSERT(a_min[index] <= a_max[index]);
	}
#endif //_DEBUG

	Rect rect;

	/*----------210803--------------------*/
	rect.m_min = new ELEMTYPE[NUMDIMS];
	rect.m_max = new ELEMTYPE[NUMDIMS];
	/*------------------------------------*/

	for (int axis = 0; axis < NUMDIMS; ++axis)
	{
		rect.m_min[axis] = a_min[axis];
		rect.m_max[axis] = a_max[axis];
	}

	// NOTE: May want to return search result another way, perhaps returning the number of found elements here.

	int foundCount = 0;
	Search(m_root, &rect, foundCount, callback);


	/*----------210803--------------------*/
	delete[] rect.m_min;
	rect.m_min = nullptr;
	delete[] rect.m_max;
	rect.m_max = nullptr;
	/*------------------------------------*/

	return foundCount;
}


RTREE_TEMPLATE
void RTREE_QUAL::PrintBranchRect(const Branch& m_branch) {
	for (int i = 0; i < NUMDIMS; i++) {
		cout << "(" << m_branch.m_rect.m_min[i] << "," << m_branch.m_rect.m_max[i] << ") ";
	}
}

RTREE_TEMPLATE
void RTREE_QUAL::PrintRect(const Rect& m_rect) {
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
void RTREE_QUAL::printRTree() {// 180921
	cout << "Root Node : sub node number = " << m_root->m_count << " Root level = : " << m_root->m_level << "\n\nBegin to build a RTree:\n";
	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << Count() << endl;
}

RTREE_TEMPLATE
template<typename T>
void RTREE_QUAL::print_tree(T& const node) {
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
			cout << "  branch id " << id_branch << " id:" << node.m_branch[id_branch].m_data;
		}
		cout << endl;
	}

}


RTREE_TEMPLATE
int RTREE_QUAL::Count()
{
	int count = 0;
	CountRec(m_root, count);
	return count;
}

RTREE_TEMPLATE
void RTREE_QUAL::CountRec(Node* a_node, int& a_count)           //  count the data of the subtree, not Node???
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
void RTREE_QUAL::count_node_size(Node* a_node, T& const count_node_internal, Y& const count_node_leaf) {

	if (a_node->IsInternalNode()) {  // not a leaf node
		assert(a_node->m_level > 0);

		if (a_node->m_level > 1 ) {
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
//void RTREE_QUAL::CountRec(Node* a_node, int& a_count)           //  count the data of the subtree, not Node???
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
bool RTREE_QUAL::Load(const char* a_fileName)
{
	RemoveAll(); // Clear existing tree

	RTFileStream stream;
	if (!stream.OpenRead(a_fileName))
	{
		return false;
	}

	bool result = Load(stream);

	stream.Close();

	return result;
}



RTREE_TEMPLATE
bool RTREE_QUAL::Load(RTFileStream& a_stream)
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
bool RTREE_QUAL::LoadRec(Node* a_node, RTFileStream& a_stream)
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
void RTREE_QUAL::CopyRec(Node* current, Node* other)
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
bool RTREE_QUAL::Save(const char* a_fileName)
{
	RTFileStream stream;
	if (!stream.OpenWrite(a_fileName))
	{
		return false;
	}

	bool result = Save(stream);

	stream.Close();

	return result;
}


RTREE_TEMPLATE
bool RTREE_QUAL::Save(RTFileStream& a_stream)
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
bool RTREE_QUAL::SaveRec(Node* a_node, RTFileStream& a_stream)
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
void RTREE_QUAL::RemoveAll()
{
	// Delete all existing nodes
	Reset();

	m_root = new Node;// AllocNode();
	m_root->m_level = 0;
}


RTREE_TEMPLATE
void RTREE_QUAL::Reset()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	// Delete all existing nodes
	RemoveAllRec(m_root);
#else // RTREE_DONT_USE_MEMPOOLS
	// Just reset memory pools.  We are not using complex types
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::RemoveAllRec(Node* a_node)
{
	ASSERT(a_node);
	ASSERT(a_node->m_level >= 0);

	if (a_node->IsInternalNode()) // This is an internal node in the tree
	{
		for (int index = 0; index < a_node->m_count; ++index) {
			RemoveAllRec(a_node->m_branch[index].m_child);
		}
	}

	FreeNode(a_node);
}


RTREE_TEMPLATE
typename RTREE_QUAL::Node* RTREE_QUAL::AllocNode()
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
void RTREE_QUAL::FreeNode(Node* a_node)
{
	ASSERT(a_node);

#ifdef RTREE_DONT_USE_MEMPOOLS
	for (int index = 0; index < a_node->m_count; ++index)
	{
		if (a_node->m_branch[index].m_rect.m_min != nullptr) {
			delete[] a_node->m_branch[index].m_rect.m_min;
			a_node->m_branch[index].m_rect.m_min = nullptr;
		}

		if (a_node->m_branch[index].m_rect.m_max != nullptr) {
			delete[] a_node->m_branch[index].m_rect.m_max;
			a_node->m_branch[index].m_rect.m_max = nullptr;
		}
	}
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
typename RTREE_QUAL::ListNode* RTREE_QUAL::AllocListNode()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	return new ListNode;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

RTREE_TEMPLATE
void RTREE_QUAL::AllocBranch(Branch& a_branch) {
	a_branch.m_rect.m_min = new ELEMTYPE[NUMDIMS];
	a_branch.m_rect.m_max = new ELEMTYPE[NUMDIMS];
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeListNode(ListNode* a_listNode)
{
#ifdef RTREE_DONT_USE_MEMPOOLS
	delete a_listNode;
#else // RTREE_DONT_USE_MEMPOOLS
	// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::InitNode(Node* a_node)
{
	a_node->m_count = 0;
	a_node->m_level = -1;
	a_node->m_branch = new Branch[TMAXNODES];
}


RTREE_TEMPLATE
void RTREE_QUAL::InitRect(Rect* a_rect)
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
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRectRec(const Branch& a_branch, Node* a_node, Node** a_newNode, int a_level)
{
	ASSERT(a_node && a_newNode);
	ASSERT(a_level >= 0 && a_level <= a_node->m_level);

	// recurse until we reach the correct level for the new record. data records
	// will always be called with a_level == 0 (leaf)
	if (a_node->m_level > a_level)
	{
		// Still above level for insertion, go down tree recursively
		Node* otherNode;// = new Node;

		//memory_account[5] += sizeof(otherNode);
		//otherNode->m_branch = new Branch[TMAXNODES];

		// find the optimal branch for this record
		int index = PickBranch(&a_branch.m_rect, a_node);

		//cout <<"a_node->m_level : " <<a_node->m_level << " Pick Branch :" << index << " a_level : "<< a_level << endl;

		// recursively insert this record into the picked branch
		bool childWasSplit = InsertRectRec(a_branch, a_node->m_branch[index].m_child, &otherNode, a_level);

		if (!childWasSplit)
		{
			// Child was not split. Merge the bounding box of the new record with the
			// existing bounding box
			//cout << "Child was not splited : \n";
			//o a_node->m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));
			CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect), &(a_node->m_branch[index].m_rect));
			return false;
		}
		else
		{
			// Child was split. The old branches are now re-partitioned to two nodes
			// so we have to re-calculate the bounding boxes of each node
			//cout << "Child was splited : \n";
			//a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
			NodeCover(a_node->m_branch[index].m_child, a_node->m_branch[index].m_rect);

			Branch branch;
			branch.m_child = otherNode;
			AllocBranch(branch);

			//branch.m_rect = NodeCover(otherNode);
			NodeCover(otherNode, branch.m_rect);

			// The old node is already a child of a_node. Now add the newly-created
			// node to a_node as well. a_node might be split because of that.
			return AddBranch(&branch, a_node, a_newNode);
		}
	}
	else if (a_node->m_level == a_level)
	{
		// We have reached level for insertion. Add rect, split if necessary
		return AddBranch(&a_branch, a_node, a_newNode);
	}
	else
	{
		// Should never occur
		ASSERT(0);
		return false;
	}
}


// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
//
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRect(const Branch& a_branch, Node** a_root, int a_level)
{
	ASSERT(a_root);
	ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);
#ifdef _DEBUG
	for (int index = 0; index < NUMDIMS; ++index)
	{
		ASSERT(a_branch.m_rect.m_min[index] <= a_branch.m_rect.m_max[index]);
	}
#endif //_DEBUG  

	//cout <<"a_branch id : "<< a_branch.m_data<< ", a_lelve : " << a_level << endl;

	Node* newNode;// = new Node;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!maybe influce result
	//newNode->m_branch = new Branch[TMAXNODES];

	//memory_account[6] += sizeof(newNode->m_branch)*TMAXNODES;

	if (InsertRectRec(a_branch, *a_root, &newNode, a_level))  // Root split
	{
		// Grow tree taller and new root
		Node* newRoot = AllocNode();
		pointer[2] = newRoot;
		newRoot->m_level = (*a_root)->m_level + 1;

		Branch branch;
		AllocBranch(branch);

		// add old root node as a child of the new root
		//branch.m_rect = NodeCover(*a_root);
		NodeCover(*a_root, branch.m_rect);
		branch.m_child = *a_root;
		bool is_split = AddBranch(&branch, newRoot, NULL);

		assert(!is_split);
		// add the split node as a child of the new root
		//branch.m_rect = NodeCover(newNode);

		Branch branch0;
		AllocBranch(branch0);

		NodeCover(newNode, branch0.m_rect);
		branch0.m_child = newNode;
		is_split = AddBranch(&branch0, newRoot, NULL);
		assert(!is_split);

		// set the new root as the root node
		*a_root = newRoot;
		//cout << "true: Root not split" << endl;
		return true;
	}

	//cout << "true: Root was split" << endl;
	return false;
}


// Find the smallest rectangle that includes all rectangles in branches of a node.
//RTREE_TEMPLATE
//typename RTREE_QUAL::Rect& RTREE_QUAL::NodeCover(Node* a_node)
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
RTREE_TEMPLATE
typename RTREE_QUAL::Rect& RTREE_QUAL::NodeCover(Node* a_node, Rect& a_rect) {
	ASSERT(a_node);

	ASSERT(a_rect.m_min != nullptr &&a_rect.m_max != nullptr);

	//assertRect(a_rect);

	/*if (nullptr == a_rect.m_min&&nullptr == a_rect.m_max) {
		a_rect.m_min= new ELEMTYPE[NUMDIMS];
		a_rect.m_max = new ELEMTYPE[NUMDIMS];
	}*/

	//a_rect = a_node->m_branch[0].m_rect;      //: What the index=0 mean?????
	copy_n(a_node->m_branch[0].m_rect.m_min, NUMDIMS, a_rect.m_min);
	copy_n(a_node->m_branch[0].m_rect.m_max, NUMDIMS, a_rect.m_max);
	//	Rect rect = a_node->m_branch[0].m_rect;

	for (int index = 1; index < a_node->m_count; ++index)
	{
		//cout << "CombineRect()" << endl;
		//a_rect = CombineRect(&a_rect, &(a_node->m_branch[index].m_rect));
		CombineRect(&a_rect, &a_node->m_branch[index].m_rect, &a_rect);
	}
	//cout << "rect: "<< rect.m_min[0]  << endl;
	return a_rect;
}


// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch* a_branch, Node* a_node, Node** a_newNode)
{
	ASSERT(a_branch);
	ASSERT(a_node);

	if (a_node->m_count < TMAXNODES)  // Split won't be necessary
	{
		a_node->m_branch[a_node->m_count] = *a_branch;
		++a_node->m_count;

		return false;//0
	}
	else//split node
	{
		ASSERT(a_newNode);
		
		SplitNode(a_node, a_branch, a_newNode);
		return true;//1
	}
}


// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed
RTREE_TEMPLATE
void RTREE_QUAL::DisconnectBranch(Node* a_node, int a_index)
{
	ASSERT(a_node && (a_index >= 0) && (a_index < TMAXNODES));
	ASSERT(a_node->m_count > 0);

	// Remove element by swapping with the last element to prevent gaps in array
	a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];

	--a_node->m_count;
}


// Pick a branch.  Pick the one that will need the smallest increase in area to accommodate the new rectangle.
// This will result in the least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get the best resolution when searching.
RTREE_TEMPLATE
int RTREE_QUAL::PickBranch(const Rect* a_rect, Node* a_node)
{
	ASSERT(a_rect && a_node);

	bool firstTime = true;
	ELEMTYPEREAL increase;
	ELEMTYPEREAL bestIncr = (ELEMTYPEREAL) - 1;
	ELEMTYPEREAL area;
	ELEMTYPEREAL bestArea;
	int best = 0;

	Rect tempRect;
	tempRect.m_min = new ELEMTYPE[NUMDIMS];
	tempRect.m_max = new ELEMTYPE[NUMDIMS];

	for (int index = 0; index < a_node->m_count; ++index)//210618 for each child branch (node/point)
	{
		Rect* curRect = &a_node->m_branch[index].m_rect;// current MBR
		area = CalcRectVolume(curRect);//current volume of MBR

		//tempRect = CombineRect(a_rect, curRect);
		CombineRect(a_rect, curRect, &tempRect);

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

	delete[] tempRect.m_min;
	tempRect.m_min = nullptr;
	delete[] tempRect.m_max;
	tempRect.m_max = nullptr;

	return best;
}


//Combine two rectangles into larger one containing both
RTREE_TEMPLATE
typename RTREE_QUAL::Rect& RTREE_QUAL::CombineRect(const Rect* a_rectA, const Rect* a_rectB)
{
	ASSERT(a_rectA && a_rectB);
	ASSERT(a_rectA->m_min != nullptr &&a_rectB->m_min != nullptr);

	Rect newRect;
	newRect.m_min = new ELEMTYPE[NUMDIMS];
	newRect.m_max = new ELEMTYPE[NUMDIMS];

	for (int index = 0; index < NUMDIMS; ++index) {

		newRect.m_min[index] = (Min)(a_rectA->m_min[index], a_rectB->m_min[index]);
		newRect.m_max[index] = (Max)(a_rectA->m_max[index], a_rectB->m_max[index]);
	}

	return newRect;
}

RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::CombineRect(const Rect* a_rectA, const Rect* a_rectB, const Rect* a_rect) {
	ASSERT(a_rectA && a_rectB);
	ASSERT(a_rectA->m_min != nullptr&&a_rectB->m_min != nullptr&&a_rect->m_min != nullptr);

	for (int index = 0; index < NUMDIMS; ++index) {

		a_rect->m_min[index] = (Min)(a_rectA->m_min[index], a_rectB->m_min[index]);
		a_rect->m_max[index] = (Max)(a_rectA->m_max[index], a_rectB->m_max[index]);
	}

	return *a_rect;
}


//RTREE_TEMPLATE
//typename RTREE_QUAL::Rect& RTREE_QUAL::CombineRect(Rect* a_rectA, Rect* a_rectB, Rect& a_rect) {
//	ASSERT(a_rectA && a_rectB);
//
//	for (int index = 0; index < NUMDIMS; ++index) {
//
//		a_rect.m_min[index] = (Min)(a_rectA->m_min[index], a_rectB->m_min[index]);
//		a_rect.m_max[index] = (Max)(a_rectA->m_max[index], a_rectB->m_max[index]);
//	}
//
//	return a_rect;
//}


RTREE_TEMPLATE
void RTREE_QUAL::assertRect(Rect& a_rect) {
	for (int index = 0; index < NUMDIMS; ++index) {
		ASSERT(a_rect.m_min[index] < a_rect.m_max[index]);
	}
}

// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node* a_node, const Branch* a_branch, Node** a_newNode)
{

	ASSERT(a_node);
	ASSERT(a_branch);

	// Could just use local here, but member or external is faster since it is reused
	PartitionVars localVars;
	localVars.m_partition = new int[TMAXNODES + 1];
	localVars.m_branchBuf = new Branch[TMAXNODES + 1];
	localVars.m_coverSplit.m_min = new ELEMTYPE[NUMDIMS];
	localVars.m_coverSplit.m_max = new ELEMTYPE[NUMDIMS];

	localVars.m_cover[0].m_min = new ELEMTYPE[NUMDIMS];
	localVars.m_cover[0].m_max = new ELEMTYPE[NUMDIMS];
	localVars.m_cover[1].m_min = new ELEMTYPE[NUMDIMS];
	localVars.m_cover[1].m_max = new ELEMTYPE[NUMDIMS];

	//memory_account[7] += sizeof(localVars.m_partition)*(TMAXNODES + 1);
	//memory_account[6] += sizeof(localVars.m_branchBuf)*(TMAXNODES + 1);

	PartitionVars* parVars = &localVars;

	// Load all the branches into a buffer, initialize old node
	GetBranches(a_node, a_branch, parVars);

	// Find partition
	ChoosePartition(parVars, TMINNODES);

	// Create a new node to hold (about) half of the branches
	*a_newNode = AllocNode();
	pointer[3] = *a_newNode;
	(*a_newNode)->m_level = a_node->m_level;

	// Put branches from buffer into 2 nodes according to the chosen partition
	a_node->m_count = 0;
	LoadNodes(a_node, *a_newNode, parVars);


	/*for (int i = 1; i < localVars.m_total; i++) {
	if (localVars.m_branchBuf[i].m_rect.m_min != nullptr) {
	delete[] localVars.m_branchBuf[i].m_rect.m_min;
	localVars.m_branchBuf[i].m_rect.m_min = nullptr;
	}
	delete[] localVars.m_branchBuf[i].m_rect.m_min;
	localVars.m_branchBuf[i].m_rect.m_min = nullptr;
	}*/


	/*delete[] localVars.m_partition;
	localVars.m_partition=nullptr;
	delete[] localVars.m_branchBuf;
	localVars.m_branchBuf = nullptr;*/


	ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);
}


// Calculate the n-dimensional volume of a rectangle      : can get the area of rectangle
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectVolume(Rect* a_rect)
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
ELEMTYPEREAL RTREE_QUAL::APCARectVolume(Rect* a_rect)
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
			//if (a_rect->m_max[index] < a_rect->m_min[index]) {
				//cout << "max: " << a_rect->m_max[index] << ", min: " << a_rect->m_min[index] << endl;
				//ASSERT(0);
			//}
			f_temp_difference = a_rect->m_max[index] - a_rect->m_min[index];
//#ifdef _DEBUG
//			cout<<" temp difference: "<< f_temp_difference <<endl;
//#endif
			//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_min[index] << endl;

			
			//if (f_temp_difference < DBL_TRUE_MIN) continue;//191128 will influence result
			volume *= f_temp_difference;

		}
		else {//even
			//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_min[index] << endl;
			//ASSERT(a_rect->m_max[index] == a_rect->m_min[index]);
			volume *= (a_rect->m_max[index] - a_rect->m_min[index - 2]);
			
		}
		/*-------------------------------------------------------------------------------------------*/
		/*--------------------------------------New version-----------------------------------------*/
		//volume += (a_rect->m_max[index] - a_rect->m_min[index]);
		/*-------------------------------------------------------------------------------------------*/
		//cout << volume << endl;
	}

	
#ifdef _DEBUG
	//if (volume < DBL_TRUE_MIN) {
	//	for (int index = 1; index < NUMDIMS; index++) {
	//		if (index & 1) {//odd
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_min[index] << endl;
	//		}
	//		else {//even
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index - 2 << "] = " << a_rect->m_min[index-2] << endl;
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
ELEMTYPEREAL RTREE_QUAL::apca_rect_accumulation_volume(Rect* a_rect) {
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
			//if (a_rect->m_max[index] < a_rect->m_min[index]) {
				//cout << "max: " << a_rect->m_max[index] << ", min: " << a_rect->m_min[index] << endl;
				//ASSERT(0);
			//}
			f_temp_difference = a_rect->m_max[index] - a_rect->m_min[index];
			//#ifdef _DEBUG
			//			cout<<" temp difference: "<< f_temp_difference <<endl;
			//#endif
						//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_min[index] << endl;


						//if (f_temp_difference < DBL_TRUE_MIN) continue;//191128 will influence result
			volume += f_temp_difference;

		}
		else {//even min id
			//cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_min[index] << endl;
			//ASSERT(a_rect->m_max[index] == a_rect->m_min[index]);
			volume += (a_rect->m_max[index] - a_rect->m_min[index - 2]);

		}
		/*-------------------------------------------------------------------------------------------*/
		/*--------------------------------------New version-----------------------------------------*/
		//volume += (a_rect->m_max[index] - a_rect->m_min[index]);
		/*-------------------------------------------------------------------------------------------*/
		//cout << volume << endl;
	}


#ifdef _DEBUG
	//if (volume < DBL_TRUE_MIN) {
	//	for (int index = 1; index < NUMDIMS; index++) {
	//		if (index & 1) {//odd
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index << "] = " << a_rect->m_min[index] << endl;
	//		}
	//		else {//even
	//			cout << "max[" << index << "] = " << a_rect->m_max[index] << " - min[" << index - 2 << "] = " << a_rect->m_min[index-2] << endl;
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
ELEMTYPEREAL RTREE_QUAL::PLARectVolume(Rect* a_rect)
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
ELEMTYPEREAL RTREE_QUAL::ChebyshevRectVolume(Rect* a_rect)
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
		//cout << "m_max: " << a_rect->m_max[index] << " m_min: " << a_rect->m_min[index] << " max-min:" << a_rect->m_max[index] - a_rect->m_min[index] << endl;
		//ASSERT(a_rect->m_max[index] != 0 && a_rect->m_min[index]!=0);
#ifdef _DEBUG
		ASSERT(a_rect->m_max[index] >= a_rect->m_min[index]);
		if (a_rect->m_max[index] == a_rect->m_min[index]) {
			//cout << "m_max: " << a_rect->m_max[index] << " m_min: " << a_rect->m_min[index] << " max-min:" << a_rect->m_max[index] - a_rect->m_min[index] << endl;
			if (a_rect->m_max[index] != 0) {
				volume0 *= fabs(a_rect->m_max[index]);
			}
			
			//cout << "   volume0: " << volume0 << endl;
			count++;
			continue;
		}
#endif
		//cout << "m_max: " << a_rect->m_max[index] << " m_min: " << a_rect->m_min[index]<<" max-min:"<< a_rect->m_max[index]- a_rect->m_min[index] << endl;
		
		/*--------------------------------------Old version-----------------------------------------*/
		//volume *= (a_rect->m_max[index] - a_rect->m_min[index]) *10.0;
		/*------------------------------------------------------------------------------------------*/
		/*--------------------------------------new version-----------------------------------------*/
		volume += (a_rect->m_max[index] - a_rect->m_min[index]);//191118 change volume calculation of PLA & Cheby from multiple to accummulation
		/*------------------------------------------------------------------------------------------*/
		//cout << "   volume: " << volume << endl;
#ifdef _DEBUG
		ASSERT(volume!=0);
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
ELEMTYPEREAL RTREE_QUAL::apla_rect_volume(Rect* a_rect) {
#ifdef _DEBUG
	ASSERT(a_rect);
#endif

	ELEMTYPEREAL volume = 0;

	for (int index = 0; index < NUMDIMS; ++index){
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

// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
ELEMTYPEREAL RTREE_QUAL::RectSphericalVolume(Rect* a_rect)
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
ELEMTYPEREAL RTREE_QUAL::APCARectSphericalVolume(Rect* a_rect)
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
			//cout << a_rect->m_max[index] << " - " << a_rect->m_min[index] << endl;
		}
		else {
			halfExtent = ((ELEMTYPEREAL)a_rect->m_max[index] - (ELEMTYPEREAL)a_rect->m_min[index - 2]) * 0.5;
			//cout << a_rect->m_max[index] << " - " << a_rect->m_min[index-2] << endl;
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
ELEMTYPEREAL RTREE_QUAL::CalcRectVolume(Rect* a_rect)
{
#ifdef RTREE_USE_APCA_SPHERICAL_VOLUME
	return APCARectSphericalVolume(a_rect); // APCA data point. NUMDIMS < 51 ! Slower but helps certain merge cases 
#else // RTREE_USE_SPHERICAL_VOLUME
#ifdef _DEBUG
	assert(representation_option);
#endif
	//return RectVolume(a_rect);
	//return PLARectVolume(a_rect);
	if (representation_option==1 || representation_option == 2) {//representation_option=1 for 1 PAA, 2 APCA 180916
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
	else {
		return RectVolume(a_rect);
		//assert(0);
	}
	
#endif // RTREE_USE_SPHERICAL_VOLUME  
}


// Load branch buffer with branches from full node plus the extra branch.
RTREE_TEMPLATE
void RTREE_QUAL::GetBranches(Node* a_node, const Branch* a_branch, PartitionVars* a_parVars)
{
	ASSERT(a_node);
	ASSERT(a_branch);

	ASSERT(a_node->m_count == TMAXNODES);

	// Load the branch buffer
	for (int index = 0; index < TMAXNODES; ++index)
	{
		a_parVars->m_branchBuf[index] = a_node->m_branch[index];
		//assertRect(a_parVars->m_branchBuf[index].m_rect);

	}
	a_parVars->m_branchBuf[TMAXNODES] = *a_branch;
	a_parVars->m_branchCount = TMAXNODES + 1;

	// Calculate rect containing all in the set
	//o a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
	
	copy(a_parVars->m_branchBuf[0].m_rect.m_min, a_parVars->m_branchBuf[0].m_rect.m_min + NUMDIMS,a_parVars->m_coverSplit.m_min);
	copy(a_parVars->m_branchBuf[0].m_rect.m_max, a_parVars->m_branchBuf[0].m_rect.m_max + NUMDIMS, a_parVars->m_coverSplit.m_max);


	for (int index = 1; index < TMAXNODES + 1; ++index)
	{
		//o a_parVars->m_coverSplit = CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
		CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect, &a_parVars->m_coverSplit);
	}
	/*cout<<"GetBranches()"<<endl;
	for (int index = 0; index < NUMDIMS; ++index) {
		cout << "(min:" << a_parVars->m_coverSplit.m_min[index] << ", max:" << a_parVars->m_coverSplit.m_max[index] << ") ";
	}
	cout << endl;*/
	a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);
	//cout << "a_parVars->m_coverSplitArea: "<< a_parVars->m_coverSplitArea << endl;
	//cout << "coverSplitArea(): "<< a_parVars->m_coverSplitArea << endl;
}


// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_QUAL::ChoosePartition(PartitionVars* a_parVars, int a_minFill)
{
	//cout << "***ChoosePartition()" << endl;
	ASSERT(a_parVars);

	ELEMTYPEREAL biggestDiff;
	int group, chosen = 0, betterGroup = 0;

	InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
	PickSeeds(a_parVars);

	while (((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
		&& (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill))
		&& (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill)))
	{
		biggestDiff = (ELEMTYPEREAL)-1;
		for (int index = 0; index < a_parVars->m_total; ++index)
		{
			if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index])
			{
				Rect* curRect = &a_parVars->m_branchBuf[index].m_rect;
				Rect rect0;
				Rect rect1;

				rect0.m_min = new ELEMTYPE[NUMDIMS];
				rect0.m_max = new ELEMTYPE[NUMDIMS];
				rect1.m_min = new ELEMTYPE[NUMDIMS];
				rect1.m_max = new ELEMTYPE[NUMDIMS];

				//Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
				//Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);

				CombineRect(curRect, &a_parVars->m_cover[0], &rect0);
				CombineRect(curRect, &a_parVars->m_cover[1], &rect1);

				ELEMTYPEREAL growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];
				ELEMTYPEREAL growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];

				delete[] rect0.m_min;
				rect0.m_min = nullptr;
				delete[] rect0.m_max;
				rect0.m_max = nullptr;
				delete[] rect1.m_min;
				rect1.m_min = nullptr;
				delete[] rect1.m_max;
				rect1.m_max = nullptr;


				ELEMTYPEREAL diff = growth1 - growth0;
				if (diff >= 0)
				{
					group = 0;
				}
				else
				{
					group = 1;
					diff = -diff;
				}

				if (diff > biggestDiff)
				{
					biggestDiff = diff;
					chosen = index;
					betterGroup = group;
				}
				else if ((diff == biggestDiff) && (a_parVars->m_count[group] < a_parVars->m_count[betterGroup]))
				{
					chosen = index;
					betterGroup = group;
				}

			}
		}
		Classify(chosen, betterGroup, a_parVars);
	}

	// If one group too full, put remaining rects in the other
	if ((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
	{
		if (a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill)
		{
			group = 1;
		}
		else
		{
			group = 0;
		}
		for (int index = 0; index < a_parVars->m_total; ++index)
		{
			if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index])
			{
				Classify(index, group, a_parVars);
			}
		}
	}

	ASSERT((a_parVars->m_count[0] + a_parVars->m_count[1]) == a_parVars->m_total);
	ASSERT((a_parVars->m_count[0] >= a_parVars->m_minFill) &&
		(a_parVars->m_count[1] >= a_parVars->m_minFill));
}


// Copy branches from the buffer into two nodes according to the partition.
RTREE_TEMPLATE
void RTREE_QUAL::LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars)
{
	ASSERT(a_nodeA);
	ASSERT(a_nodeB);
	ASSERT(a_parVars);

	for (int index = 0; index < a_parVars->m_total; ++index)
	{
		ASSERT(a_parVars->m_partition[index] == 0 || a_parVars->m_partition[index] == 1);

		int targetNodeIndex = a_parVars->m_partition[index];
		Node* targetNodes[] = { a_nodeA, a_nodeB };

		// It is assured that AddBranch here will not cause a node split. 
		bool nodeWasSplit = AddBranch(&a_parVars->m_branchBuf[index], targetNodes[targetNodeIndex], NULL);
		ASSERT(!nodeWasSplit);
	}
}


// Initialize a PartitionVars structure.
RTREE_TEMPLATE
void RTREE_QUAL::InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill)
{

	ASSERT(a_parVars);

	a_parVars->m_count[0] = a_parVars->m_count[1] = 0;
	a_parVars->m_area[0] = a_parVars->m_area[1] = (ELEMTYPEREAL)0;
	a_parVars->m_total = a_maxRects;
	a_parVars->m_minFill = a_minFill;
	//a_parVars->m_branchBuf = new Branch[TMAXNODES+1];
	//a_parVars->m_partition = new int[TMAXNODES+1];
	for (int index = 0; index < a_maxRects; ++index)
	{
		a_parVars->m_partition[index] = PartitionVars::NOT_TAKEN;
	}
}


RTREE_TEMPLATE
void RTREE_QUAL::PickSeeds(PartitionVars* a_parVars)
{
	//cout << "**** PickSeeds()" << endl;
	int seed0 = 0, seed1 = 0;
	ELEMTYPEREAL worst, waste;
	ELEMTYPEREAL *area = new ELEMTYPEREAL[TMAXNODES + 1];

	for (int index = 0; index < a_parVars->m_total; ++index)
	{
		area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
		//cout << "  area[index]: " << area[index] << endl;
	}

	worst = -a_parVars->m_coverSplitArea - 1;
	//cout << "****** worst: " << worst << endl;
	for (int indexA = 0; indexA < a_parVars->m_total - 1; ++indexA)    // Bubble sort????????
	{
		for (int indexB = indexA + 1; indexB < a_parVars->m_total; ++indexB)
		{
			Rect oneRect;
			oneRect.m_min = new ELEMTYPE[NUMDIMS];
			oneRect.m_max = new ELEMTYPE[NUMDIMS];

			//Rect oneRect = CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect);
			CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect, &oneRect);
			waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];

#ifdef _DEBUG
			//cout << "combined rectangle: \n";
			//for (int segment_id = 0; segment_id < NUMDIMS; segment_id++) {
				//cout << oneRect.m_min[segment_id]<<" , "<< oneRect.m_max[segment_id] <<endl;
			//}
#endif
			delete[] oneRect.m_min;
			oneRect.m_min = nullptr;
			delete[] oneRect.m_max;
			oneRect.m_max = nullptr;
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

	delete[] area;
	area = nullptr;

	//cout << "******** seed0: "   << seed0 << endl;
	Classify(seed0, 0, a_parVars);
	//cout << "******** seed1:  " << seed1 << endl;
	Classify(seed1, 1, a_parVars);
}


// Put a branch in one of the groups.
RTREE_TEMPLATE
void RTREE_QUAL::Classify(int a_index, int a_group, PartitionVars* a_parVars)
{
	//cout << "*********Classify()" << endl;

	//cout << "seed: " << a_index << " group:"<< a_group<<"a_parVars->m_partition[a_index]: "<<a_parVars->m_partition[a_index] << endl;

	ASSERT(a_parVars);

	ASSERT(PartitionVars::NOT_TAKEN == a_parVars->m_partition[a_index]);
	a_parVars->m_partition[a_index] = a_group;

	// Calculate combined rect
	if (a_parVars->m_count[a_group] == 0)
	{
		//a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
		copy_n(a_parVars->m_branchBuf[a_index].m_rect.m_min, NUMDIMS, a_parVars->m_cover[a_group].m_min);
		copy_n(a_parVars->m_branchBuf[a_index].m_rect.m_max, NUMDIMS, a_parVars->m_cover[a_group].m_max);
	}
	else
	{
		//a_parVars->m_cover[a_group] = CombineRect(&a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group]);
		CombineRect(&a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group], &a_parVars->m_cover[a_group]);
	}

	// Calculate volume of combined rect
	a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);

	++a_parVars->m_count[a_group];
}


// Delete a data rectangle from an index structure.   // delete the data, not the Node,internalNode??????
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns 1 if record not found, 0 if success.
// RemoveRect provides for eliminating the root.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRect(Rect* a_rect, const DATATYPE& a_id, Node** a_root)
{
	ASSERT(a_rect && a_root);
	ASSERT(*a_root);

	ListNode* reInsertList = NULL;

	if (!RemoveRectRec(a_rect, a_id, *a_root, &reInsertList))  // Returns true if record not found, false if success.
	{
		// Found and deleted a data item
		// Reinsert any branches from eliminated nodes
		while (reInsertList)
		{
			Node* tempNode = reInsertList->m_node;

			for (int index = 0; index < tempNode->m_count; ++index)
			{
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
bool RTREE_QUAL::RemoveRectRec(Rect* a_rect, const DATATYPE& a_id, Node* a_node, ListNode** a_listNode)
{
	ASSERT(a_rect && a_node && a_listNode);
	ASSERT(a_node->m_level >= 0);

	if (a_node->IsInternalNode())  // not a leaf node
	{
		for (int index = 0; index < a_node->m_count; ++index)
		{
			if (Overlap(a_rect, &(a_node->m_branch[index].m_rect)))
			{
				if (!RemoveRectRec(a_rect, a_id, a_node->m_branch[index].m_child, a_listNode))
				{
					if (a_node->m_branch[index].m_child->m_count >= TMINNODES)
					{
						// child removed, just resize parent rect
						//a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
						NodeCover(a_node->m_branch[index].m_child, a_node->m_branch[index].m_rect);
					}
					else
					{
						// child removed, not enough entries in node, eliminate node
						ReInsert(a_node->m_branch[index].m_child, a_listNode);
						DisconnectBranch(a_node, index); // Must return after this call as count has changed
					}
					return false;
				}
			}
		}
		return true;
	}
	else // A leaf node
	{
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
bool RTREE_QUAL::Overlap(Rect* a_rectA, Rect* a_rectB) const
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
void RTREE_QUAL::ReInsert(Node* a_node, ListNode** a_listNode)
{
	ListNode* newListNode;

	newListNode = AllocListNode();
	newListNode->m_node = a_node;
	newListNode->m_next = *a_listNode;
	*a_listNode = newListNode;
}


// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Search(Node* a_node, Rect* a_rect, int& a_foundCount, std::function<bool(const DATATYPE&)> callback) const
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


#undef RTREE_TEMPLATE
#undef RTREE_QUAL

#endif //RTREE_H

