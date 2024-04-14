#include <cassert>
#include <cmath> // for infinity
#include <stdio.h>

#define BIG_INT 2147483647

template<typename Key, typename Value>
struct Node 
{
	__host__ __device__ Node() {}

	__host__ __device__ Node(Key _key, Value _val)
		: key (_key), val (_val) {}

	__device__ __host__
  Node<Key, Value>& operator=(Node<Key, Value> const &other)
	{
    key = other.key;
    val = other.val;
    return *this;
  }

	Key   key;
	Value val;
};

template<typename Key, typename Value>
struct Heap /* minHeap */
{
  using Node_t = Node<Key, Value>;

	public: 

    __host__ __device__ Heap(int gridX, int gridY)
  	{
  		heapSize_ = 0;
			array_  = new Node_t[gridX * gridY];
  	}

		__host__ __device__ ~Heap()
		{
			//printf("destructor called\n");
			delete array_;
		}

  	__host__ __device__ Key extractMin()      
  	{
      if(heapSize_ < 1)
  			assert(0);
  
  		Node_t minNode = array_[0];
  		array_[0] = array_[heapSize_ - 1];
  		heapSize_--;
  		heapify(0);
  
  		return minNode.key;
  	}
  
  	__host__ __device__ void print() const
  	{
  		for(size_t i = 0; i < heapSize_; i++)
  		  printf("Key : %d Value : %d\n", array_[i].key, array_[i].val);
  	}
  
  	__host__ __device__ void insert(Node_t node)
  	{
  	  heapSize_++;
  		array_[heapSize_ - 1] = Node_t(-1, BIG_INT);
      decreaseKey(heapSize_ - 1, node);
  	}

  	__host__ __device__ void insert(Key key, Value val)
  	{
			Node_t newNode(key, val);
			insert(newNode);
  	}

    __host__ __device__ bool empty() const
		{
			return heapSize_ == 0;
		}

		__host__ __device__ void clear() 
		{
			heapSize_ = 0;
		}

	private:

 	  __host__ __device__ int  getParent(int  childIdx) const { return  childIdx / 2;     }
  	__host__ __device__ int  getLeft  (int parentIdx) const { return parentIdx * 2 + 1; }
  	__host__ __device__ int  getRight (int parentIdx) const { return parentIdx * 2 + 2; }
    __host__ __device__ void swap(int idx1, int idx2) 
  	{
			Node_t temp  = array_[idx1];
  		array_[idx1] = array_[idx2];
  		array_[idx2] = temp;
  	}
  
    __host__ __device__ void heapify(int idx)
  	{
  		int left  = getLeft(idx);
  		int right = getRight(idx); 
      int minIdx;
  
  		if(left < heapSize_ && array_[left].val < array_[idx].val)
        minIdx = left;
  		else
  			minIdx = idx;
  
  		if(right < heapSize_ && array_[right].val < array_[minIdx].val)
  			minIdx = right;
  
  		if(minIdx != idx)
  		{
  			swap(minIdx, idx);
  			heapify(minIdx);
  		}
  	}

  	__host__ __device__ void decreaseKey(int idx, Node_t node)
  	{
  		if(array_[idx].val < node.val)
  			assert(0);
  
  		array_[idx] = node;
  
  		int _idx = idx;
  		while(_idx > 0 && array_[getParent(_idx)].val > array_[_idx].val)
  		{
        swap(_idx, getParent(_idx));
  			_idx = getParent(_idx);
  		}
  	}

    Node_t* array_;
	  int heapSize_;
};
