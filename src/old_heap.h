#include <cassert>
#include <cmath> // for infinity
#include <stdio.h>

#define BIG_INT 2147483647

struct Heap /* minHeap */
{
	public: 

    __host__ __device__ Heap(int gridX, int gridY, float* cost)
  	{
			gridX_        = gridX;
			gridY_        = gridY;
  		heapSize_     = 0;
			key2cost_     = cost;
      key2heapIdx_  = new  int[gridX * gridY];
			heapIdx2key_  = new  int[gridX * gridY];
			keyExist_     = new bool[gridX * gridY];

			for(int i = 0; i < gridX_ * gridY_; i++)
			  keyExist_[i] = false;
  	}

		__host__ __device__ ~Heap()
		{
			//printf("destructor called\n");
			delete keyExist_;
			delete key2heapIdx_;
			delete heapIdx2key_;
		}

  	__host__ __device__ int extractMin()      
  	{
      if(heapSize_ < 1)
  			assert(0);
  
  		int minNode = heapIdx2key_[0];
  		heapIdx2key_[0] = heapIdx2key_[heapSize_ - 1];
  		heapSize_--;
  		heapify(0);
  
      keyExist_[minNode] = false;
  		return minNode;
  	}
  
  	__host__ __device__ void print() const
  	{
  		for(size_t i = 0; i < heapSize_; i++)
  		  printf("[%d] Key : %d Value : %f\n", int(i), heapIdx2key_[i], key2cost_[heapIdx2key_[i]]);
  	}
  
  	__host__ __device__ void insert(int key, float val)
  	{
  	  heapSize_++;
			key2heapIdx_[key] = heapSize_ - 1;
  		heapIdx2key_[heapSize_ - 1] = key;
			key2cost_[key] = float(BIG_INT);
      decreaseVal(heapSize_ - 1, val);
			keyExist_[key] = true;
  	}

    __host__ __device__ bool keyExist(int key) const
		{
			return keyExist_[key];
		}

		__host__ __device__ void clear() 
		{
			heapSize_ = 0;
			for(int i = 0; i < gridX_ * gridY_; i++)
			  keyExist_[i] = false;
		}

		__host__ __device__ bool empty() const
		{
			return heapSize_ == 0;
		}

		__host__ __device__ void setCost(int key, float newCost)
		{
			assert(key2heapIdx_[key] < heapSize_);
			if(key2cost_[key] < newCost) 
				increaseVal(key2heapIdx_[key], newCost);
			else                         
				decreaseVal(key2heapIdx_[key], newCost);
		}
	
		__host__ __device__ float getCost(int key) const
		{
			// assert(key2heapIdx_[key] < heapSize_);
		  return key2cost_[key];
		}

	private:

  	__host__ __device__ void decreaseVal(int idx, float val)
  	{
  		key2cost_[heapIdx2key_[idx]] = val;
  
  		int _idx = idx;
  		while(_idx > 0 && !compare(getParent(_idx), _idx))
  		{
        swap(_idx, getParent(_idx));
  			_idx = getParent(_idx);
  		}
  	}

  	__host__ __device__ void increaseVal(int idx, float val)
  	{
  		key2cost_[heapIdx2key_[idx]] = val;
			heapify(idx);
  	}

 	  __host__ __device__ int  getParent(int  childIdx) const { return  childIdx / 2;     }
  	__host__ __device__ int  getLeft  (int parentIdx) const { return parentIdx * 2 + 1; }
  	__host__ __device__ int  getRight (int parentIdx) const { return parentIdx * 2 + 2; }
    __host__ __device__ void swap(int idx1, int idx2) 
  	{
			int tempKey = heapIdx2key_[idx1];
  		heapIdx2key_[idx1] = heapIdx2key_[idx2];
  		heapIdx2key_[idx2] = tempKey;

			int key1 = heapIdx2key_[idx1];
			int key2 = heapIdx2key_[idx2];
			key2heapIdx_[key1] = idx1;
			key2heapIdx_[key2] = idx2;
  	}
  
    __host__ __device__ void heapify(int idx)
  	{
  		int left  = getLeft(idx);
  		int right = getRight(idx); 
      int minIdx;
  
  		if(left < heapSize_ && compare(left, idx))
        minIdx = left;
  		else
  			minIdx = idx;
  
  		if(right < heapSize_ && compare(right, minIdx))
  			minIdx = right;
  
  		if(minIdx != idx)
  		{
  			swap(minIdx, idx);
  			heapify(minIdx);
  		}
  	}
 
		__host__ __device__ bool compare(int heapIdx1, int heapIdx2) const
		{
      return key2cost_[heapIdx2key_[heapIdx1]] < key2cost_[heapIdx2key_[heapIdx2]];
		}

		bool* keyExist_;
		int*  heapIdx2key_;
		int*  key2heapIdx_;
		float*  key2cost_;
	  int   heapSize_;

		int gridX_;
		int gridY_;
};