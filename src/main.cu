#include <cstdio>
#include <vector>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include "heap.h"
#include "flute.h"

__device__ inline int getManhattanDist(int x1, int y1, int x2, int y2)
{
  return abs(x1 - x2) + abs(y1 - y2);
}

__device__ inline float getHeuristicCost(int cur_x, int cur_y, int dst_x, int dst_y)
{
  return getManhattanDist(cur_x, cur_y, dst_x, dst_y);
}

__device__ inline int loc2id(int x, int y, int gridX)
{
  return x + y * gridX;
}

__device__ inline int id2loc(int v, int& x, int& y, int gridX)
{
  x = v % gridX;
  y = v / gridX;
}

__device__ inline int getNearestTerm(int v_cur, 
                                     int numTerm,
                                     int gridX_local,
                                     const bool* isVisitedTerm,
                                     const  int* termX, 
                                     const  int* termY)
{
  int v_nearest_term = -1;
  int min_dist = BIG_INT;

  int x_cur, y_cur;
  id2loc(v_cur, x_cur, y_cur, gridX_local);

  for(int i = 0; i < numTerm; i++)
  {
    int x_term = termX[i];
    int y_term = termY[i];
    int v_term = loc2id(x_term, y_term, gridX_local);

    printf("v_cur : %d v_term : %d\n", v_cur, v_term);
    if(v_cur != v_term && isVisitedTerm[v_term] == false)
    {
      int dist = getManhattanDist(x_cur, y_cur, x_term, y_term);
      if(dist < min_dist)
      {
        min_dist = dist;
        v_nearest_term = v_term;
        printf("near term : %d\n", v_nearest_term);
      }
    }
  }

  return v_nearest_term;
}

__device__ inline void getNextSrcDst(int x_min, 
                                     int y_min, 
                                     int numTerm, 
                                     int gridX_local,
                                     const bool* isVisitedTerm,
                                     const  int* termX,
                                     const  int* termY,
                                     int& v_src,
                                     int& v_dst)
{
  int min_dist = BIG_INT;
  v_src = -1;
  v_dst = -1;

  for(int i = 0; i < numTerm; i++)
  {
    // candidate for v_src -> must be visited already
    int x_term1 = termX[i];
    int y_term1 = termY[i];
    int v_term1 = loc2id(x_term1, y_term1, gridX_local);

    if(isVisitedTerm[v_term1] == false)
      continue;

    for(int j = 0; j < numTerm; j++)
    {
      // candidate for v_dst -> must not be visited already
      int x_term2 = termX[j];
      int y_term2 = termY[j];
      int v_term2 = loc2id(x_term2, y_term2, gridX_local);

      if(isVisitedTerm[v_term2] == true || i == j)
        continue;

      int dist = getManhattanDist(x_term1, y_term1, x_term2, y_term2);
      if(dist < min_dist)
      {
        min_dist = dist;
        v_src = v_term1;
        v_dst = v_term2;
      }
    }
  }
}

__device__ inline void getNeighbors(int* neighbors, int v_cur, int gridX, int gridY)
{
  int x, y;
  id2loc(v_cur, x, y, gridX);
  // v_cur must be local id

  neighbors[0] = -1;
  neighbors[1] = -1;
  neighbors[2] = -1;
  neighbors[3] = -1;

  if(x < gridX - 1) neighbors[0] = v_cur + 1; 
  if(x >         0) neighbors[1] = v_cur - 1;  
  if(y < gridY - 1) neighbors[2] = v_cur + gridX; 
  if(y >         0) neighbors[3] = v_cur - gridX; 

  printf("neighbors[0] : %d\n", neighbors[0]);
  printf("neighbors[1] : %d\n", neighbors[1]);
  printf("neighbors[2] : %d\n", neighbors[2]);
  printf("neighbors[3] : %d\n", neighbors[3]);
}

__device__ inline void getBias(int numTerm, 
                               int gridX_local,
                               int x_min,
                               int y_min,
                               const int*  xpos,
                               const int*  ypos,
                               const bool* isVisitedTerm, 
                               int v_src,
                               int v_dst,
                               int& x_bias,
                               int& y_bias)
{
  int sumBiasX = 0;
  int sumBiasY = 0;
  int numBiasToCount = 0;
  
  int x_src, y_src;
  id2loc(v_src, x_src, y_src, gridX_local);

  int x_dst, y_dst;
  id2loc(v_dst, x_dst, y_dst, gridX_local);

  bool isRight = x_dst > x_src ? true : false;
  bool isUp    = y_dst > y_src ? true : false;

  for(int i = 0; i < numTerm; i++)
  {
    int x_term_bias = xpos[i];
    int y_term_bias = ypos[i];
    int v_term_bias = loc2id(x_term_bias, y_term_bias, gridX_local);
    if(isVisitedTerm[v_term_bias] == false && v_term_bias != v_src && v_term_bias != v_dst)
    {
      if(isRight == true && isUp == true)
      {
        if(x_term_bias <= x_src && y_term_bias <= y_src)
          continue;
        if(x_term_bias >= x_dst && y_term_bias >= y_dst)
          continue;
      }
      else if(isRight == true && isUp == false)
      {
        if(x_term_bias <= x_src && y_term_bias >= y_src)
          continue;
        if(x_term_bias >= x_dst && y_term_bias <= y_dst)
          continue;
      }
      else if(isRight == false && isUp == true)
      {
        if(x_term_bias >= x_src && y_term_bias <= y_src)
          continue;
        if(x_term_bias <= x_dst && y_term_bias >= y_dst)
          continue;
      }
      else if(isRight == false && isUp == false)
      {
        if(x_term_bias >= x_src && y_term_bias >= y_src)
          continue;
        if(x_term_bias <= x_dst && y_term_bias <= y_dst)
          continue;
      }
      numBiasToCount++;
      sumBiasX += x_term_bias - x_min;
      sumBiasY += y_term_bias - y_min;
      printf("SumBiasX : %d SumBiasY : %d\n", sumBiasX, sumBiasY);
    }
  }

  if(numBiasToCount == 0)
  {
    x_bias = -1;
    y_bias = -1;
  }
  else
  {
    x_bias = sumBiasX / numBiasToCount;
    y_bias = sumBiasY / numBiasToCount;
  }
}

__device__ inline float getViaCost(int v_temp_prev, int v_prev, int v_cur, int v_new)
{
  if(v_temp_prev == -1 && v_prev == -1)
		return 0.0;
	else if(v_prev - v_cur == v_cur - v_new || v_temp_prev - v_cur == v_cur - v_new)
		return 0.0;
	else
		return 1.0;
}

__device__ inline float getBiasCost(int x_new, int y_new, int x_bias, int y_bias)
{
  if(x_bias == -1 || y_bias == -1)
    return 0.0;
  else
    return getManhattanDist(x_new, y_new, x_bias, y_bias);
}

__device__ inline void ensurePath(int v_src, 
                                  int v_dst, 
                                  const int* temp_prev,
                                        int* prev)
{
  int v_iter = v_dst;
  while(v_iter != v_src)
  {
    printf("[ensure] prev of %d is %d\n", v_iter, temp_prev[v_iter]);
    if(prev[v_iter] == -1)
      prev[v_iter] = temp_prev[v_iter];
    v_iter = temp_prev[v_iter];
  }
}

template<typename T>
__device__ inline void initArray(T* arr, size_t arr_size, T val)
{
	for(size_t i = 0; i < arr_size; i++)
		arr[i] = val;
}

__global__ void mazeKernel(const int   x_max,
                           const int   x_min,
                           const int   y_max,
                           const int   y_min,
                           const int   numTerm,
                           const int*  xpos, /* x coordinates of terms */
                           const int*  ypos, /* y coordinates of terms */
                                 int*  prev)
{
  const int gridX_local = x_max - x_min + 1;
  const int gridY_local = y_max - y_min + 1;

  float*  cost        = new float[gridX_local * gridY_local];
  int*  temp_prev     = new   int[gridX_local * gridY_local];
  bool* isVisitedTerm = new  bool[gridX_local * gridY_local];

	initArray<float>(cost, gridX_local * gridY_local, float(BIG_INT));
	initArray<bool>(isVisitedTerm, gridX_local * gridY_local, false);

  Heap<int, int> minHeap(gridX_local, gridY_local);

	// Assign the most lower-left terminal as a root
	int v_src = BIG_INT;
	int x_src, y_src;
  for(int i = 0; i < numTerm; i++)
	{
	  int x_temp = xpos[i] - x_min;
		int y_temp = ypos[i] - y_min;
    int v_temp = loc2id(x_temp, y_temp, gridX_local);

		if(v_temp < v_src)
		{
			x_src = x_temp;
			y_src = y_temp;
			v_src = v_temp;
		}
	}

  prev[v_src] = -1;
  isVisitedTerm[v_src] = true;

  int v_dst = getNearestTerm(v_src, numTerm, gridX_local, isVisitedTerm, xpos, ypos);
  int x_dst, y_dst;
  id2loc(v_dst, x_dst, y_dst, gridX_local);

	int cost_src = getHeuristicCost(x_src, y_src, x_dst, y_dst);
	cost[v_src] = cost_src;
  minHeap.insert(v_src, cost_src);

  int v_cur; // local id 
  while(v_dst != -1)
  {
    bool isSuccess = false;
    
    int x_bias, y_bias;
    getBias(numTerm, gridX_local, x_min, y_min, xpos, ypos, isVisitedTerm, v_src, v_dst, x_bias, y_bias);

    printf("Source : %d Destination : %d\n", v_src, v_dst);
    printf("BiasX : %d BiasY : %d\n", x_bias, y_bias);

		initArray<int>(temp_prev, gridX_local * gridY_local, -1);

    while(!minHeap.empty())
    {
      v_cur = minHeap.extractMin(); // local_id
      int x_cur, y_cur;
      id2loc(v_cur, x_cur, y_cur, gridX_local);
      printf("1) v_cur : %d (%d, %d)\n", v_cur, x_cur, y_cur);

      if(v_cur == v_dst)
      {
        isSuccess = true;
        printf("Found term : %d\n", v_cur);
        isVisitedTerm[v_cur] = true;
        ensurePath(v_src, v_dst, temp_prev, prev);

        getNextSrcDst(x_min, y_min, numTerm, gridX_local, isVisitedTerm, xpos, ypos, v_src, v_dst);
      
        id2loc(v_src, x_src, y_src, gridX_local);
        id2loc(v_dst, x_dst, y_dst, gridX_local);

        printf("Source : %d (%d, %d) Destination : %d (%d, %d)\n", v_src, x_src, y_src, v_dst, x_dst, y_dst);

        minHeap.clear();

				initArray<float>(cost, gridX_local * gridY_local, float(BIG_INT));

				int cost_src = getHeuristicCost(x_src, y_src, x_dst, y_dst);
				cost[v_src] = cost_src;
        minHeap.insert(v_src, cost_src);
        break;
      }

      printf("2) v_cur : %d (%d, %d) cost : %f\n", v_cur, x_cur, y_cur, cost[v_cur]);

      int neighbors[4];
      getNeighbors(neighbors, v_cur, gridX_local, gridY_local);

      for(int i = 0; i < 4; i++)
      {
        int v_new = neighbors[i];
        if(v_new == -1)
          continue;

        int x_new, y_new;
        id2loc(v_new, x_new, y_new, gridX_local);
   
        float gCost = (prev[v_new] == v_cur || prev[v_cur] == v_new)
                    ? cost[v_cur]
                    : cost[v_cur] + getManhattanDist(x_cur, y_cur, x_new, y_new);
        float hCost = getHeuristicCost(x_new, y_new, x_dst, y_dst);
        float vCost = getViaCost(temp_prev[v_cur], prev[v_cur], v_cur, v_new);
        float bCost = getBiasCost(x_new, y_new, x_bias, y_bias);

        printf("g : %f h : %f v: %f b : %f\n", gCost, hCost, vCost, bCost);

        float newCost = gCost + hCost + bCost * 0.1 + vCost * 0.001;

        printf("  v_new : %d (%d, %d) oldCost %f newCost : %f\n", v_new, x_new, y_new, cost[v_new], newCost);

        if(newCost < cost[v_new])
        {
					cost[v_new] = newCost;
          minHeap.insert(v_new, newCost);
          temp_prev[v_new] = v_cur;
          printf("prev of %d is %d\n", v_new, v_cur);
        }
      }
    }

    if(!isSuccess)
      assert(0);
  }

  delete temp_prev;
  delete cost;
  delete isVisitedTerm;

  printf("Maze Finished!\n");
}


void readInput(const std::string& fileName, 
               int& gridX,
               int& gridY,
               std::vector<int>& termX, 
               std::vector<int>& termY)
{
  std::ifstream fileStream(fileName);
  std::string line;

  std::getline(fileStream, line);
  std::istringstream iss(line);

  iss >> gridX >> gridY;

  int xpos;
  int ypos;

  while(getline(fileStream, line))
  {
    iss.clear();
    iss.str(line);
    iss >> xpos >> ypos;
    termX.push_back(xpos);
    termY.push_back(ypos);
    //printf("x, y : (%d, %d)\n", xpos, ypos);
  }
}

int main(int argc, char* argv[]) 
{
  if(argc != 2)
  {
    printf("No Input!\n");
    assert(0);
  }

  //////////////////////////////////////////////////
  /* Initialization */
  int gridX   = 0;
  int gridY   = 0;
  int numTerm = 0;

  std::vector<int> termX;
  std::vector<int> termY;

  readInput(std::string(argv[1]), gridX, gridY, termX, termY);
  numTerm = termX.size();
  //////////////////////////////////////////////////
  
  int*  d_xpos;
  int*  d_ypos;
  int*  d_prev;
  bool* d_pass;

  bool* h_pass = new bool[gridX * gridY];
  for(int i = 0; i < gridX * gridY; i++)
    h_pass[i] = false;

  std::vector<int> h_prev(gridX * gridY, -1);

  cudaMalloc((void**)&d_xpos, numTerm * sizeof(int));
  cudaMalloc((void**)&d_ypos, numTerm * sizeof(int));
  cudaMalloc((void**)&d_prev, gridX * gridY * sizeof(int));
  cudaMalloc((void**)&d_pass, gridX * gridY * sizeof(bool));

  cudaMemcpy(d_xpos,  termX.data(), numTerm * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ypos,  termY.data(), numTerm * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_pass, h_pass       , gridX * gridY * sizeof(bool), cudaMemcpyHostToDevice);
  cudaMemcpy(d_prev, h_prev.data(), gridX * gridY * sizeof(int) , cudaMemcpyHostToDevice);

  mazeKernel<<<1, 1>>>(gridX - 1,
                       0,
                       gridY - 1, 
                       0, 
                       numTerm, 
                       d_xpos, 
                       d_ypos, 
                       d_prev);

  cudaDeviceSynchronize();

  cudaMemcpy(h_pass       , d_pass, gridX * gridY * sizeof(bool), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_prev.data(), d_prev, gridX * gridY * sizeof(int), cudaMemcpyDeviceToHost);

  std::ofstream output;
  output.open("output.txt");

  output << gridX << " " << gridY << std::endl;

  int wl = 0;
  for(int i = 0; i < h_prev.size(); i++)
  {
    //printf("[output] prev of %d = %d\n", i, h_prev[i]);
    if(h_prev[i] == -1) //|| h_pass[i] == false)
      continue;
    else
    {
      int p1 = h_prev[i];
      int x1 = p1 % gridX;
      int y1 = p1 / gridX;

      int x2 = i % gridX;
      int y2 = i / gridX;
      wl++;

      output << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
    }
  }
  
  printf("      wl : %d\n", wl);
  printf("flute wl : %d\n", flt::flute_wl(numTerm, termX, termY, 3));

	//flt::plottree(flt::flute(termX, termY, 6));

  auto flt_tree = flt::flute(termX, termY, 3);

	for(auto branch : flt_tree.branch)
	{
		int x1 = branch.x;
		int y1 = branch.y;

		int x2 = flt_tree.branch[branch.n].x;
		int y2 = flt_tree.branch[branch.n].y;

    //output << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
	}

  return 0;
}
