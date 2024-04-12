#include <cstdio>
#include <vector>
#include <cassert>

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

#include "heap.h"
#include "flute.h"

inline int getManhattanDist(int x1, int y1, int x2, int y2)
{
  return abs(x1 - x2) + abs(y1 - y2);
}

inline int getHeuristicCost(int cur_x, int cur_y, int dst_x, int dst_y)
{
  return getManhattanDist(cur_x, cur_y, dst_x, dst_y);
}

inline int loc2id(int x, int y, int gridX)
{
  return x + y * gridX;
}

inline void id2loc(int v, int& x, int& y, int gridX)
{
  x = v % gridX;
  y = v / gridX;
}

inline int getNextDst(int v_cur, 
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

 inline void getNeighbors(int* neighbors, int v_cur, int gridX, int gridY)
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

inline void getBias(int numTerm, 
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
				//if(x_term_bias <= x_src && y_term_bias <= y_src)
				if(x_term_bias <= x_src)
					continue;
				if(x_term_bias >= x_dst && y_term_bias >= y_dst)
					continue;
			}
			else if(isRight == true && isUp == false)
			{
				//if(x_term_bias <= x_src && y_term_bias >= y_src)
				if(x_term_bias <= x_src)
					continue;
				if(x_term_bias >= x_dst && y_term_bias <= y_dst)
					continue;
			}
			else if(isRight == false && isUp == true)
			{
				//if(x_term_bias >= x_src && y_term_bias <= y_src)
				if(x_term_bias >= x_src)
					continue;
				if(x_term_bias <= x_dst && y_term_bias >= y_dst)
					continue;
	    }
			else if(isRight == false && isUp == false)
			{
				//if(x_term_bias >= x_src && y_term_bias >= y_src)
				if(x_term_bias >= x_src)
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

inline void ensurePath(int v_src, 
                       int v_dst, 
                       const int* temp_prev,
                             int* prev, 
                       bool* pass)
{
  int v_iter = v_dst;
  while(v_iter != v_src)
  {
    pass[v_iter] = true;
    printf("[ensure] prev of %d is %d\n", v_iter, temp_prev[v_iter]);
    if(prev[v_iter] == -1)
      prev[v_iter] = temp_prev[v_iter];
    v_iter = temp_prev[v_iter];
  }
}

void mazeKernel(const int   x_max,
                const int   x_min,
                const int   y_max,
                const int   y_min,
                const int   numTerm,
                const int*  xpos, /* x coordinates of terms */
                const int*  ypos, /* y coordinates of terms */
                      int*  prev,      
                     bool*  pass)
{
  const int gridX_local = x_max - x_min + 1;
  const int gridY_local = y_max - y_min + 1;

  int*  cost          = new  int[gridX_local * gridY_local];
  int*  temp_prev     = new  int[gridX_local * gridY_local];
  bool* isTerm        = new bool[gridX_local * gridY_local];
  bool* isVisitedTerm = new bool[gridX_local * gridY_local];

  for(int i = 0; i < gridX_local * gridY_local; i++)
  {
    cost[i]          = BIG_INT;
    isTerm[i]        = false;
    isVisitedTerm[i] = false;
  }

  Heap minHeap(gridX_local, gridY_local, cost);

  for(int i = 0; i < numTerm; i++)
  {
    int v_local = loc2id(xpos[i] - x_min, ypos[i] - y_min, gridX_local);
    printf("v_local = %d\n", v_local);
    isTerm[v_local] = true;
  }

  // this should be local_id 
  int root = 0;
  int root_id = loc2id(xpos[root] - x_min, ypos[root] - y_min, gridX_local);
  prev[root_id] = -1;
  isVisitedTerm[root_id] = true;

  int x_src, y_src;
  x_src = xpos[root];
  y_src = ypos[root];
  int v_src = root_id;
  int v_dst = getNextDst(v_src, numTerm, gridX_local, isVisitedTerm, xpos, ypos);
  int x_dst, y_dst;
  id2loc(v_dst, x_dst, y_dst, gridX_local);

  minHeap.insert(root_id, getHeuristicCost(x_src, y_src, x_dst, y_dst));

  int v_cur; // local id 
  while(v_dst != -1)
  {
    bool isSuccess = false;
    
    int x_bias, y_bias;
    getBias(numTerm, gridX_local, 
            x_min, y_min, xpos, ypos, isVisitedTerm, v_src, v_dst, x_bias, y_bias);

    printf("Source : %d Destination : %d\n", v_src, v_dst);
    printf("BiasX : %d BiasY : %d\n", x_bias, y_bias);

    for(int i = 0; i < gridX_local * gridY_local; i++)
      temp_prev[i] = -1;

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
        minHeap.clear();
        for(int i = 0; i < gridX_local * gridY_local; i++)
          cost[i] = BIG_INT;
        minHeap.insert(v_cur, getHeuristicCost(x_src, y_src, x_dst, y_dst));
        ensurePath(v_src, v_dst, temp_prev, prev, pass);
        v_src = v_dst;
        x_src = x_dst;
        y_src = y_dst;
        v_dst = getNextDst(v_src, numTerm, gridX_local, isVisitedTerm, xpos, ypos);
        id2loc(v_dst, x_dst, y_dst, gridX_local);
        break;
      }

      printf("2) v_cur : %d (%d, %d) cost : %d\n", v_cur, x_cur, y_cur, cost[v_cur]);

      int neighbors[4];
      getNeighbors(neighbors, v_cur, gridX_local, gridY_local);

      for(int i = 0; i < 4; i++)
      {
        int v_new = neighbors[i];
        if(v_new == -1)
          continue;

        int x_new, y_new;
        id2loc(v_new, x_new, y_new, gridX_local);
   
        int gCost = (pass[v_new] && pass[v_cur]) ? 0 : minHeap.getCost(v_cur) + getManhattanDist(x_cur, y_cur, x_new, y_new);
        int hCost = getHeuristicCost(x_new, y_new, x_dst, y_dst);
        int bCost = (x_bias == -1 && y_bias == -1) ? 0 : getManhattanDist(x_new, y_new, x_bias, y_bias);

        int newCost = gCost + hCost + bCost;
        // int newCost = gCost + hCost;

        printf("  v_new : %d (%d, %d) oldCost %d newCost : %d\n", v_new, x_new, y_new, minHeap.getCost(v_new), newCost);
        if(i == 0) assert(v_new == v_cur + 1);
        if(i == 1) assert(v_new == v_cur - 1);
        if(i == 2) assert(v_new == v_cur + gridX_local);
        if(i == 3) assert(v_new == v_cur - gridX_local);

        if(newCost < minHeap.getCost(v_new))
        {
          if(minHeap.keyExist(v_new))
            minHeap.setCost(v_new, newCost);
          else
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
  delete isTerm;
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

  d_xpos = ( int*)malloc(numTerm * sizeof(int));
  d_ypos = ( int*)malloc(numTerm * sizeof(int));
  d_prev = ( int*)malloc(gridX * gridY * sizeof(int));
  d_pass = (bool*)malloc(gridX * gridY * sizeof(bool));

  memcpy(d_xpos,  termX.data(), numTerm * sizeof(int));
  memcpy(d_ypos,  termY.data(), numTerm * sizeof(int));
  memcpy(d_pass, h_pass       , gridX * gridY * sizeof(bool));
  memcpy(d_prev, h_prev.data(), gridX * gridY * sizeof(int) );

  mazeKernel(gridX - 1,
             0,
             gridY - 1, 
             0, 
             numTerm, 
             d_xpos, 
             d_ypos, 
             d_prev,
             d_pass);

  memcpy(h_pass       , d_pass, gridX * gridY * sizeof(bool));
  memcpy(h_prev.data(), d_prev, gridX * gridY * sizeof(int) );

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
  printf("flute wl : %d\n", flt::flute_wl(numTerm, termX, termY, 6));

  return 0;
}
