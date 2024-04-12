import numpy as np
import matplotlib.pyplot as plt

with open("output.txt", 'r') as INPUT_OPEN:
    INPUT = INPUT_OPEN.readlines()

pointX = []
pointY = []
edgeList = []

dimX = int(INPUT[0].split(" ")[0])
dimY = int(INPUT[0].split(" ")[1])

for line in INPUT[1:]:
    point = line.split(" ")
    x1 = int(point[0])
    y1 = int(point[1])
    x2 = int(point[2])
    y2 = int(point[3])

    pointX.append(x1)
    pointX.append(x2)
    pointY.append(y1)
    pointY.append(y2)

    edge = [x1, y1, x2, y2]
    edgeList.append(edge)

print(dimX)
print(dimY)

plt.xticks(range(dimX + 1))
plt.yticks(range(dimY + 1))
plt.xlim(-1, dimX + 1)
plt.ylim(-1, dimY + 1)
plt.grid(color = 'black', linestyle = '-', linewidth = 0.2)
plt.scatter(pointX, pointY, color="black", s=100)

for edge in edgeList:
    x1 = edge[0]
    y1 = edge[1]

    x2 = edge[2]
    y2 = edge[3]

    print(f"x1 : {x1} y1 : {y1} x2 : {x2} y2 : {y2}")

    plt.arrow(x1, y1, x2 - x1, y2 - y1)

plt.show()
