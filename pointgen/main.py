import os
import argparse
import pointGeneration as pg

parser = argparse.ArgumentParser()
#parser.add_argument("fileName", help="file name of point set")
parser.add_argument("-n", "--pointNum", type=int, default=100, help="number of points")
parser.add_argument("-g", "--gridSize", type=int, default=100, help="range of the points 0 ~ gridSize")
parser.add_argument("-s", "--sigma", type=int, default=15, help="std deviation of normal distribution")
args = parser.parse_args()
    
pointFilePath = "./"

fileNumList = []
for root, dirs, files in os.walk(pointFilePath):
    fileNumList = []
    for name in files:
        if(name[:5] == "point"):
            tokens = name.split("_")
            if(tokens[0][5:] == str(args.pointNum)):
                if(tokens[1].split("x")[0] == str(args.gridSize)):
                    if(tokens[2] == str(args.sigma)):
                        fileNum = int(tokens[3].split(".")[0])
                        fileNumList.append(fileNum)

if(len(fileNumList) == 0):
    fileNum = 1
else:
    fileNum = max(fileNumList) + 1
fileName = "point" + str(args.pointNum) + "_" + str(args.gridSize) + "x" + str(args.gridSize) + "_" + str(args.sigma) + "_" + str(fileNum) + ".pt"
fileName = os.path.join(pointFilePath, fileName)

# random number of points [low, high]
low = high = args.pointNum
# number of point set
num = 1
# std deviation of normal distribution
sigma = args.sigma
# grid size
width = args.gridSize

pointset = pg.point_generate(low, high, num, sigma, width)

wfile = open(fileName, 'w')
for points in pointset:
    num_str = "%d\n" % (len(points)) 
    wfile.write(num_str)
    for point in points:
        point_str = "%d %d\n" % (point[0], point[1])
        wfile.write(point_str)
wfile.close()
