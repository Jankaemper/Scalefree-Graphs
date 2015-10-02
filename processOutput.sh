#! /bin/bash
dot -Tps Output/graph_visu.dot -o Output/graph_visu.ps
evince Output/graph_visu.ps&
cat Output/distance.dat
