This repo have two files containing submission to PACE Challenge in 2018 in Track-A and Track-C.

trackA.cpp contains a code which accepts an undirected, edge-weighted graph and set of terminals as input and returns optimums Steiner Tree. It is implementation of an exact algorithm by Erickson, Monma, Veinott. One can run the code using following commands.

g++ -std=c++11 trackA.cpp -o trackA

./trackA < ../Path-to-instances/


trackC.py contains a code which accepts the same input and solves Steiner Tree problem using heuristics methods. The code does not terminate by itself and needs to be killed by user. To run code for 10 sec, use following command. 

timeout 10 python trackC.py < ../Path-to-instances


PACE Challenge: https://pacechallenge.wordpress.com/pace-2018/
