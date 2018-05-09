#include <bits/stdc++.h>

using namespace std;

void SO(){cout << endl;}
template<typename S, typename ... Strings>
void SO(S x, const Strings&... rest){cout<<x<<' ';SO(rest...);}

#define pb(x) push_back(x)
#define mp(x, y) make_pair(x, y)

const int AM_MAX = 10001, MAX_NODES = 10000001;

int GRAPH_TYPE; //0 denotes AdjacencyMatrix. 1 denotes AdjacencyHash

int adjacencyMatrix[AM_MAX][AM_MAX];
unordered_map<int, int> adjacencyHash[MAX_NODES];
set<int> adjacency[MAX_NODES]; //to-vertex

int cnodes, cedges, cterminals, totalNodes;

vector<int>terminals;
//int terminalNeighbour[MAX_NODES];
int marked[MAX_NODES]; //marked[i] = i if steiner. = posn(i) in vector<int>terminals if terminal.. -1 otherwise
int degree[MAX_NODES]; //degree -1 means deleted vertex

const int MAX_DEGREE_SET = 5;
set<int> degreeVertex[MAX_DEGREE_SET];

int auxillary_node = 0; //The last node is created as the auxillary_NODE. This is for the "shortest path" in EMV.

int parent[MAX_NODES]; //Parent of the node X when Y is the source : parent[X][Y] 
int shortest_paths[MAX_NODES][2]; //Updated everytime shortestPaths function runs. Second index is the dijkstra-ID for which dijkstra was run on
int dijktraID = 1;

int EMV_dp[MAX_NODES][4];//DP for Dreyfus Wagner Acc to the thesis... 0 -> Fv, 1 -> Bv[0], 2 -> Bv[1], 3-> Gv
set<pair<int, int> > solution, toPrint;

map< pair<int, int>, vector<pair<int, int> > >alias;

//****************************************** Graph **************************************************

bool valid(int a)
{
	if(a > 0 and a < 1e9+7)
		return true;
	return false;
}

int graphGet(int a, int b)
{
	if(GRAPH_TYPE == 0)
		return adjacencyMatrix[a][b];
	if(adjacencyHash[a].find(b) == adjacencyHash[a].end())
		return 0;
	return adjacencyHash[a][b];
}

void graphSet(int a, int b, int c)
{
	if(GRAPH_TYPE == 0)
		adjacencyMatrix[a][b] = c;
	else
		adjacencyHash[a][b] = c;
}

void getAdjacency()
{
	for(int i = 1; i <= totalNodes; i++)
	{
		adjacency[auxillary_node].insert(i);
		graphSet(i, auxillary_node, 1e9+7);
		graphSet(auxillary_node, i, 1e9+7);
		for(int j = 1; j <= totalNodes; j++)
		{
			if(graphGet(i, j) > 0)
			{
				adjacency[i].insert(j);
			}
		}
	}
	totalNodes = cnodes;
}

void getInput()
{
	string temp;
	int a, b, c;
	cin>>temp>>temp>>temp;
	cin>>cnodes>>temp>>cedges;
	
	if(cnodes<AM_MAX)
	{
		GRAPH_TYPE = 0;
	}
	else
	{
		GRAPH_TYPE = 1;
	}
		
	//Undirected graph Input
	for(int i = 0; i < cedges; i++)
	{
		cin>>temp>>a>>b>>c;
		graphSet(a, b, c);
		graphSet(b, a, c);
		degree[a]++;
		degree[b]++;
		marked[a] = -1;
		marked[b] = -1;
	}
	cin>>temp>>temp>>temp>>temp>>cterminals;
	//Terminals are stored in "terminals". And each terminal is marked as i. So when u want the index of the terminal, query "marked"
	//Note: Even a non-terminal can have "marked[a] = i
	for(int i = 0; i < cterminals; i++)
	{
		cin>>temp>>a;
		terminals.push_back(a);
		marked[a] = i;
	}
	
	totalNodes = cnodes;
	
	for(int i = 1; i <= totalNodes; i++)
	{
		if(degree[i] < MAX_DEGREE_SET)
			degreeVertex[degree[i]].insert(i);
	}
	
	getAdjacency();
}

void printAlias(int a, int b)
{
	if(alias.find(mp(min(a, b), max(a, b))) == alias.end())
		toPrint.insert(mp(a, b));
	else
	{
		for(auto i : alias[mp(min(a, b), max(a, b))])
			printAlias(i.first, i.second);
	}
}

void kernalizationDeg1();
void kernalizationDeg2(int x = -1)
{
	for(auto i : degreeVertex[2])
	{
		if(marked[i] == -1 or x != -1)
		{
			if(x != -1)
				i = x;
			
			set<int>::iterator it = adjacency[i].begin();
			int u = *it;
			++it;
			int v = *it;
			if(graphGet(u, v) > 0 and graphGet(u, v) < graphGet(u, i) + graphGet(i, v))
			{
				if(x != -1)
					return;
				else
					continue;
			}
			else
			{
				if(graphGet(u, v) > 0)
				{
					if(degree[u] < MAX_DEGREE_SET)
						degreeVertex[degree[u]].erase(u);
					degree[u]--;
					if(degree[u] < MAX_DEGREE_SET)
						degreeVertex[degree[u]].insert(u);
					
					if(degree[v] < MAX_DEGREE_SET)
						degreeVertex[degree[v]].erase(v);
					degree[v]--;
					if(degree[v] < MAX_DEGREE_SET)
						degreeVertex[degree[v]].insert(v);
				}
				
				alias[mp(min(u, v), max(u, v))].push_back(mp(min(u, i), max(u, i)));
				alias[mp(min(u, v), max(u, v))].push_back(mp(min(v, i), max(v, i)));
				
				graphSet(u, v, graphGet(u, i) + graphGet(i, v));
				graphSet(v, u, graphGet(u, i) + graphGet(i, v));
				graphSet(u, i, 0);
				graphSet(i, u, 0);
				graphSet(v, i, 0);
				graphSet(i, v, 0);
				
				adjacency[v].erase(i);
				adjacency[u].erase(i);
				adjacency[v].insert(u);
				adjacency[u].insert(v);
				adjacency[i].erase(u);
				adjacency[i].erase(v);
				
				marked[i] = -2;
				
				if(degree[v] == 1 or degree[u] == 1)
					kernalizationDeg1();
					
				if(x != -1)
					return;
			}
		}
	}
}

void kernalizationDeg1()
{
	for(auto i : degreeVertex[1])
	{
		if(marked[i] == -1)
		{
			marked[i] = -2;
			int u = *(adjacency[i].begin());
			
			adjacency[u].erase(i);
			if(degree[u] < MAX_DEGREE_SET)
				degreeVertex[degree[u]].erase(u);
			degree[u]--;
			if(degree[u] < MAX_DEGREE_SET)
				degreeVertex[degree[u]].insert(u);
			
			graphSet(u, i, 0);
			graphSet(i, u, 0);
			
			degreeVertex[1].erase(i);
			kernalizationDeg1();
			if(degree[u] == 2)
				kernalizationDeg2(u);
				
			return;
		}
	}
}

void kernalization()
{
	kernalizationDeg2();
	kernalizationDeg1();
}


//****************************************** EMV **************************************************

set<pair<int, int> > shortestPath(int start, int param1) //param1 is DEST if param2 is -1. Else Param1 is the "Sub division" as per the new dp table. See its usage for more sense.
{
	priority_queue<pair<int, int>> Q; 
	set<pair<int, int> > ans;
	Q.push(make_pair(0, start));
	shortest_paths[start][0] = 0;
	shortest_paths[start][1] = dijktraID;
	parent[start] = -1;
	while(!Q.empty()) 
	{
		pair<int, int> top = Q.top();
		Q.pop();
		//Pop the Value, Vertex from the queue (regular dijkstra)
		int v = top.second, d = top.first;
		if(shortest_paths[v][1] == dijktraID and d <= shortest_paths[v][0])
			for(auto v2 : adjacency[v])
			{
				int edgeValue = graphGet(v, v2);
				//Either the ID for that vertex is different or the vertex was already visited when run for this vertex and the value is less
				if(valid(edgeValue) and (shortest_paths[v2][1] != dijktraID) or (shortest_paths[v2][0] > shortest_paths[v][0] + edgeValue))
				{
					shortest_paths[v2][0] = shortest_paths[v][0] + edgeValue;
					shortest_paths[v2][1] = dijktraID;
					parent[v2] = v;
					Q.push(make_pair(shortest_paths[v2][0], v2));
				}
			}
	}
	//param not equal to -1 means you want the path from start to param1 to be stored in the global variable "ans"
	if(param1 != -1)
	{
		int dest = param1;
		while(dest != start)
		{
			ans.insert(mp(min(dest, parent[dest]), max(dest, parent[dest])));
			dest = parent[dest];
		}
	}
	dijktraID++;
	return ans;
}

int getEMVDP(int a, int b, int c)
{
	int l1 = totalNodes + 5;
	if(a*l1 + b >= MAX_NODES)
		while(1);
	return EMV_dp[a*l1 + b][c];
}

void setEMVDP(int a, int b, int c, int d)
{
	int l1 = totalNodes + 5;
	if(a*l1 + b >= MAX_NODES)
		while(1);
	EMV_dp[a*l1 + b][c] = d;
}

void EMV_GetSolution(int D, int l) //Recurse on the two subtrees and retrieve the solution
{
	if(D == 0 or  l == 0)
		return;
	if(getEMVDP(D, l, 2) == D or getEMVDP(D, l, 2) == 0)
	{
		set<pair<int, int> > S2 = shortestPath(l, getEMVDP(D, l, 1));
		solution.insert(S2.begin(), S2.end());
		EMV_GetSolution(getEMVDP(D, l, 2), getEMVDP(D, l, 1));
	}
	else
	{
		EMV_GetSolution(getEMVDP(D, l, 2), getEMVDP(D, l, 1));
		EMV_GetSolution(getEMVDP(D, l, 2) ^ D, getEMVDP(D, l, 1));
	}
}

int Solver()
{
	kernalization();
	int K = (1<<cterminals) - 1, ans  = 1e9+7;
	for(int i = 0; i < cterminals; i++)
	{
		int posn = 1<<i;
		shortestPath(terminals[i], -1);
		for(int j = 1; j <= totalNodes; j++)
			if(marked[j] != -2)
			{
				setEMVDP(posn, j, 0, shortest_paths[j][0]);
				setEMVDP(posn, j, 1, terminals[i]);
				setEMVDP(posn, j, 2, 0);
				setEMVDP(posn, j, 3, 1e9+7);
			}
	}
	for(int D = 3; D <= K; D++)
	{
		if(__builtin_popcount(D) > 1)
		{
			for(int v = 1; v <= totalNodes; v++)
				if(marked[v] != -2)
				{
					setEMVDP(D, v, 0, 1e9+7);
					setEMVDP(D, v, 3, 1e9+7);
					for(int D1 = D; D1 > 0; D1 = (D1-1)&D)
					{
					 	if(getEMVDP(D1, v, 0) + getEMVDP(D ^ D1, v, 0) < getEMVDP(D, v, 3))
					 	{
							setEMVDP(D, v, 3, getEMVDP(D1, v, 0) + getEMVDP(D ^ D1, v, 0));
							setEMVDP(D, v, 2, D1);
							setEMVDP(D, v, 1, v);
						}
					}
				}
				for(int v = 1; v <= totalNodes; v++)
					if(marked[v] != -2)
					{
						graphSet(auxillary_node, v, getEMVDP(D, v, 3));
					}
				for(int i = 0; i < cterminals; i++)
				{
					if((D & (1<<i)) != 0) //i is present in the set denoted by D
					{
						graphSet(auxillary_node, terminals[i], getEMVDP(D^(1<<i), terminals[i], 0));
					}
				}
				shortestPath(auxillary_node, -1);
				for(int v = 1; v <= totalNodes; v++)
					if(marked[v] != -2)
					{
						if((marked[v] == -1) or ((D&(1<<marked[v]))==0)) //either its not a terminal or it doesnt belong to D.
						{
							setEMVDP(D, v, 0, shortest_paths[v][0]);
							if(parent[v] != auxillary_node)
							{
								setEMVDP(D, v, 1, parent[v]);
								setEMVDP(D, v, 2, D);
							}
						}
					}
				
		}
	 }
	ans = getEMVDP(K^1, terminals[0], 0);
	 EMV_GetSolution(K^1, terminals[0]);
	 cout<<"VALUE "<<ans<<"\n";
	 for(auto i : solution)
	{
	 	printAlias(i.first, i.second);
	}
	for(auto i : toPrint)
	{
		cout<<i.first<<" "<<i.second<<endl;
	}
}	


int main()
{
	getInput();
	Solver();
	return 0;	
}
