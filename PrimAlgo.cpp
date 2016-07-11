#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <list>
#include <sstream>
#include <iomanip>

using namespace std;

#define SIZE     50     // Number of Vertices
#define TRUE      1
#define FALSE     0
#define MAX    9999

/* Grraph class */
class Graph{
	private:
		int numV;
		int numE;
		int matrix[SIZE][SIZE];
		int node[SIZE];
	public:
		Graph();                         // Constructor 1
		Graph(int numV, int numE);       // Constructor 2
		int getV();                      // Get number of vertices
		int getE();                      // Get number of edges
		void setE(int E);                // Set number of edges
		void setV(int V);
		bool adjacent(int src, int dest);  // Check if the two vertices have an edge between them
		void neighbors(int src, vector<int>& neigh);  // Return all the neighbors of given vertex
		void addEdge(int src, int dest);        // Add an edge between two vertices
		void deleteEdge(int src, int dest);     // Delete the edge between two vertices
		int get_node_value(int src);            // Get node value of vertex
		void set_node_value(int src, int value); // Set node value of vertex
		int get_edge_value(int src, int dest);   // Get edge or distance between two vertices
		void set_edge_value(int src, int dest, int value);  // Set edge or distance between two vertices
		void neighborsList(vector<int> &neigh, int u);
		bool isExists(vector<int> &neigh, int u);
};

// Constructor 1
Graph::Graph()
{
	numV = SIZE;
	numE = (SIZE) * (SIZE-1) / 2;
	for(int i = 0; i < numV; i++){
		for (int j = 0; j < numV; j++){
			matrix[i][j] = 0;
		}
		node[i] = MAX;
	}
}

// Constructor 2
Graph::Graph(int V, int E)
{
	numV = V;
	numE = E;
	for(int i = 0; i < numV; i++){
		for (int j = 0; j < numV; j++){
			matrix[i][j] = 0;
		}
		node[i] = MAX;
	}		
}

// Get number of vertices
int Graph::getV()
{
	return numV;
}

// Get number of edges
int Graph::getE()
{
	return numE;
}

// Set number of edges
void Graph::setE(int E)
{
	numE = E;
}

// Check if the two vertices have an edge between them
bool Graph::adjacent(int src, int dest)
{
	if (matrix[src][dest] == 0)
		return FALSE;
	else 
		return TRUE;
}

// Return all the neighbors of given vertex
void Graph::neighbors(int src, vector<int> &neigh)
{
	for(int i = 0; i < numV; i++){
		if (matrix[src][i])
			neigh.push_back(i);
	}
}

// Add an edge between two vertices
void Graph::addEdge(int src, int dest)
{
	matrix[src][dest] = 101;
}

// Delete the edge between two vertices
void Graph::deleteEdge(int src, int dest)
{
	matrix[src][dest] = 0;	
}

// Get node value of vertex
int Graph::get_node_value(int src)
{
	return node[src];
}

// Set node value of vertex
void Graph::set_node_value(int src, int value)
{
	node[src] = value;
}

// Get edge or distance between two vertices
int Graph::get_edge_value(int src, int dest)
{
	return matrix[src][dest];
}

// Set edge or distance between two vertices
void Graph::set_edge_value(int src, int dest, int value)
{
	matrix[src][dest] = value;
	matrix[dest][src] = value;
}

void Graph::setV(int V)
{
	numV = V;
}

bool Graph::isExists(vector<int> &neigh, int v)
{
	for(int i = 0; i < numV; i++){
		cout << i << " ";
		if (v == neigh[i])
			return 1;
	}	
	return 0;
}

void Graph::neighborsList(vector<int> &neigh, int u)
{
	int exists = 0;
	for(int i = 0; i < numV; i++){
		if((matrix[u][i] != 0)){ 
			for(int j = 0; j < neigh.size(); j++){
				//cout << neigh[j];
				if (i == neigh[j])
					exists = 1;
			}
			if (exists == 0)
				neigh.push_back(i);
			exists = 0;
		}
	}		
}

/* Structure to represent node in priority queue */
struct node{
	int vertex;   // vertex
	int cost;	  // node value
	node *next;
};

/* Priority Queue Class: Implemented similiar to min heap */
class priorityQ{
	private:
		struct node *root;     //head of priority queue
		int size;              //number of elements in the queue
	public:
		priorityQ();                   // Constructor 1
		priorityQ(int src, int cost);  // Constructor 2
		~priorityQ();                  // Destructor
		void chgPriority(int src, int cost);   // Change priority or node value of give vertex
		void minPriority();             // Delete the first element in the priority queue
		bool contains(int src);         // Check if the priority queue contains the given vertex
		void insert(int src, int cost);   // Insert the vertex in the priority queue
		int top();     // Return the top priority vertex in the queue
		int sizeQ();    // Return the size of priority queue
		bool isEmpty();  // Check if the priority queue is empty or not
};

// Constructor 1
priorityQ::priorityQ()
{
	root = NULL;	
	size = 0;
}

// Constructor 2
priorityQ::priorityQ(int src, int cost)
{
	root = new node;
	root->vertex = src;
	root->cost = cost;	
	root->next = NULL;
	size = 1;
}

// Destructor
priorityQ::~priorityQ()
{
	struct node *temp = root;
	if(root){
		temp = root;
		root = root->next;
		delete temp;
    }
    root = NULL;
    size = 0;
}

// Change priority or node value of give vertex
void priorityQ::chgPriority(int src, int cost)
{
	node *temp = root;
	while (temp != NULL){
		if (temp->vertex == src){
			temp->cost = cost;
			break;
		}
		temp = temp->next;
	}
	
	if(temp == NULL)
		return;
	
	/* After changing priority, compute min heap again to find the new root */
	temp = root;
	node *min = root, *prev = NULL, *p = NULL;
	int count = 0;
	while(temp != NULL && (count < size)){
		if(min->cost > temp->cost){
			p = prev;
			min = temp;
		}
		prev = temp;
		temp = temp->next;
		count++;
	}
	
	if (min != root){
		p->next = min->next;
		min->next = root;
		root =  min;
	}
}

// Delete the first element in the priority queue
void priorityQ::minPriority()
{
	node *temp1 = root;
			
	if (root){
		temp1 = root;
		size = size - 1;
        if (size == 0){
			delete root;
			root = NULL;
			return;
		}
		root = root->next;
		delete temp1;
		temp1 = NULL;
	}
	
	if(size <= 0){
		root = NULL;
		return;
	}
	
	if(root == NULL){
		size = 0;
		return;
	}
	
	if (root->next == NULL)
		return;
		
	/* After deleteing the root, compute min heap again to find the new root */
	node *temp = root;
	node *min = root, *prev = NULL, *p = NULL;
	int count = 0;
	while(temp != NULL && (count < size)){
		if(min->cost >= temp->cost){
			p = prev;
			min = temp;
		}
		prev = temp;
		temp = temp->next;
		count++;
	}
	
	if (min != root){
		p->next = min->next;
		min->next = root;
		root =  min;
	}
}

// Check if the priority queue contains the given vertex
bool priorityQ::contains(int src)
{
	node *temp = root;
	while(temp != NULL){
		if (temp->vertex == src)
			return TRUE;
	}	
	return FALSE;
}

// Insert the vertex in the priority queue
void priorityQ::insert(int src, int cost)
{
	node *newNode = new node;
    newNode->vertex = src;
	newNode->cost = cost;
    if (root != NULL){
		if (root->cost < cost){
			newNode->next = root->next;
			root->next = newNode;
		} else{
			newNode->next = root;
			root = newNode;
		}
		size += 1;
	} else{
		root = newNode;
		size = 1;
	}
}

// Return the top priority vertex in the queue
int priorityQ::top()
{
	if (root)
		return root->vertex;
	else 
		return -1;
}

// Return the size of priority queue
int priorityQ::sizeQ()
{
	return size;
}

// Check if the priority queue is empty or not
bool priorityQ::isEmpty()
{
	if (size <= 0)
		return TRUE;
	return FALSE;
}

/* Shortest Path Class: Compute Dijstra's Algorithm for the Graph */
class ShortestPath{
	public: 
	    void path(Graph& G, int u, int w); // Return the shortest distance between two vertices and also the path to reach each other
		float path_size(Graph G);    // Return average shortest distance for graph
};

// Return the shortest distance between two vertices and also the path to reach each other
void ShortestPath::path(Graph& G, int u, int w)
{
	int path[50];
	
	G.set_node_value(u, 0);  // Mark source nod evalue as zero

    // Insert all vertices in the queue
	priorityQ Q;
	for(int i = 0; i < G.getV(); i++){
		Q.insert(i, G.get_node_value(i));
	}

	int vertex;
	while(!Q.isEmpty()){
		vertex = Q.top();    // Will be source for the first time
		Q.minPriority();
		vector<int> neigh;
		G.neighbors(vertex, neigh);  // Find all neighbors
		
		/* For each neighbor check if the node has the shortest possible value */
		for(int i = 0; i < neigh.size(); i++){
			int temp = G.get_node_value(vertex) + G.get_edge_value(vertex, neigh[i]);
			if (temp < G.get_node_value(neigh[i])){
				G.set_node_value(neigh[i], temp);
				Q.chgPriority(neigh[i], temp);
				path[neigh[i]] = vertex;
			}
		}
		
		/* If we already found the destination node, Return */
		if (vertex == w)
			break;
	}
	
	cout << "Path is: ";
	
	int dest = w, src = u;
	int newPath[1000], cnt = 0;
	// Edge exists
	if (G.get_edge_value(src, dest) != 0){ 
		while(dest != src){
			int temp = path[dest];
			newPath[cnt] = temp;
			cnt++;
			dest = temp;
		}
		for(int k = cnt-1; k >= 0; k--)
			cout << newPath[k] << " ";
		cout << w;
		cout << endl << "Shortest Distance between " << u << " and " << w << " is: " << G.get_node_value(w) << endl;
	} else {     //Edge does not exist
		cout << "No path exists";
	}
}

// Return average shortest distance for graph
float ShortestPath::path_size(Graph G)
{
	G.set_node_value(0, 0);  // Will be zero(th node) for the first time
	
	priorityQ Q;
	for(int i = 0; i < G.getV(); i++){
		Q.insert(i, G.get_node_value(i));
	}
	
	int vertex;
	while(!Q.isEmpty()){
		vertex = Q.top();
		Q.minPriority();

		vector<int> neigh;
		G.neighbors(vertex, neigh);   // Find all neighbors
		
		/* For each neighbor check if the node has the shortest possible value */
		for(int i = 0; i < neigh.size(); i++){
			int temp = G.get_node_value(vertex) + G.get_edge_value(vertex, neigh[i]);
			if (temp < G.get_node_value(neigh[i])){
				G.set_node_value(neigh[i], temp);
				Q.chgPriority(neigh[i], temp);
			}
		}
	}
	
	/* Find the average cost for the graph depending on the number of edges */
	float avgEcost = 0;
	for(int i = 0; i < G.getV(); i++){
		for(int j = 0; j < G.getV(); j++){
			int tempCost = G.get_edge_value(i, j);
			avgEcost += tempCost;
		}
	}
	
	avgEcost = avgEcost / (2 * G.getE());
	return avgEcost;
}

/* Monte Carlo Simulation to generate Graph */
void generateGraph(Graph &G, float edgeDensity, int distanceRange)
{
	int v = G.getV();
	int E = (v * (v - 1)) / 2;
	int numE = E * edgeDensity;
	G.setE(numE);
	
	int rand1, rand2;
	int i = 0;
	while(i < numE){
		rand1 = rand() % v;
		rand2 = rand() % v;
		if(rand1 == rand2)
			continue;
		if(G.get_edge_value(rand1, rand2) != 0)
			continue;
		int cost = rand() % distanceRange + 1;
		/* Print Random Edges */
		//cout << "i: " << rand1 << " j: " << rand2 << " cost: " << cost << endl;
		G.set_edge_value(rand1, rand2, cost);
		i += 1;
	}
}

void dijikstraAlgo()
{
	ShortestPath S;
	Graph G;
	float density = 0;
	float distanceRange = 100.0;
	
	cout << "***********DIJIKSTRA'S ALGORITHM*************" << endl;
	
	/* Generate graph using monte carlo simulation 
	   for given density with distance range of 1.0-100.0 */
	cout << "Enter the density (eg: 0.1 or 0.2 or etc. .) ";
	cin >> density;
	cout << "The distance range for the edge cost is: 1.0 to 100.0" << endl;
	cout << "The number of vertices is always considered 50." << endl;
	cout << "If the number of vertices need to be changed, then change the value of SIZE in the program" << endl;
	generateGraph(G, density, distanceRange);
		
	/* Get the average shortest distance for the graph using Dijstra's algorithm */
	int dist = S.path_size(G);
	cout << endl;
	cout << "Average edge value after Dijkstra's algorithm: " << dist << "\n";
	
	/* Get the shortest distance for the given two vertices using Dijstra's algorithm */
	int src = 0, dest = 0;
	cout << endl;
	cout << "Enter the src and dest for which shortest distance and path is needed: ";
	cin >> src;
	cin >> dest;
	S.path(G, src, dest);
}

/* MST class */
class minimumSpanning{
	public:
		void primSpanning(Graph G);	
        int findMin(vector<int> visited, vector<int> dist, Graph G);
        void printMst(vector<int> visited, vector<int> parent, Graph G);
};

/* Find edge with minimum distance */
int minimumSpanning::findMin(vector<int> visited, vector<int> dist, Graph G)
{
	int numV = G.getV();
	int min = MAX;
	int idx = -1;
	
	for(int i = 0; i < numV; i++){
		if((visited[i] == FALSE) && (dist[i] < min)){
			min = dist[i];
			idx = i;
		}
	}
	
	return idx;
}

/* Printing the MST */
void minimumSpanning::printMst(vector<int> visited, vector<int> parent, Graph G)
{
	cout << "   Edge           Weight" << endl;
	for(int i = 1; i < visited.size(); i++){
		cout << setfill('0') << setw(2) << i << "    "; 
		cout << setfill('0') << setw(2) << parent[i] << "            "; 
		cout << setfill('0') << setw(2) << G.get_edge_value(parent[i], i);
		cout << endl;
	}	
}

/* Find the Miniumum Spanning Tree for a given graph using Prim's Algorithm */
void minimumSpanning::primSpanning(Graph G)
{
	vector<int> visited;
	vector<int> parent;
	vector<int> dist;
	int numV = G.getV();
	for(int i = 0; i < numV; i++){
		visited.push_back(FALSE);
		parent.push_back(-1);
		dist.push_back(MAX);
	}
	
	dist[0] = 0;
	for(int i = 0; i < numV; i++){
		int u = findMin(visited, dist, G);
		visited[u] = TRUE;
		
		for(int j = 0; j < numV; j++){
			int val = G.get_edge_value(u, j);
			if((visited[j] == FALSE) && (val != 0) && (val < dist[j])){
				dist[j] = val;
				parent[j] = u;
			}
		}
	}
	
	printMst(visited, parent, G);
}

/* Generate graph from given input file */
void generateMinGraph(Graph &G, minimumSpanning &m)
{
	ifstream fd("sample.txt");
	string line;
	if(fd.is_open()){
		getline(fd, line);
		stringstream stream(line);
		int num;
		stream>>num;
		G.setV(num);
		while(getline(fd, line)){
			stringstream stream(line);
			int arr[3];
			for(int i = 0; i < 3; i++){
				stream>>arr[i];
			}
			G.set_edge_value(arr[0], arr[1], arr[2]);
		}
		fd.close();
	}
}

int main() {
	
	//dijikstraAlgo();
	
	Graph G;
	minimumSpanning m;
	
	cout << "***********PRIM'S ALGORITHM*************" << endl << endl;
	generateMinGraph(G, m);
	m.primSpanning(G);
	
	return 0;
}

