#include <iostream>
#include <list>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// This class represents a directed graph using
// adjacency list representation
class Graph
{
    int V; // No. of vertices

    // Pointer to an array containing adjacency
    // lists
    list<int> *adj;
public:
    Graph(int V); // Constructor

    //read the adjacency list
    void readFile();

    //store the adjacency list
    void storeFile();

    // function to add an edge to graph
    void addEdge(int v, int w);

    //generate the graph
    void generateGraph();

    // prints BFS traversal from a given source s
    void BFS(int s);

    //print the adjacency list
    void printList();
};

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

void Graph::readFile()
{
    ifstream infile = ifstream("data.txt");
    string str;

    for(int i = 0; i < V; i++)
    {
        getline(infile, str);
        istringstream s(str);
        int data = 0;
        while(s >> data)
        {
            addEdge(i, data);
        }
    }

    infile.close();
}

void Graph::storeFile()
{
    generateGraph();
    ofstream out = ofstream("data.txt");

    //write file
    for(int i = 0 ; i < V; i++)
    {
        list<int>::iterator j;
        for(j = adj[i].begin(); j != adj[i].end(); j++)
        {
            out<<*j<<" ";
        }
        out<<endl;
    }

    out.close();
    adj->clear();
}

void Graph::generateGraph()
{
    for(int i = 0; i < V; i++)
    {
        //generate the adjacency list
        for(int j = 0; j < V; j++)
        {
            if(j == i)
            {}
            else {
                int randomNum = rand() % 10 + 1;
                //cout<<randomNum<<" ";
                if (randomNum > 5) {
                    addEdge(i, j);
                }
            }
        }
    }
}

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w); // Add w to v’s list.

    //generate the graph adjacency list

}

void Graph::BFS(int s)
{
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for(int i = 0; i < V; i++)
        visited[i] = false;

    // Create a queue for BFS
    list<int> queue;

    // Mark the current node as visited and enqueue it
    visited[s] = true;
    queue.push_back(s);

    // 'i' will be used to get all adjacent
    // vertices of a vertex
    list<int>::iterator i;

    while(!queue.empty())
    {
        // Dequeue a vertex from queue and print it
        s = queue.front();
        //cout << s << " ";
        queue.pop_front();

        // Get all adjacent vertices of the dequeued
        // vertex s. If a adjacent has not been visited,
        // then mark it visited and enqueue it
        for (i = adj[s].begin(); i != adj[s].end(); ++i)
        {
            if (!visited[*i])
            {
                visited[*i] = true;
                queue.push_back(*i);
            }
        }
    }
}

void Graph::printList()
{
    list<int>::iterator j;
    for(int i = 0; i < V; i++)
    {
        cout<<"Node is :"<< i <<". And following is :";
        for(j = adj[i].begin(); j != adj[i].end(); j++)
        {
            cout<<*j<<" ";
        }
        cout<<endl;
    }
}

int main()
{
    clock_t start, end;
    start = clock();
    // Create a graph given in the above diagram
    Graph g(1000);

    g.storeFile();

    g.readFile();

    cout << "Following is Breadth First Traversal "
         << "(starting from vertex 0) \n";
    g.BFS(0);

    cout<<endl;

   //g.printList();

    end = clock();
    cout<<"The running time is :"<<(double)(end - start) / CLOCKS_PER_SEC<<endl;

    return 0;
}
