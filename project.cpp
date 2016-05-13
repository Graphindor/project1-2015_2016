#include <iostream>
#include <vector>
#include <string>
#include <list>

#include <limits> // for numeric_limits
 
#include <set>
#include <utility> // for pair
#include <algorithm>
#include <iterator>
#include <new>
#include <fstream>

using namespace std;
 
typedef int vertex_t;
typedef double weight_t;
 
const weight_t max_weight = numeric_limits<double>::infinity();
 
struct neighbor {
    vertex_t target;
    weight_t weight;
    neighbor(vertex_t arg_target, weight_t arg_weight)
        : target(arg_target), weight(arg_weight) { }
};
 
typedef vector<vector<neighbor> > adjacency_list_t;
 
 
void DijkstraComputePaths(vertex_t source,
                          const adjacency_list_t &adjacency_list,
                          vector<weight_t> &min_distance,
                          vector<vertex_t> &previous)
{
    int n = adjacency_list.size();
    min_distance.clear();
    min_distance.resize(n, max_weight);
    min_distance[source] = 0;
    previous.clear();
    previous.resize(n, -1);
    set<pair<weight_t, vertex_t> > vertex_queue;
    vertex_queue.insert(make_pair(min_distance[source], source));
 
    while (!vertex_queue.empty()) 
    {
        weight_t dist = vertex_queue.begin()->first;
        vertex_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());
 
        // Visit each edge exiting u
    const vector<neighbor> &neighbors = adjacency_list[u];
        for (vector<neighbor>::const_iterator neighbor_iter = neighbors.begin();
             neighbor_iter != neighbors.end();
             neighbor_iter++)
        {
            vertex_t v = neighbor_iter->target;
            weight_t weight = neighbor_iter->weight;
            weight_t distance_through_u = dist + weight;
        if (distance_through_u < min_distance[v]) {
            vertex_queue.erase(make_pair(min_distance[v], v));
 
            min_distance[v] = distance_through_u;
            previous[v] = u;
            vertex_queue.insert(make_pair(min_distance[v], v));
 
        }
 
        }
    }
}
 
 
vector<vertex_t> DijkstraGetShortestPathTo(
    vertex_t vertex, const vector<vertex_t> &previous)
{
    vector<vertex_t> path;
    for ( ; vertex != -1; vertex = previous[vertex])
        path.insert(path.begin(), vertex);
    return path;
}

bool ContainsPair(vector<pair<pair<int,int>,int>> &v, pair<int,int> arc)
{
    bool retval = false;

    for(int i = 0; i < v.size() && retval == false; i++)
        if(v[i].first.first == arc.first && v[i].first.second == arc.second)
            retval = true;

    cout << "v[] ContainsPair of " << arc.first << " " << arc.second << " is " << retval << endl;

    return retval;
}

void AddToAdjacencyList(adjacency_list_t &adjacency_list, vector<pair<pair<int,int>,int>> &v, int index)
{
    bool found = false;
    
    for(int i = 0; i < adjacency_list[v[index].first.first].size() && found == false; i++)
        if(adjacency_list[v[index].first.first][i].target == v[index].first.second)
            found = true;

    if(found == false)
        adjacency_list[v[index].first.first].push_back(neighbor(v[index].first.second,v[index].second)); 
}

vector <int> state;
vector <neighbor> parent;
bool t = 1; 
int theNodeInTheCycle;
int theWeightInTheCycle;

void dfs(int x, adjacency_list_t &ls)
{
    state[x] = 1;
    cout << "x==>> " << x << endl;
    cout << "I'm visiting " << x << " with size => " << ls[x].size() << endl;
    if(ls[x].size() > 0)
    for(int j = 0; j < ls[x].size(); j++)
    {
        cout << "ls[" << x << "][" << j << "].target => " << ls[x][j].target << endl;
        cout << "parent[" << x << "].target => " << parent[x].target << endl;
        if(state[ls[x][j].target] == 1 && parent[x].target != ls[x][j].target)
        {   
            cout << "Closed cycle since state[ls["<< x <<"]["<<j<<"].target] => " << state[ls[x][j].target] << endl;
            parent[ls[x][j].target].target = x;
            theNodeInTheCycle = ls[x][j].target; //ls[x][j] belongs to the cycle since state[ls[x][j]]==1
            t = 0;
        }
        if(state[ls[x][j].target] == 0)// && parent[x].target != -1)
        {
            parent[ls[x][j].target].target = x;
            dfs(ls[x][j].target, ls);
        }
    }
}

vector <neighbor> GetCycle ()
{
    vector <neighbor> cycle;
    int firstNodeInTheCycle = theNodeInTheCycle;
    do 
    {
            theNodeInTheCycle = parent[theNodeInTheCycle].target;
            theWeightInTheCycle = parent[theNodeInTheCycle].weight;
            cycle.push_back (neighbor(theNodeInTheCycle,theWeightInTheCycle));
            cout << "theNodeInTheCycle => " << theNodeInTheCycle << " firstNodeInTheCycle => " << firstNodeInTheCycle << endl;
    } while (theNodeInTheCycle != firstNodeInTheCycle && theNodeInTheCycle != -1);

    reverse(cycle.begin(), cycle.end()); //to get them in the right order
    if(theNodeInTheCycle == -1)
        cycle.clear();
    
    return cycle;
}

int main()
{
    ifstream in("input.txt");
    
    ofstream out("output.txt");

    int N, M, L, K;

    in >> N >> M >> L >> K;
    
    vector<int> entrances;

    for(int i = 0; i < L; i++)
    {
        int entry;
        in >> entry;
        entrances.push_back(entry);
    }

    // remember to insert edges both ways for an undirected graph
    adjacency_list_t adjacency_list(M);
    state.resize(N);

    for(int i = 0; i < N; i++)
        parent.push_back(neighbor(-1,-1));

    for(int i = 0; i < M; i++)
    {
        int from, to, weight;
        in >> from >> to >> weight;
        adjacency_list[from].push_back(neighbor(to,weight));
    }

    vector<int> answers;
    cout << "starting finding paths " << endl;
    for(int i = 0; i < entrances.size(); i++)
    {
        int evil_moves = K;
        vector<weight_t> min_distance;
        vector<vertex_t> previous;
        
        DijkstraComputePaths(entrances[i], adjacency_list, min_distance, previous);

        cout << "Distance from " << entrances[i] <<" to " << N - 1 << ": " << min_distance[N - 1] << endl;
        vector<vertex_t> best_path = DijkstraGetShortestPathTo(N - 1, previous);
        int best_dist = min_distance[N - 1];
        cout << "Best distance => " << best_dist << endl;
        cout << "Path : ";
        copy(best_path.begin(), best_path.end(), ostream_iterator<vertex_t>(cout, " "));
        cout << endl;

        vector<pair<pair<int,int>,int>> v;
        //int residual_path = 0;
        int minimum_dist = min_distance[N - 1];
        bool evil_found = false;

        if(K == 1 && evil_found == false)
        {
            cout << "Searching for evil moves" << endl;
            for(int j = best_path.size() - 2; j >= 0; j--)
            {
                vector<pair<pair<int,int>,int>> backup;
                cout << "Trying on " << best_path[j] << " neighbors => " << adjacency_list[best_path[j]].size() << endl;
                if(adjacency_list[best_path[j]].size() > 1) /*&& path[j] != entrances[i]*/// && ContainsPair(v,make_pair(best_path[j], best_path[j + 1])) == false)
                {
                    for(int z = 0; z < adjacency_list[best_path[j]].size(); z++)
                    {
                        int uno = adjacency_list[best_path[j]][z].weight;
                     
                        cout << "Saving " << best_path[j] << " and " << adjacency_list[best_path[j]][z].target << " with weight => " << uno << endl;
                        backup.push_back(make_pair(make_pair(best_path[j], adjacency_list[best_path[j]][z].target), uno));
                    }

                    int weight;
                    //cout << "\tNeighbors count => " << adjacency_list[best_path[j]].size() << endl;
                    for(int k = 0; k < adjacency_list[best_path[j]].size() && evil_found == false; k++)
                    {
                        cout << "Isolating target on " << best_path[j] << " => " << adjacency_list[best_path[j]][k].target << endl;
                        
                        bool found_isolated = false;
                        for(int m = 0; m < adjacency_list[best_path[j]].size() && found_isolated == false; m++)
                            if(m == k)
                            {
                                found_isolated = true;
                                weight = adjacency_list[best_path[j]][m].weight;
                     
                                cout << "Pushing back " << best_path[j] << " and " << adjacency_list[best_path[j]][m].target << " with weight => " << weight << endl;
                                v.push_back(make_pair(make_pair(best_path[j], adjacency_list[best_path[j]][m].target), weight));
                            }

                        adjacency_list[best_path[j]].clear();

                        cout << "\tAdding => " << v[0].first.first << " to " << v[0].first.second << endl;
                        adjacency_list[v[0].first.first].push_back(neighbor(v[0].first.second, v[0].second));
                            
                        /*
                        for(int z = 0; z < adjacency_list[best_path[j]].size(); z++)
                            if(z != k)
                                adjacency_list[best_path[j]].erase(adjacency_list[best_path[j]].begin() + z);
                        */
                        for(int u = 0; u < adjacency_list[best_path[j]].size(); u++)
                            cout << "\t\tRemaining neighbors " << adjacency_list[best_path[j]][u].target << endl;
                        
                        DijkstraComputePaths(entrances[i], adjacency_list, min_distance, previous);
                        vector<vertex_t> my_path = DijkstraGetShortestPathTo(N - 1, previous);
                        
                        cout << "New distance from " << entrances[i] <<" to " << N - 1 << ": " << min_distance[N - 1] << endl;
                        /*cout << "Path : ";
                        copy(my_path.begin(), my_path.end(), ostream_iterator<vertex_t>(cout, " "));
                        cout << endl;
                        */
                        if(min_distance[N - 1] > minimum_dist)
                        {
                            minimum_dist = min_distance[N - 1];
                            cout << "\t$$$New minimum_dist => " << minimum_dist << endl;
                        }

                        if(min_distance[N - 1] == max_weight)
                        {
                            evil_found = true;
                            //TODO: compute circle circumference
                            cout << "Searching for a cycle " << endl;
                            dfs(entrances[i], adjacency_list);
                            
                            if (t==0) 
                            {
                                cout << "Cycle Cycle Cycle Cycle Cycle Cycle" << endl;
                                vector <neighbor> cycle = GetCycle ();
                                for (int c = 0; c < cycle.size (); ++c)
                                    cout << cycle[c].target << " ";
                                cout << "\n";

                                int circumference = 0;
                                for (int c = 0; c < cycle.size (); ++c)
                                    circumference += cycle[c].weight;
                                circumference -= 1;
                                circumference *= -1;
                                cout << "circumference => " << circumference << endl;
                                if(circumference != 1)
                                    minimum_dist = circumference + best_dist;
                                else
                                    minimum_dist = -1;
                            }
                            else 
                            {
                                cout << "No cycle\n";
                                minimum_dist = -1;
                            }
                            t = 1;
                        }

                        //TODO: Create a backup vector<pair<pair<int,int>,int> to prevent the re-ordering of the adjacencies
                        cout << "\tv.size() => " << v.size() << endl;
                        adjacency_list[best_path[j]].clear();
                        
                        for(int l = 0; l < backup.size(); l++)
                        {
                            //cout << "\tRestoring => " << backup[l].first.first << " to " << backup[l].first.second << endl;
                            adjacency_list[backup[l].first.first].push_back(neighbor(backup[l].first.second, backup[l].second));
                            //AddToAdjacencyList(adjacency_list, v, l);
                        }
                        v.clear();
                    }
                }
                backup.clear(); 
            }
        }

        v.clear();

        if(minimum_dist == max_weight || minimum_dist < 0)
            answers.push_back(-1);
        else
            answers.push_back(minimum_dist);
        
        cout << ">>>>>>>>>>>>>>>>>>>Pushing back answer => " << min_distance[N - 1] << "<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    }
    
    copy(answers.begin(), answers.end(), ostream_iterator<int>(out, "\n"));
    
    cout << "answers : " << endl;
    copy(answers.begin(), answers.end(), ostream_iterator<int>(cout, "\n"));
    cout << endl;
    
    return 0;
}