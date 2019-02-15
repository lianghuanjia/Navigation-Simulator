// Digraph.hpp
//
// ICS 46 Spring 2018
// Project #5: Rock and Roll Stops the Traffic
//
// This header file declares a class template called Digraph, which is
// intended to implement a generic directed graph.  The implementation
// uses the adjacency lists technique, so each vertex stores a linked
// list of its outgoing edges.
//
// Along with the Digraph class template is a class DigraphException
// and a couple of utility structs that aren't generally useful outside
// of this header file.
//
// In general, directed graphs are all the same, except in the sense
// that they store different kinds of information about each vertex and
// about each edge; these two types are the type parameters to the
// Digraph class template.

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <exception>
#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <queue>
#include <limits>

// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException : public std::runtime_error
{
public:
    DigraphException(const std::string& reason);
};


inline DigraphException::DigraphException(const std::string& reason)
    : std::runtime_error{reason}
{
}



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a struct template.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a struct template.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
};

/*
DigraphVertex<std::string, RoadSegment>
*/


// Digraph is a class template that represents a directed graph implemented
// using adjacency lists.  It takes two type parameters:
//
// * VertexInfo, which specifies the kind of object stored for each vertex
// * EdgeInfo, which specifies the kind of object stored for each edge
//
// You'll need to implement the member functions declared here; each has a
// comment detailing how it is intended to work.
//
// Each vertex in a Digraph is identified uniquely by a "vertex number".
// Vertex numbers are not necessarily sequential and they are not necessarily
// zero- or one-based.

template <typename VertexInfo, typename EdgeInfo>
class Digraph
{
public:
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The move constructor initializes a new Digraph from an expiring one.
    Digraph(Digraph&& d) noexcept;

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph() noexcept;

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    Digraph& operator=(Digraph&& d) noexcept;

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    std::vector<int> vertices() const;

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    std::vector<std::pair<int, int>> edges() const;

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    std::vector<std::pair<int, int>> edges(int vertex) const;

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    VertexInfo vertexInfo(int vertex) const;

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    EdgeInfo edgeInfo(int fromVertex, int toVertex) const;

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    void addVertex(int vertex, const VertexInfo& vinfo);

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    void addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo);

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead.
    void removeVertex(int vertex);

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    void removeEdge(int fromVertex, int toVertex);

    // vertexCount() returns the number of vertices in the graph.
    int vertexCount() const noexcept;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const noexcept;

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    int edgeCount(int vertex) const;

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    bool isStronglyConnected() const;

    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the precedessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    std::map<int, int> findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;
       

private:
    // Add whatever member variables you think you need here.  One
    // possibility is a std::map where the keys are vertex numbers
    // and the values are DigraphVertex<VertexInfo, EdgeInfo> objects.
    
    std::map< int, DigraphVertex<VertexInfo, EdgeInfo>> all_info;
    bool check_vertex_exists(int vertex) const;
    bool check_edge_exists(int fromVertex, int toVertex) const;
/*    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double,int>>, std::greater<std::pair<double, int>>> find_cloest_graph(
        int current_vertex, 
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;
*/    
    bool DFT(int vertex) const;
    bool DFTr(int vertex, std::map<int, bool>& visited_map) const;
    // You can also feel free to add any additional member functions
    // you'd like (public or private), so long as you don't remove or
    // change the signatures of the ones that already exist.
};



// You'll need to implement the member functions below.  There's enough
// code in place to make them compile, but they'll all need to do the
// correct thing instead.

template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph()
{
    
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(const Digraph& d)
{
    all_info = d.all_info;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(Digraph&& d) noexcept
{
    all_info.swap(d.all_info);

}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::~Digraph() noexcept
{
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(const Digraph& d)
{

    all_info = d.all_info;
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(Digraph&& d) noexcept
{
    all_info.swap(d.all_info);
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<int> Digraph<VertexInfo, EdgeInfo>::vertices() const
{
    std::vector<int> all_vec_num;
    for (typename std::map<int,DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it = all_info.begin(); it != all_info.end(); ++it)
    {
    all_vec_num.push_back(it -> first);
    }
    return all_vec_num;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges() const
{
    std::vector<std::pair<int, int>> vec_edges_info;
    for(typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
    {
        for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator it2 = it1->second.edges.begin(); it2 != it1->second.edges.end(); ++it2)
        {
        std::pair<int, int> single_pair_of_edge = std::make_pair(it2->fromVertex, it2->toVertex);
        vec_edges_info.push_back(single_pair_of_edge);
        }
    }
    return vec_edges_info;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges(int vertex) const
{
    std::vector<std::pair<int, int>> vec_edges_info;
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
        {
        if (it1 -> first == vertex)
            {
            for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator it2 = it1 -> second.edges.begin(); it2 != it1 -> second.edges.end(); ++it2)
                {
                std::pair<int,int> single_pair_of_edge = std::make_pair(it2->fromVertex, it2->toVertex);
                vec_edges_info.push_back(single_pair_of_edge);
                }
            }
        }
    if (vec_edges_info.size() == 0)
        {
        throw DigraphException("There is no such a vertex inside the all_info.(edges(int vertex)) function");
        }
    else
        {
        return vec_edges_info;
        }
}


template <typename VertexInfo, typename EdgeInfo>
VertexInfo Digraph<VertexInfo, EdgeInfo>::vertexInfo(int vertex) const
{   
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
    {
        if (it1 -> first == vertex)
            {
            return (it1 -> second.vinfo);
            }
    }
    throw DigraphException("There is no such a vertex inside the all_info.(vertexInfo(int vertex)) function"); 
}

template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::check_vertex_exists(int vertex) const
{
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
    {
        if (it1 -> first == vertex)
            {
            return true;
            }
    }
    return false;
}

template <typename VertexInfo, typename EdgeInfo>
EdgeInfo Digraph<VertexInfo, EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const
{
    for(typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
    {
        for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator it2 = it1 -> second.edges.begin(); it2 != it1 -> second.edges.end(); ++it2)
        {
            if (it2 -> fromVertex == fromVertex && it2 -> toVertex == toVertex && (check_vertex_exists(it2->fromVertex) == true) && (check_vertex_exists(it2 -> toVertex) == true))
            {
                return it2 -> einfo;
            }
        }
    }
    throw DigraphException("There is no such an EdgeInfo inside the all_info.(edgeInfo(int fromVertex, int toVertex) function)");
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo)
{
    if (check_vertex_exists(vertex) == true)
        {
        throw DigraphException("In ::addVertex, the vertex is already existed");
        }
    else
        {
        std::pair<int, DigraphVertex<VertexInfo, EdgeInfo>> add_object = std::make_pair(vertex, DigraphVertex<VertexInfo, EdgeInfo>{vinfo});
        // the line above, when I make pair, why can I just write std::make_pair(vertex, DigraphVertex{vinfo}) ???
        all_info.insert(add_object);
        //std::cout << "adding" << std::endl;
        }

}

template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::check_edge_exists(int fromVertex, int toVertex) const
{
for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
    {
    for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator it2 = it1->second.edges.begin() ; it2 != it1->second.edges.end(); ++it2 )
        {
            if (it2 -> fromVertex == fromVertex && it2 -> toVertex == toVertex)
            // should I also check the einfo? like "it2 -> einfo == einfo"
                {
                return true;
                }
        }
    }
return false;
}

template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo)
{
    if (check_vertex_exists(fromVertex) == false || check_vertex_exists(toVertex) == false)
        throw DigraphException("In ::addEdge, the fromVertex or the toVertex does not exist in the map");
    else if (check_edge_exists(fromVertex, toVertex) == true)
        throw DigraphException("In ::addEdge, this edge exists already");
    else
        {
        DigraphEdge<EdgeInfo> new_edge = DigraphEdge<EdgeInfo>{fromVertex, toVertex, einfo};
        all_info.at(fromVertex).edges.push_back(new_edge);
        }
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeVertex(int vertex)
{
    if (check_vertex_exists(vertex) == false)
        {
        throw DigraphException("In ::removeVertex, the vertex you want to remove does NOT exist");
        }
    else
        {
        for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1) // NOTE: because we are going to delete some element, so we are going to modify the map, therefore, we use iterate, NOT const_iterator;
            {
            for(typename std::list<DigraphEdge<EdgeInfo>>::iterator it2 = it1->second.edges.begin(); it2 != it1->second.edges.end();)
                {
                if (it2 -> toVertex == vertex)
                    {
                    it2 = it1 -> second.edges.erase(it2);
                    }
                else
                    {
                    ++it2;
                    }
                }
            }
        }
        all_info.erase(vertex);
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeEdge(int fromVertex, int toVertex)
{
    if (all_info.count(fromVertex) == 0 || all_info.count(toVertex) == 0)
        {
        throw DigraphException("In ::removeEdge, the fromVertex or toVertex does not exist");
        }
    else if (check_edge_exists(fromVertex, toVertex) == false)
        {
        throw DigraphException("In ::removeEdge, the edge you want to delete does not exist");
        }
    else
        {
        for (typename std::list<DigraphEdge<EdgeInfo>>::iterator it1 = all_info.at(fromVertex).edges.begin(); it1 != all_info.at(fromVertex).edges.end();)
            {
            if (it1 -> fromVertex == fromVertex && it1 -> toVertex == toVertex)
                {
                it1 = all_info.at(fromVertex).edges.erase(it1);
                }
            else
                {
                ++it1; // If I delete the current it1, will it lose and cannot iterate to the next element?
                }
            }
        }
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::vertexCount() const noexcept
{
    return all_info.size();
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount() const noexcept
{
    int count = 0;
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
    {
        count += it1 -> second.edges.size();
    }
    return count;
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount(int vertex) const
{
    if (all_info.count(vertex) == 0)
        throw DigraphException("In ::edgeCount(int vertex), the target vertex does not exist");
    else
        return all_info.at(vertex).edges.size();

}


template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::DFT(int vertex) const
{
    std::map<int, bool> visited_map;
    std::vector<int> vec_of_vertives = vertices();
    
    for (int i = 0; i < vec_of_vertives.size(); ++i)
    {
        visited_map[vec_of_vertives[i]] = false;
    }
    return DFTr(vertex, visited_map);

}




template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::DFTr(int vertex, std::map<int, bool>& visited_map) const
{
    visited_map[vertex] = true;
    for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator it1 = all_info.at(vertex).edges.begin(); it1 != all_info.at(vertex).edges.end(); ++it1)
    {
        if (visited_map.at(it1 -> toVertex) == false)
            DFTr(it1 -> toVertex, visited_map);
    //the DFTr return a boolean value, above, if I don't make a boolean variable to get the DFTr result, will it be bad?
    }
    for (std::map<int, bool>::const_iterator it1 = visited_map.begin(); it1 != visited_map.end(); ++it1)
    {
        if (it1 -> second == false)
            return false;
    }
    return true;
}




template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::isStronglyConnected() const
{
for(typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator it1 = all_info.begin(); it1 != all_info.end(); ++it1)
    {
    if (DFT(it1->first) == false)
        return false;
    }
return true;
}



template <typename VertexInfo, typename EdgeInfo>
std::map<int, int> Digraph<VertexInfo, EdgeInfo>::findShortestPaths(
    int startVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc) const
{

    std::map<int, bool> known;
    std::map<int, int> predecessor;
    std::map<int, double> distance;

    std::vector<int> vec_of_vertives = vertices();
    for (int i = 0; i < vec_of_vertives.size(); ++i)
    {
        known[vec_of_vertives[i]] = false;
        predecessor[vec_of_vertives[i]] = -1; // -1 ,means unknown
        distance[vec_of_vertives[i]] = std::numeric_limits<double>::infinity();
    }
    predecessor[startVertex] = startVertex; // the start vertex's predecessor is just itself.
    distance[startVertex] = 0;

    std::priority_queue<std::pair<double,int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double,int>>> pri_que;
    pri_que.push(std::make_pair(0,startVertex));
    
    while (pri_que.size() != 0)
    {
    std::pair<double, int> shortest_vertex = pri_que.top();
    pri_que.pop();
    int v = shortest_vertex.second;

    if (known[v] == false)        
        {
        known[v] = true;
        for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator it1 = all_info.at(v).edges.begin(); it1 != all_info.at(v).edges.end(); ++it1)
            {
            if (distance[(it1 -> toVertex)] > (distance[v] + edgeWeightFunc(it1 -> einfo)))
                {
                    distance[(it1 -> toVertex)] = distance[v] + edgeWeightFunc(it1 -> einfo);
                    predecessor[it1 -> toVertex] = v;
                    pri_que.push(std::make_pair(distance[(it1 -> toVertex)], it1 -> toVertex));
                }
            }
    
        }
    }
   // when you tarce back, remember some predecessor are -1. Maybe have seg fault
    return predecessor;
}



#endif // DIGRAPH_HPP

