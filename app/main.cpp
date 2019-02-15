// main.cpp
//
// ICS 46 Spring 2018
// Project #5: Rock and Roll Stops the Traffic
//
// This is the program's main() function, which is the entry point for your
// console user interface.
#include "InputReader.hpp"
#include "RoadMap.hpp"
#include "RoadMapReader.hpp"
#include "RoadMapWriter.hpp"
#include "RoadSegment.hpp"
#include "Trip.hpp"
#include "TripMetric.hpp"
#include "TripReader.hpp"
#include "Digraph.hpp"
#include <iomanip>
#include <stack>
#include <sstream>

std::function<double(const RoadSegment&)> get_distance = [](RoadSegment rs){return rs.miles;};
std::function<double(const RoadSegment&)> get_time = [](RoadSegment rs){return rs.miles/rs.milesPerHour;};

//std::function< return type ( parameter type )> function name = [] (parameter){return xxx;};

//final

double print_and_return_time(double miles, double mph)
{
    int hour = miles/mph;
    int min = (miles/mph - hour)*60;
    double sec = ((miles/mph - hour)*60 - min)*60;
    std::stringstream stream;
    if (hour == 0)
        {
        if (min == 0)
            if (sec == 0)
                {
                std::cout << "0.0 " << "secs)" << std::endl;
                }
            else
                {
                stream<<std::fixed<<std::setprecision(1)<<sec;
                std::cout << stream.str() << " secs)" << std::endl;
                }
        else
        {
            if (sec == 0)
                {
                std::cout << min << " mins " << "0.0" << " secs)" << std::endl;
                }
            else
                {
                stream<<std::fixed<<std::setprecision(1)<<sec;
                std::cout << min << " mins " << stream.str()<< " secs)" << std::endl;
                }
        }
        }
    else
        if (sec == 0)
            {
            std::cout << hour << " hrs " << min << " mins " << "0.0" << " sec)" << std::endl;
            }
        else
            {
            stream<<std::fixed<<std::setprecision(1)<<sec;
            std::cout << hour << " hrs " << min << " mins " << stream.str() << " sec)" << std::endl;
            }
    return (hour*3600 + min*60 + sec);
}




void time_conversion(double all_sec)
{
    int hour = all_sec/3600;
    int min = (all_sec-(hour*3600))/60;
    double sec = (all_sec-(hour*3600)-min*60);
    std::stringstream stream;
    if (hour == 0)
        {
        if (min == 0)
            {
            if (sec == 0 )
                std::cout << "Total time: 0.0 sec" << std::endl;
            else
                stream << std::fixed << std::setprecision(1)<<sec;
                std::cout << "Total time: " << stream.str() << " secs" << std::endl;
            }
        else
            {
            if (sec == 0)
                {
                std::cout << "Total time: " << min << " mins " << "0.0" << " sec" << std::endl;  
                }    
            else
                {
                stream << std::fixed << std::setprecision(1) << sec;
                std::cout << "Total time: " << min << " mins " << stream.str() << " sec" << std::endl;
                }
            }
        }
    else
        {
        if (sec == 0)
            std::cout << "Total time: " << hour << " hrs " << min << " mins " << "0.0" << " secs" << std::endl;
        else
            {
            stream << std::fixed << std::setprecision(1) << sec;
            std::cout << "Total time: " << hour << " hrs " << min << " mins " << stream.str() << " secs" << std::endl;
            }
        }
        
}

void print_current_to_current(RoadMap road_map, int from, int to)
{
std::cout << "  Continue to " << road_map.vertexInfo(to) << " (0.0 mile)" << std::endl;
std::cout << "Total distance: 0.0 mile" << std::endl;
}

void print_choose_distance(RoadMap road_map, int start, int end)
{

    std::map<int, int> shortest_map;
    std::stack<std::pair<int, int>> trace_back;
    shortest_map = road_map.findShortestPaths(start, get_distance );
    std::cout << "Shortest distance from " << road_map.vertexInfo(start) << " to " << road_map.vertexInfo(end) << std::endl; 
    std::cout << "  Begin at " << road_map.vertexInfo(start) << std::endl;
    if (end == start)
        {
        print_current_to_current(road_map, start, end);   
        std::cout << "\n\n\n" << std::endl;
        }
    else
        {
        int p = shortest_map[end];
        trace_back.push(std::make_pair(p,end));
        double total_mile = 0.0;
        if (p != start)
            {
            while (p != start)
                {
                int temp = p;
                p = shortest_map[p];
                trace_back.push(std::make_pair(p, temp));
                }
            }
        while (!trace_back.empty())
            {
            double mile = road_map.edgeInfo(trace_back.top().first, trace_back.top().second).miles;
            std::cout << "  Continue to " << road_map.vertexInfo(trace_back.top().second) << " (" << mile << " miles)"
                << std::endl;
            total_mile += mile;
            trace_back.pop();
            }
        std::cout << "Total time: " << total_mile << " miles" << "\n\n\n" << std::endl; 
        }
}

void print_choose_time(RoadMap road_map, int start, int end)
{
    std::map<int, int> shortest_map;
    std::stack<std::pair<int, int>> trace_back;
    shortest_map = road_map.findShortestPaths(start, get_time );
    std::cout << "Shortest driving time from " << road_map.vertexInfo(start) << " to " << road_map.vertexInfo(end) << std::endl; 
    std::cout << "  Begin at " << road_map.vertexInfo(start) << std::endl;
    if (start == end)
        {
        std::cout << "  Continue to " << road_map.vertexInfo(end) <<
        "(0.0 mile @ 0.0mph = 0.0 sec)" << std::endl;
        std::cout << "Total time: 0.0 sec" << "\n\n" << std::endl;
        }
    else
        {
        int p = shortest_map[end];
        trace_back.push(std::make_pair(p, end));
        double total_time = 0.0;
        if ( p != start )
            {
            while (p != start)
                {
                int temp = p;
                p = shortest_map[p];
                trace_back.push(std::make_pair(p, temp));
                }
            }
        while (!trace_back.empty())
            {
            double mile = road_map.edgeInfo(trace_back.top().first, trace_back.top().second).miles;
            double mph = road_map.edgeInfo(trace_back.top().first, trace_back.top().second).milesPerHour;
            std::cout  << "  Continue to " << road_map.vertexInfo(trace_back.top().second) << " ("
                        <<   mile  << " miles"<< " @ " << mph  << "mph "<< " = ";
            double single_time = print_and_return_time(mile, mph);
            total_time += single_time;
            trace_back.pop();
            }
        time_conversion(total_time);
        std::cout << "\n\n"   << std::endl;
        }
}

void get_client_output()
{
    InputReader read_input{std::cin};   
    RoadMapReader read_map{};
    RoadMap road_map = read_map.readRoadMap(read_input);
    TripReader read_trip{};
    std::vector<Trip> target_vector = read_trip.readTrips(read_input);
    if (road_map.isStronglyConnected() == false)
        {
        std::cout << "Disconnected Map" << std::endl;
        return;
        }
    else    
        {
        for (int i = 0; i < target_vector.size(); ++i)
            {
            if (target_vector[i].metric == TripMetric::Distance)//the TripMetric is a class, so we use ::, not .
                {
                    print_choose_distance(road_map, target_vector[i].startVertex,  target_vector[i].endVertex);
                }
            else
                {
                print_choose_time(road_map,target_vector[i].startVertex, target_vector[i].endVertex);
                }
            }
           
        }
}

int main()
{
    get_client_output();

    return 0;
}

