#ifndef _KDTREE_H_
#define _KDTREE_H_

/** -------------------------------------------------------------- 
 * @file   kdtree.h
 * 
 * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
 * @date   06/28/2014
 * 
 * @brief  
 * Implementation of the Kd-Tree data structure
 * 
 * ------------------------------------------------------------ */


#include "math/vector.hpp"
#include "math/color.hpp"
#include "scene/ray.hpp"
#include <vector>
#include <cmath>        // std::abs

namespace _462 
{


struct Photon 
{
    /// x, y, z position of the photon
    Vector3 m_position;
    /// color intensity
    Color3 m_intensity;
    /// direction of incidence
    Ray m_incidentRay;
    /// flag to indicate if this Photon is a valid node
    bool nodeFlag;

    /// default constructor
    Photon() { nodeFlag = false; }

    
    /// returns number of dimensions for data point
    static uint numberDimensions () { return 3; }

    /** ----------------------------------------------------------------------
     * static functions to compare photon positions for sorting
     * 
     * @param i 
     * @param j 
     * 
     * @return 
     ---------------------------------------------------------------------- */
    static bool comparator_x (Photon i, Photon j)
    { return (i.m_position.x < j.m_position.x); }

    static bool comparator_y (Photon i, Photon j)
    { return (i.m_position.y < j.m_position.y); }

    static bool comparator_z (Photon i, Photon j)
    { return (i.m_position.z < j.m_position.z); }
    
    /// ---------------------------------------------------------------------- 

    static void sort (std::vector<Photon> * data,
                     uint dimension);

    real_t distance_normal(const Vector3& point, uint dimension);

    real_t euclidean_squared(const Vector3& point);

    real_t position(uint dimension);


    
};



/** ----------------------------------------------------------------------
 * Overloaded ostream operators to aid in printing Photons
 * 
 --------------------------------------------------------------------- */
std::ostream& operator<<(std::ostream &out, const Photon &p);

std::ostream& operator<<(std::ostream &out,
                         const std::vector<Photon> &p);
    
/* ------------------------------------------------------------------- */ 
    
        

/** -------------------------------------------------------------- 
 * @class  KdTree
 * 
 * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
 * @date   06/28/2014
 * 
 * @brief  
 * 
 * ------------------------------------------------------------ */
template <typename Node>
class KdTree
{
    
  public:
    /// default constructor
    KdTree(){}
    KdTree(std::vector<Node> * nodes) { this->create_tree(nodes); }
    
    ~KdTree(){}
    

    /** ----------------------------------------------------------------------
     * Tree Constructor
     * 
     * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
     * @date 07/03/2014 
     * @param nodes 
     ---------------------------------------------------------------------- */
    void create_tree (std::vector<Node> * nodes)
    {

        std::cout << "nodes ptr: " << nodes << std::endl;
        std::cout << "nodes size: " << nodes->size() << std::endl;
        
        /// instantiate the nodeTree to the size of data
        _nodeTree.resize(nodes->size());

        std::vector<Node> * allNodes =
            new std::vector<Node>(nodes->begin(), nodes->end());

        /// destroy the original set of nodes
        nodes->clear();
        
        /// recursively build the tree
        this->recursiveTree (allNodes, 0, 0);

        /// -------------!!!!!!!!!    HACK :(   !!!!!!!!!------------- ///
        /// iterate through and print the nodes in tree
        std::cout << _nodeTree << std::endl;

        Vector3 point;
        point.x = 7;
        point.y = 1;
        point.z = 0;

        std::vector<Node> neighbors = nearestNeighborSearch(point, 2);

        std::cout << "Closest Nodes to Point:  " << point << "\n" << neighbors;
        
        
        
    }



    /** ----------------------------------------------------------------------
     * Computes the k nearest photons from photon map to the input node
     * 
     * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
     * @date 07/03/2014 
     * @param point 
     * @param k 
     * 
     * @return 
     ---------------------------------------------------------------------- */
    std::vector<Node> nearestNeighborSearch(Vector3 point, uint k)
    {

        std::vector< std::pair<Node, real_t> > list;

        recursiveSearch(point, k, list, 0, 0);

        /// allocate a container for closest nodes
        std::vector<Node> neighbors(list.size());

        /// strip out the paired distances and return the nodes
        for (uint i = 0; i < list.size(); i++)
            neighbors[i] = list[i].first;
        
        return neighbors;

    }
    
    

    /// protected:
  private:
    /// kd-tree: vector of nodes
    /// since the tree is supposed to be balanced -> an array can be used
    /// since the median is the larger of 2 nodes -> left node always exists
    /// if not a leaf
    std::vector<Node> _nodeTree;

    /** ----------------------------------------------------------------------
     * 
     * 
     * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
     * @date 07/03/2014 
     * @param point 
     * @param k 
     * @param list 
     * @param dimension 
     ---------------------------------------------------------------------- */
    void recursiveSearch (const Vector3& point,
                          uint k,
                          std::vector< std::pair<Node, real_t> >& list,
                          uint dimension,
                          uint index)
    {

        uint leftIndex = (2*index) + 1;
        uint rightIndex = (2*index) + 2;
        Node self = _nodeTree[index];
        
        /// base condition: check if it is a leaf node
        if (leftIndex >= _nodeTree.size())
        {
            /// when hitting a leaf there should always be room 
            if (list.size() >= k) {
                std::cerr << "List ran out of space on leaf!" << std::endl;
                assert(0);
            }

            real_t distance = self.euclidean_squared(point);
            
            std::pair<Node, real_t> candidate;

            candidate = std::make_pair(self, distance);

            add_to_list(list, candidate);
            
            return;
        }

        /// roll over the dimension
        dimension = dimension % Node::numberDimensions();

        /// initialize indices to move forward
        /// go left unless instructed otherwise
        uint recurseIndex = leftIndex;
        uint otherIndex = rightIndex;

        uint splittingValue = self.position(dimension);
        
        /// check which way to recurse
        if (splittingValue < point[dimension]) {
            /// the point lies to the right of splitting place
            /// make sure that a right node exists
            if (rightIndex < _nodeTree.size())
            {
                /// go right instead
                recurseIndex = rightIndex;
                otherIndex = leftIndex;
            }
                
        }
        
        
        /// go in the recursion direction
        recursiveSearch(point, k, list, dimension+1, recurseIndex);

        /// prepare yourself to be added if needed
        std::pair<Node, real_t> candidate;
        real_t my_distance = self.euclidean_squared(point);
        candidate = std::make_pair(self, my_distance);

        /// fetch the farthest photon so far
        /// the list should always be sorted from closest to farthest
        std::pair<Node, real_t> farthestCandidate = list[list.size()-1];

        /// if there is space add yourself
        if (list.size() < k) {
            add_to_list(list, candidate);
        }
        /// check if you are closer than the farthest node
        else if (farthestCandidate.second > my_distance) {                    
            /// remove it and add self instead
            list.pop_back();
            add_to_list(list, candidate);
        }

        
        /// if opposing side does not exist we are done
        /// since the tree is balanced this will always be the right side
        if (otherIndex >= _nodeTree.size())
            return;
        
        
        /// check if the opposite side needs to be recursed
        /// recurse the other way if: (1) there is still room, (2) distance
        /// of point to splitting plane is closer than the farthest photon

        if (farthestCandidate.second >=
            std::abs(splittingValue - point[dimension]))
        {
            /// make some roon
            list.pop_back();
        }
        /// check if there is room
        if (list.size() < k)
        {
            /// recurse the other way
            recursiveSearch(point, k, list, dimension+1, otherIndex);
        }
        
        
    }
    

    /** ----------------------------------------------------------------------
     * Pushes new element to list and the sorts it according to the 
     * corresponding distance values in pair
     * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
     * @date 07/03/2014 
     * @param list 
     * @param candidate 
     ---------------------------------------------------------------------- */
    void add_to_list(std::vector< std::pair<Node, real_t> >& list,
                     std::pair<Node, real_t> candidate)
    {
        list.push_back(candidate);
        std::sort(list.begin(), list.end(), listCompare);        
    }
    

    /** ----------------------------------------------------------------------
     * Compares the distance values in the photon pair and return if 
     * a is smaller than b
     * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
     * @date 07/03/2014 
     * @param a 
     * @param b 
     * 
     * @return 
     ---------------------------------------------------------------------- */
    static bool listCompare (std::pair<Node, real_t> a, std::pair<Node, real_t> b)
    {
        return (a.second < b.second);
    }
    

    /** ----------------------------------------------------------------------
     * Recursive Function to build the Kd-Tree
     * 
     * @author: Ammar Husain <ahusain@nrec.ri.cmu.edu>
     * @date 07/03/2014 
     * @param nodes 
     * @param dimension 
     * @param index 
     ---------------------------------------------------------------------- */
    void recursiveTree (std::vector<Node> * nodes,
                        uint dimension,
                        uint index)
    {
        

        std::cout << "nodes: " << nodes->size()
                  << " index: " << index
                  << std::endl;

        /// sanity check on whether any nodes were passed in recursion
        if (nodes->size() == 0)
        {
            /// use the default photon with a false flag
            return;
        }    

        /// sanity check on index
        if (index >= _nodeTree.size())
        {
            std::cerr
                << "Error: tree index exceeds # of nodes: index = "
                << index
                << std::endl;
            assert(0);
        }

        /// base condition
        if (nodes->size() == 1)
        {
            _nodeTree[index] = nodes->at(0);
            return;
        }

        Node::sort(nodes, (dimension % Node::numberDimensions()) );

        /// get the median index
        uint median = std::floor(nodes->size()/2);


        /// set the median as the current node
        _nodeTree[index] = nodes->at(median);

        std::vector<Node> * leftNodes =
            new std::vector<Node>(nodes->begin(), nodes->begin() + median);

        std::vector<Node> * rightNodes =
            new std::vector<Node>(nodes->begin() + (median+1), nodes->end());

        /// destroy the input nodes to avoid stack overflow during recursion
        delete nodes;

        /// create the left subtree
        recursiveTree (leftNodes, dimension + 1, (2*index)+1);

        /// create the right subtree
        recursiveTree (rightNodes, dimension + 1, (2*index)+2);
    }    


    
};

} /* _462 */


#endif
