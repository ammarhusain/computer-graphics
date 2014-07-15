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
    bool m_nodeFlag;

    /// default constructor
    Photon() { m_nodeFlag = false; }

    
    /// returns number of dimensions for data point
    static uint numberDimensions () { /// -------------!!!!!!!!!    HACK :(   !!!!!!!!!------------- /// return 3; }
        return 2; }
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

    bool isNode() { return m_nodeFlag; }
    
    
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
     * Tree Constructor: Uses an array to hold the nodes
     * Rounds up the number of nodes passed in to the closest 2^n -1 size
     * container. Since we control the number of photons, make sure they are
     * close to 2^n -1 for an integer n. This optimizes the container and
     * saves space since creating pointers can be costly
     * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
     * @date 07/03/2014 
     * @param nodes 
     ---------------------------------------------------------------------- */
    void create_tree (std::vector<Node> * nodes)
    {

        std::cout << "nodes ptr: " << nodes << std::endl;
        std::cout << "nodes size: " << nodes->size() << std::endl;

        /// compute the size of tree to instantiate
        /// there are 2^n - 1 nodes in a tree
        double log2_size = log2(nodes->size());
        uint n = (unsigned int) std::ceil(log2_size);

        /// check whether the nodes size was not a power of 2
        /// in that case we have to step up one larger; ex 8 nodes
        /// n must equal 4 to since n=3 only stores 7 nodes
        if (fabs((double) n - log2_size) == 0.0)
            n++;
        
        /// instantiate the nodeTree
        /// _nodeTree.resize(nodes->size())
        _nodeTree.resize(pow(2,n) - 1);


        std::cout << "_nodeTree size: " << _nodeTree.size() << std::endl;
        
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
        point.x = 5.5;
        point.y = 5.5;
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

        /// check if this index is a valid node
        if (!self.isNode())
            return;
        
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

        /// initialize variable to stored the farthest photon so far
        std::pair<Node, real_t> farthestCandidate;
        
        /// check if the opposite side needs to be recursed
        /// recurse the other way if: (1) there is still room, (2) distance
        /// of point to splitting plane is closer than the farthest photon
        /// if there is space add yourself
        if (list.size() < k) {
            /// recurse the other way
            recursiveSearch(point, k, list, dimension+1, otherIndex);
        } else {
            /// fetch the farthest photon 
            /// the list should always be sorted from closest to farthest
            farthestCandidate = list[list.size()-1];

            /// check if you are closer than the farthest node
            if (farthestCandidate.second >
                std::abs(splittingValue - point[dimension])) {
                /// remove it and add self instead
                list.pop_back();
                /// recurse the other way
                recursiveSearch(point, k, list, dimension+1, otherIndex);
            }
        }

        /// if there is space add yourself
        if (list.size() < k) {
            add_to_list(list, candidate);
        } else {
            /// fetch the farthest photon 
            /// the list should always be sorted from closest to farthest
            farthestCandidate = list[list.size()-1];

            /// check if you are closer than the farthest node
            if (farthestCandidate.second > my_distance) {                    
                /// remove it and add self instead
                list.pop_back();
                add_to_list(list, candidate);
            }
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
            std::cout << "hitting base" << std::endl;
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

        std::cout << "left: " << leftNodes->size()
                  << "\t right: " << rightNodes->size()
                  << std::endl;
        
        
        /// create the left subtree
        recursiveTree (leftNodes, dimension + 1, (2*index)+1);

        /// create the right subtree
        recursiveTree (rightNodes, dimension + 1, (2*index)+2);
    }    


    
};

} /* _462 */


#endif
