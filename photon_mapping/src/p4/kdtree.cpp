/*! -------------------------------------------------------------- 
 * @file   kdtree.cpp
 *
 * @author Ammar Husain (ahusain@nrec.ri.cmu.edu)
 * @date   06/28/2014
 * ---------------------------------------------------------------- 
 */

#include "kdtree.h"

namespace _462 
{

/** ----------------------------------------------------------------------
 * Sorts the data vector passed as argument along the dimension specified
 * 
 * @param data 
 * @param dimension 
 ---------------------------------------------------------------------- */
void Photon::sort (std::vector<Photon> * data,
                   unsigned int dimension)
{

    if (dimension == 0)
        std::sort (data->begin(), data->end(), Photon::comparator_x);
    else if (dimension == 1)
        std::sort (data->begin(), data->end(), Photon::comparator_y);
    else if (dimension == 2)
        std::sort (data->begin(), data->end(), Photon::comparator_z);
    else {
        std::cerr << "Non-existent dimension: " << dimension << std::endl;
        assert(0);
    }
    
}


/** ----------------------------------------------------------------------
 * Returns the normal distance to splitting plane along the
 * dimension specified (splitting dimension)
 * @author: Ammar Husain <ahusain@nrec.ri.cmu.edu>
 * @date 07/03/2014 
 * @param point 
 * @param dimension 
 * 
 * @return 
 ---------------------------------------------------------------------- */
real_t Photon::distance_normal(const Vector3& point, uint dimension)
{
    /// sanity check
    if (dimension >=3) {
        std::cerr << "Non-existent dimension: " << dimension << std::endl;
        assert(0);
    }

    return std::abs(point[dimension] - this->m_position[dimension]);
}


/** ----------------------------------------------------------------------
 * 
 * 
 * @author: Ammar Husain <ahusain@nrec.ri.cmu.edu>
 * @date 07/04/2014 
 * @param point 
 * 
 * @return 
 ---------------------------------------------------------------------- */
real_t Photon::euclidean_squared(const Vector3& point)
{
    return squared_length(point - this->m_position);
}


/** ----------------------------------------------------------------------
 * 
 * 
 * @author: Ammar Husain <ahusain@nrec.ri.cmu.edu>
 * @date 07/04/2014 
 * @param dimension 
 * 
 * @return 
 ---------------------------------------------------------------------- */
real_t Photon::position(uint dimension)
{
    return this->m_position[dimension];
}


/** ----------------------------------------------------------------------
 *
 *
 * @author Ammar Husain <ahusain@nrec.ri.cmu.edu>
 * @date 07/04/2014
 * @param out
 * @param p
 *
 * @return
 ---------------------------------------------------------------------- */
std::ostream& operator<<(std::ostream &out, const Photon &p)
{
    if (p.m_nodeFlag)
        out << "Position: " << p.m_position ;
    else
        out << "Empty Node";

    return out;
}


/** ----------------------------------------------------------------------
 * 
 * 
 * @author: Ammar Husain <ahusain@nrec.ri.cmu.edu>
 * @date 07/04/2014 
 * @param out 
 * @param p 
 * 
 * @return 
 ---------------------------------------------------------------------- */
std::ostream& operator<<(std::ostream &out,
                                const std::vector<Photon> &p)
{

    for (uint i = 0; i < p.size(); i++)
        std::cout << "index: " << i << "\t"
                  << p[i] << "\n --------" << std::endl;

    return out;
}
        

} /* _462 */

