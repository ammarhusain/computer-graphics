/**
 * @file model.hpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 *         Ammar Husain (ammarh)
 */

#ifndef _462_SCENE_MODEL_HPP_
#define _462_SCENE_MODEL_HPP_

#include "scene/scene.hpp"
#include "scene/mesh.hpp"

namespace _462 {

/**
 * Class to store intersection parameters
 */
class ModelIntersection : public Intersection
{
  public: 
    real_t beta;
    real_t gamma;
    unsigned int triangle_id;

    // default constructor   
    ModelIntersection() : Intersection() { triangle_id = -1; }

};

/**
 * A mesh of triangles.
 */
class Model : public Geometry
{
  public:

    const Mesh* mesh;
    const Material* material;

    Model();
    virtual ~Model();

    virtual void render() const;

    // abstract function to compute hit parameters of a given ray
    virtual Intersection* hasHit( Ray& r );
    virtual void populateHit( Intersection* hit );

};


} /* _462 */

#endif /* _462_SCENE_MODEL_HPP_ */

