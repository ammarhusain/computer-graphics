/**
 * @file sphere.hpp
 * @brief Class defnition for Sphere.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_SPHERE_HPP_
#define _462_SCENE_SPHERE_HPP_

#include "scene/scene.hpp"

namespace _462 {

/**
 * A sphere, centered on its position with a certain radius.
 */
class Sphere : public Geometry
{
  public:

    real_t radius;
    const Material* material;

    Sphere();
    virtual ~Sphere();
    virtual void render() const;

    // abstract function to compute hit parameters of a given ray
    virtual Intersection* hasHit( Ray& r );
    virtual void populateHit( Intersection* hit );

  private:
    Vector2 ComputeSphereTextureCoord(Vector3 hitPosition);
    
        
};

} /* _462 */

#endif /* _462_SCENE_SPHERE_HPP_ */

