/**
 * @file raytacer.hpp
 * @brief Raytracer class
 *
 * Implement these functions for project 3.
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#ifndef _462_RAYTRACER_HPP_
#define _462_RAYTRACER_HPP_

#define MAX_DEPTH 5

#include "math/color.hpp"
#include "math/random462.hpp"
//#include "scene/scene.hpp"

namespace _462 {

    class Scene;
    class Ray;
    class Mesh;
    struct Intersection;

/*
    // struct to store intersection parameters                
    struct TriangleIntersection
    {
	real_t beta;
	real_t gamma;
	real_t t;
	real_t epsilon_t;

	TriangleIntersection(): t(1000.0), epsilon_t(0.001) {}
    };
*/

    class Raytracer
    {

    public:

	Raytracer();

	~Raytracer();

	bool initialize(Scene* scene, size_t num_samples,
			size_t width, size_t height);

	bool raytrace(unsigned char* buffer, real_t* max_time);

    private:

	Color3 trace_pixel(const Scene* scene,
			   size_t x,
			   size_t y,
			   size_t width,
			   size_t height);

	// the scene to trace
	Scene* scene;

	// the dimensions of the image to trace
	size_t width, height;

	// the next row to raytrace
	size_t current_row;

	unsigned int num_samples;

	/* helper functions */
	Color3 ShootRay( const Scene* scene, Ray& r );

	Color3 shader( const Scene* scene, Intersection* closestGeomIntersection );

	int IntersectGeometry( const Scene* scene, Ray& r,
			       Intersection*& closestGeomIntersection );

    };

} /* _462 */

#endif /* _462_RAYTRACER_HPP_ */
