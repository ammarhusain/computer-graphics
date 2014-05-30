/**
 * @file raytacer.cpp
 * @brief Raytracer class
 *
 * Implement these functions for project 4.
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "raytracer.hpp"
#include "scene/scene.hpp"

#include <SDL_timer.h>
#include <iostream>
#include <random>

#ifdef OPENMP // just a defense in case OpenMP is not installed.

#include <omp.h>

#endif
namespace _462 {

// max number of threads OpenMP can use. Change this if you like.
#define MAX_THREADS 1

static const unsigned STEP_SIZE = 8;

Raytracer::Raytracer()
    : scene(0), width(0), height(0) { }

// random real_t in [0, 1)
static inline real_t random()
{
    return real_t(rand())/RAND_MAX;
}

Raytracer::~Raytracer() { }

/**
 * Initializes the raytracer for the given scene. Overrides any previous
 * initializations. May be invoked before a previous raytrace completes.
 * @param scene The scene to raytrace.
 * @param width The width of the image being raytraced.
 * @param height The height of the image being raytraced.
 * @return true on success, false on error. The raytrace will abort if
 *  false is returned.
 */
bool Raytracer::initialize(Scene* scene, size_t num_samples,
			   size_t width, size_t height)
{
    /*
     * omp_set_num_threads sets the maximum number of threads OpenMP will
     * use at once.
     */
#ifdef OPENMP
    omp_set_num_threads(MAX_THREADS);
#endif
    this->scene = scene;
    this->num_samples = num_samples;
    this->width = width;
    this->height = height;

    current_row = 0;

    Ray::init(scene->camera);
    scene->initialize();

    const SphereLight* lights = scene->get_lights();

    // TODO any initialization or precompuation before the trace

    // Num samples = 1, width = 800, height = 600
    return true;
}

/**
 * Performs a raytrace on the given pixel on the current scene.
 * The pixel is relative to the bottom-left corner of the image.
 * @param scene The scene to trace.
 * @param x The x-coordinate of the pixel to trace.
 * @param y The y-coordinate of the pixel to trace.
 * @param width The width of the screen in pixels.
 * @param height The height of the screen in pixels.
 * @return The color of that pixel in the final image.
 */
Color3 Raytracer::trace_pixel(const Scene* scene,
			      size_t x,
			      size_t y,
			      size_t width,
			      size_t height)
{
    assert(x < width);
    assert(y < height);

    real_t dx = real_t(1)/width;
    real_t dy = real_t(1)/height;

    Color3 res = Color3::Black();

    unsigned int iter;
    for (iter = 0; iter < num_samples; iter++)
    {
        // pick a point within the pixel boundaries to fire our
        // ray through.
        real_t i = real_t(2)*(real_t(x)+random())*dx - real_t(1);
        real_t j = real_t(2)*(real_t(y)+random())*dy - real_t(1);

        Ray r = Ray(scene->camera.get_position(), Ray::get_pixel_dir(i, j));

        // TODO return the color of the given pixel
        // you don't have to use this stub function if you prefer to
        // write your own version of Raytracer::raytrace.
    
	res = ShootRay( scene, r );
    }


    return res*(real_t(1)/num_samples);
}

/**
 * Raytraces some portion of the scene. Should raytrace for about
 * max_time duration and then return, even if the raytrace is not copmlete.
 * The results should be placed in the given buffer.
 * @param buffer The buffer into which to place the color data. It is
 *  32-bit RGBA (4 bytes per pixel), in row-major order.
 * @param max_time, If non-null, the maximum suggested time this
 *  function raytrace before returning, in seconds. If null, the raytrace
 *  should run to completion.
 * @return true if the raytrace is complete, false if there is more
 *  work to be done.
 */
bool Raytracer::raytrace(unsigned char* buffer, real_t* max_time)
{
    // TODO Add any modifications to this algorithm, if needed.

    std::cout << "# of -- \n \t meshes = " << scene->num_meshes() 
	      << " \n \t geometries = " << scene->num_geometries()
	      << " \n \t materials = " << scene->num_materials()
 	      << " \n \t lights = " << scene->num_lights()
	      << std::endl;


    static const size_t PRINT_INTERVAL = 64;

    // the time in milliseconds that we should stop
    unsigned int end_time = 0;
    bool is_done;

    if (max_time)
    {
        // convert duration to milliseconds
        unsigned int duration = (unsigned int) (*max_time * 1000);
        end_time = SDL_GetTicks() + duration;
    }

    // until time is up, run the raytrace. we render an entire group of
    // rows at once for simplicity and efficiency.
    for (; !max_time || end_time > SDL_GetTicks(); current_row += STEP_SIZE)
    {
        // we're done if we finish the last row
        is_done = current_row >= height;
        // break if we finish
        if (is_done) break;

        int loop_upper = std::min(current_row + STEP_SIZE, height);

        // This tells OpenMP that this loop can be parallelized.
        #pragma omp parallel for
        for (int c_row = current_row; c_row < loop_upper; c_row++)
        {
            /*
             * This defines a critical region of code that should be
             * executed sequentially.
             */
             #pragma omp critical
            {
                if (c_row % PRINT_INTERVAL == 0)
                    printf("Raytracing (Row %d)\n", c_row);
            }

            for (size_t x = 0; x < width; x++)
            {
                // trace a pixel
                Color3 color = trace_pixel(scene, x, c_row, width, height);
                // write the result to the buffer, always use 1.0 as the alpha
                color.to_array(&buffer[4 * (c_row * width + x)]);
            }
        }
    }

    if (is_done) printf("Done raytracing!\n");

    return is_done;
}

    
/** ----------------------------------------------------------------
 * Recursive function to create 3 rays: shadow, reflective
 * and refractive
 * @param scene 
 * @param r 
 * 
 * @return 
 ---------------------------------------------------------------- */
    Color3 Raytracer::ShootRay( const Scene* scene, Ray& r )
    {

	// iterate through all the geometries to find a hit
	Intersection* closestGeomIntersection = new Intersection();
	// ray intersection test
	int closestGeom_ind = IntersectGeometry(scene,
						r,
						closestGeomIntersection);

	Color3 res = scene->background_color;

	// request geometry for color
	if( closestGeom_ind >= 0 ) {
	    // fetch geometries in the scene
	    Geometry* const* sceneGeometries = scene->get_geometries();

	    // get the details of the ray intersection to geometry
	    sceneGeometries[ closestGeom_ind ]->populateHit( closestGeomIntersection );
    
	    // compute shading color
	    res = shader( scene, closestGeomIntersection );

	}

	// destroy the intersection object
	delete closestGeomIntersection;
	return res;
    }


/** ----------------------------------------------------------------
 * 
 * 
 * @param scene 
 * @param intersection 
 * 
 * @return 
 ---------------------------------------------------------------- */
Color3 Raytracer::shader( const Scene* scene, Intersection* intersection )
{

    // sanity check for an instantiated pointer
    if (intersection == NULL)
        return Color3::Black();
        
    
    // instantitate required diffuse and ambient material/ light colors
    Color3 k_d = intersection->int_material.diffuse;
    Color3 k_a = intersection->int_material.ambient;
    Color3 c_a = scene->ambient_light;
    // get texture color
    Color3 t_p = intersection->int_material.texture;

    // instantiate the normal vector for the intersection
    Vector3 normal = intersection->int_point.normal;
    // instantiate a color accumulator to add light color contributions
    Color3 lightAccumulator(0.0, 0.0, 0.0);

    // fetch all the light sources
    const SphereLight* sceneLights = scene->get_lights();

    // iterate through the light sources
    for ( unsigned int light_ctr = 0; light_ctr < scene->num_lights(); light_ctr++ )
    {
        SphereLight currentLight = sceneLights[ light_ctr ];

        // store color of light
        Color3 c = currentLight.color;
        // store attenuation terms
        real_t a_c = currentLight.attenuation.constant;
        real_t a_l = currentLight.attenuation.linear;
        real_t a_q = currentLight.attenuation.quadratic;

        // monte carlo sampling = 10 samples / light
        for ( unsigned int sample_ctr = 0; sample_ctr < 10; sample_ctr++ )
        {
            Vector3 lightSample;
            lightSample.x = random_gaussian();
            lightSample.y = random_gaussian();
            lightSample.z = random_gaussian();

            // normalize the vector
            lightSample = normalize( lightSample );
            // scale the light sample
            lightSample = currentLight.radius * lightSample;
            //transform the light sample
            lightSample += currentLight.position;

            // instantiate light ray
            Ray L = Ray( intersection->int_point.position, lightSample );

            Intersection* lightIntersection = new Intersection();
            // check for obstruction in path to light
            int intersection_ind = IntersectGeometry( scene,
                                                      L,
                                                      lightIntersection );
            // if the pointer points to a valid intersection container
            if ( intersection_ind < 0 ) {
                // light is not blocked by any geomtery

                // compute distance to light source
                real_t d = length( intersection->int_point.position - lightSample );

                Color3 c_i = c * ( 1 / (a_c + (d*a_l) + (pow(d,2)*a_q)) );
                real_t normalDotLight = dot( normal, L.d );
                if ( normalDotLight > 0 )
                    lightAccumulator += (c_i * k_d * normalDotLight);

            }
            // free the intersection container
            delete lightIntersection;

        }
    }
    
    //std::cout << "texture: " << t_p <<std::endl;
    return t_p*((c_a*k_a) + lightAccumulator);

}


/** ----------------------------------------------------------------
 * 
 * 
 * @param scene 
 * @param r 
 * @param closestGeomIntersection 
 * 
 * @return 
 ---------------------------------------------------------------- */
int Raytracer::IntersectGeometry( const Scene* scene,
				  Ray& r, 
				  Intersection*& closestGeomIntersection )
{

    // initialize closest index in scene geometries
    int closestGeom_ind = -1;
    Geometry* const* sceneGeometries = scene->get_geometries();

    // initialize input container if not already done so
    if (closestGeomIntersection == NULL)
	closestGeomIntersection = new Intersection();

    for (unsigned int geom_ctr = 0; geom_ctr < scene->num_geometries(); geom_ctr++)
    {
	Intersection* currIntersection = sceneGeometries[ geom_ctr ]->hasHit( r );
	// check if this geometry is closer to ray
	if (currIntersection->t < closestGeomIntersection->t) {
	    delete closestGeomIntersection;
	    closestGeomIntersection = currIntersection;
	    closestGeom_ind = geom_ctr;
	}
	else
	    delete currIntersection;
    }

    return closestGeom_ind;
}

} /* _462 */
