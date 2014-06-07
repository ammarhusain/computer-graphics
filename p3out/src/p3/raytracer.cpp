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
#include <ctime>
/*
#include "scene/scene.hpp"

#include <SDL_timer.h>
#include <iostream>
#include <random>
*/
#ifdef OPENMP // just a defense in case OpenMP is not installed.

#include <omp.h>

#endif
namespace _462 {

// max number of threads OpenMP can use. Change this if you like.
#define MAX_THREADS 1

#define ANTI_ALIASING_SAMPLES 1//8

#define MONTE_CARLO_LIGHT_SAMPLES 1//5

#define PRINT_TIMING 0

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
    /// increase the number of samples for antialiasing
    this->num_samples = ANTI_ALIASING_SAMPLES;
    
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

    /// stack container to keep track of refractive indices
    std::stack<real_t> refractive_indices;
    refractive_indices.push(scene->refractive_index);        

    /// std::cout << "W: " << width << "\t H: " << height << std::endl;
    
    /// std::cout << "Pixel: " << x << " , " << y << std::endl;

    std::clock_t start;

    start = std::clock();
    
    /// !!!!---- TODO: Increase num_samples to get rid of antialiasing
    unsigned int iter;
    for (iter = 0; iter < num_samples; iter++)
    {
        // pick a point within the pixel boundaries to fire our
        // ray through.
        real_t i = real_t(2)*(real_t(x)+random())*dx - real_t(1);
        real_t j = real_t(2)*(real_t(y)+random())*dy - real_t(1);

        Ray r = Ray(scene->camera.get_position(), Ray::get_pixel_dir(i, j));

        /// std::cout << "\t sample- " << iter << "\t";
        
        
        // TODO return the color of the given pixel
        // you don't have to use this stub function if you prefer to
        // write your own version of Raytracer::raytrace.
        res += RecursiveRayTrace(scene, r, 0, refractive_indices);

        /// std::cout << " " << std::endl;
    }

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    if (PRINT_TIMING)
        std::cout << "Pixel: " << x << "," << y << "\t processing time: " << duration << std::endl;
    
    
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
            // This tells OpenMP that this loop can be parallelized.

#pragma omp parallel for

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
Color3 Raytracer::RecursiveRayTrace(const Scene* scene, Ray r, int depth,
                                    std::stack<real_t> refractive_indices)
{


    /// std::cout << depth << "\t";
    
    /// !!!!---- TODO: Need to account for refraction
    
    /// get background color
    Color3 backgroundColor = scene->background_color;

    
	// iterate through all the geometries to find a hit
	Intersection* closestGeomIntersection = new Intersection();

    /// check if ray hits an object in scene
	int closestGeom_ind = RayCast(scene,
                                  r,
                                  closestGeomIntersection);

    Color3 finalColor(0.0, 0.0, 0.0);
    
    /// instantiate the sources of color
    Color3 lightContribution(0.0, 0.0, 0.0),
        recursiveContribution(0.0, 0.0, 0.0);
    

	// request geometry for color
	if( closestGeom_ind >= 0 ) {

        /// shadow rays => visibility of light
        lightContribution =
            SampleShadowRays(scene, closestGeomIntersection);

        /**/
        /// check for recursion termination condition
        if (depth < 5) {
            /// REFLECTION COMPONENT
            /// create a reflected ray
            Vector3 intersectionNormal =
                normalize(closestGeomIntersection->int_point.normal); // dont need to normalize
            Vector3 intersectionPoint =
                closestGeomIntersection->int_point.position;
            Vector3 reflectionVector =
                normalize(r.d - (2*dot(r.d, intersectionNormal)*intersectionNormal));
            Ray reflectedRay(intersectionPoint, reflectionVector);
            
            /// FRESNEL EFFECT
            if(closestGeomIntersection->int_material.refractive_index > 0)
            {
                Vector3 incomingDirection =
                    closestGeomIntersection->ray.d;
                /// normalize incoming ray: just in case
                incomingDirection = normalize(incomingDirection);
                
                real_t n = refractive_indices.top();
                real_t n_t =
                    closestGeomIntersection->int_material.refractive_index;

                /// -------------!!!!!!!!!    HACK :(   !!!!!!!!!------------- ///
                /**/
                /// figure out if entering or exiting the medium
                /// should ideally just be done once: needs to be cleaned out
                if (dot(incomingDirection, intersectionNormal) < 0)
                {}/// nothing to do here  
                else /// exiting dielectric
                {
                    /// swap the n and n_t;
                    real_t tmp;
                    tmp = n;
                    n = n_t;
                    n_t = tmp;
                }
                /**/
                
                /// compute the refracted ray direction: Shirley 10.7
                real_t squareRootTerm = 
                    real_t(1.0) - (pow(n,2)/pow(n_t,2)*(real_t(1.0) - pow(dot(incomingDirection, intersectionNormal),2)));

                /// randomly pick reflection vs refraction
                if (squareRootTerm >= 0)
                {
                    Vector3 outgoingDirection =
                        ((n/n_t)*(incomingDirection - intersectionNormal*dot(incomingDirection, intersectionNormal)))
                        - (intersectionNormal*sqrt(squareRootTerm));

                    /// normalize outgoing ray
                    outgoingDirection = normalize(outgoingDirection);
                    
                    real_t R = getFresnelCoefficient(incomingDirection,
                                                     outgoingDirection,
                                                     intersectionNormal,
                                                     std::make_pair(n,n_t));

                    /// do a coin toss to figure out which ray to emit
                    real_t unbiased_coin = random_gaussian();
                    /// Probability: reflection = R, refraction = 1-R
                    if (unbiased_coin < R) {
                        /// REFLECTION
                        recursiveContribution =
                            RecursiveRayTrace(scene, reflectedRay,
                                              depth+1, refractive_indices);
                        /// multiply this with specular color of material and texture color
                        recursiveContribution =
                            recursiveContribution *
                            closestGeomIntersection->int_material.specular *
                            closestGeomIntersection->int_material.texture;    
                    }
                    else {
                        /// REFRACTION
                        /// entering dielectric
                        if (dot(incomingDirection, intersectionNormal) < 0)
                        {
                            //std::cout << "Entering dielectric" << std::endl;
                            refractive_indices.push(n_t);
                        }
                        
                        else /// exiting dielectric
                        {
                            //std::cout << "Exiting dielectric" << std::endl;
                            refractive_indices.pop();
                        }
                            
                        //std::cout << "Refracting babay" << std::endl;
                        
                        Ray refractedRay(intersectionPoint, outgoingDirection);
                        recursiveContribution =
                            RecursiveRayTrace(scene, refractedRay,
                                              depth+1, refractive_indices);
                        
                    }
                    
                }
                /// Total Internal Reflection
                else {
                    
                    recursiveContribution =
                        RecursiveRayTrace(scene, reflectedRay,
                                          depth+1, refractive_indices);
                    /// multiply this with specular color of material and texture color
                    recursiveContribution =
                        recursiveContribution *
                        closestGeomIntersection->int_material.specular *
                        closestGeomIntersection->int_material.texture;    
                }

            } /// Fresnel effect
            else
            {
                recursiveContribution =
                    RecursiveRayTrace(scene, reflectedRay,
                                      depth+1, refractive_indices);
            /// multiply this with specular color of material and texture color
            recursiveContribution =
                recursiveContribution *
                closestGeomIntersection->int_material.specular *
                closestGeomIntersection->int_material.texture;    

            } /// pure reflection
            
            
        }
        /**/
        /// add up the various colors
        finalColor = lightContribution + recursiveContribution;
	}
    else {
        finalColor = backgroundColor;
    }
    

	// destroy the intersection object
	delete closestGeomIntersection;
	return finalColor;
}


/** ----------------------------------------------------------------------
 * Computes the contribution of light at the point of intersection
 * 
 * @param scene 
 * @param intersection 
 * 
 * @return 
 ---------------------------------------------------------------------- */
Color3
Raytracer::SampleShadowRays(const Scene* scene, Intersection* intersection)
{

    // sanity check for an instantiated pointer
    if (intersection == NULL) {
        std::cerr << "SAMPLE-SHADOW-RAYS: Passed in a NULL Intersection! WTF!"
                  << std::endl;
        return Color3::Black();
    }
    
    // instantitate required diffuse and ambient material/ light colors
    Color3 k_d = intersection->int_material.diffuse;
    Color3 k_a = intersection->int_material.ambient;
    Color3 c_a = scene->ambient_light;
    // get texture color
    Color3 t_p = intersection->int_material.texture;

    // instantiate the normal vector for the intersection
    Vector3 normal = normalize(intersection->int_point.normal);
    /// variable to sum over all light sources
    Color3 avgLightColor(0.0, 0.0, 0.0);

    // fetch all the light sources
    const SphereLight* sceneLights = scene->get_lights();

    int numSamples = MONTE_CARLO_LIGHT_SAMPLES;
    
    // iterate through the light sources
    for (unsigned int light_ctr = 0;
         light_ctr < scene->num_lights(); light_ctr++)
    {
        SphereLight currentLight = sceneLights[light_ctr];

        // store color of light
        Color3 c = currentLight.color;
        // store attenuation terms
        real_t a_c = currentLight.attenuation.constant;
        real_t a_l = currentLight.attenuation.linear;
        real_t a_q = currentLight.attenuation.quadratic;

        // instantiate a color accumulator to add light color contributions
        Color3 lightAccumulator(0.0, 0.0, 0.0);

        // monte carlo sampling = 10 samples / light
        for (unsigned int sample_ctr = 0;
             sample_ctr < numSamples; sample_ctr++)
        {
            Vector3 lightSample;
            lightSample.x = random_gaussian() - 0.50;
            lightSample.y = random_gaussian() - 0.50;
            lightSample.z = random_gaussian() - 0.50;
            lightSample.x = random_gaussian();
            lightSample.y = random_gaussian();
            lightSample.z = random_gaussian();

            // normalize the vector
            lightSample = normalize( lightSample );
            // scale the light sample
            lightSample = (currentLight.radius * lightSample);
            //transform the light sample
            lightSample += currentLight.position;
            
            Vector3 sampleDirection =
                normalize(lightSample - intersection->int_point.position);
            
            //lightSample =
            //    currentLight.position - intersection->int_point.position;
                
            // instantiate light ray
            Ray L = Ray(intersection->int_point.position, sampleDirection);

            /// compute t till the light intersection
            real_t t_light = length(lightSample - intersection->int_point.position);
            
            Intersection* lightIntersection = new Intersection();

            /// set max t that ray is allowed to travel
            lightIntersection->t = t_light;
            
            // check for obstruction in path to light
            int intersection_ind = RayCast( scene,
                                            L,
                                            lightIntersection );
            
            // free the intersection container
            // do not need to know where the light got blocked
            delete lightIntersection;

            // if the pointer points to a valid intersection container
            if (intersection_ind >= 0)
                continue;
                    
            // light is not blocked by any geomtery
            // compute distance to light source
            real_t d = length(intersection->int_point.position - lightSample);

            Color3 c_i = c * (real_t(1.0)/(a_c + (d*a_l) + (pow(d,2.0)*a_q)));

            real_t normalDotLight = dot(normal, L.d);
            if (normalDotLight > 0) {
                /// sum over all samples sent to a light source
                lightAccumulator += (c_i * k_d * normalDotLight);
            }            
        }
        /// sum over all lights in scene
        avgLightColor += lightAccumulator*(real_t(1.0)/real_t(numSamples));
    }

    /// average over the various lights in the scene
    avgLightColor = avgLightColor*(real_t(1.0)/real_t(scene->num_lights()));
    
    Color3 finalLightColor = t_p*((c_a*k_a) + avgLightColor);
    
    return finalLightColor;
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
int Raytracer::RayCast( const Scene* scene,
                        Ray& r, 
                        Intersection*& closestGeomIntersection )
{

    // initialize closest index in scene geometries
    int closestGeom_ind = -1;
    Geometry* const* sceneGeometries = scene->get_geometries();

    // initialize input container if not already done so
    //if (closestGeomIntersection == NULL)
    //    closestGeomIntersection = new Intersection();

    
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

    /// if a valid hit occured => get details of intersection
    if (closestGeom_ind >= 0)
        sceneGeometries[closestGeom_ind]->populateHit(closestGeomIntersection);
    
    return closestGeom_ind;
}


/** ----------------------------------------------------------------------
 * Comoutes the Fresnel coefficient R for reflectivity of dielectric
 * Implementation: Shirley 10.7
 * @param incoming 
 * @param outgoing 
 * @param normal 
 * 
 * @return 
 ---------------------------------------------------------------------- */
real_t Raytracer::getFresnelCoefficient(Vector3 incoming,
                                        Vector3 outgoing, Vector3 normal,
                                        std::pair<real_t, real_t> r_ind)
{

    /// for sanity: normalize everything
    incoming = normalize(incoming);
    outgoing = normalize(outgoing);
    normal = normalize(normal);

    Vector3 incidenceVector;
    real_t cos_theta;
    real_t n_t, n;

    n = r_ind.first;
    n_t = r_ind.second;
    cos_theta = dot(incoming, normal);
    
    
    
    /// figure out the ray with larger incidence angle
    if(dot(incoming, normal) < dot(outgoing, normal))
    {
        cos_theta = dot(incoming, normal);
    } else {
        cos_theta = dot(outgoing, normal);
    }
    

    /// direction does not matter here: take absolute value of angle
    cos_theta = fabs(cos_theta);
    
    real_t R_o = pow(((n_t-n)/(n_t+n)), 2);

    real_t R = R_o + ((1.0-R_o)*pow(1-cos_theta, 5));
/*
    std::cout << "R_o:" << R_o << "\tcos:" << cos_theta << "\tpow:"
              << pow(real_t(1.0-cos_theta), 5) << "\tR:" << R << std::endl;
*/  
    return R;
    
}

} /* _462 */
