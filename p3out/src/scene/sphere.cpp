/**
 * @file sphere.cpp
 * @brief Function defnitions for the Sphere class.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 */

#include "scene/sphere.hpp"
#include "application/opengl.hpp"

namespace _462 {

#define SPHERE_NUM_LAT 80
#define SPHERE_NUM_LON 100

#define SPHERE_NUM_VERTICES ( ( SPHERE_NUM_LAT + 1 ) * ( SPHERE_NUM_LON + 1 ) )
#define SPHERE_NUM_INDICES ( 6 * SPHERE_NUM_LAT * SPHERE_NUM_LON )
// index of the x,y sphere where x is lat and y is lon
#define SINDEX(x,y) ((x) * (SPHERE_NUM_LON + 1) + (y))
#define VERTEX_SIZE 8
#define TCOORD_OFFSET 0
#define NORMAL_OFFSET 2
#define VERTEX_OFFSET 5

#define PI 3.14159265

static unsigned int Indices[SPHERE_NUM_INDICES];
static float Vertices[VERTEX_SIZE * SPHERE_NUM_VERTICES];

static void init_sphere()
{
    static bool initialized = false;
    if ( initialized )
        return;

    for ( int i = 0; i <= SPHERE_NUM_LAT; i++ ) {
        for ( int j = 0; j <= SPHERE_NUM_LON; j++ ) {
            real_t lat = real_t( i ) / SPHERE_NUM_LAT;
            real_t lon = real_t( j ) / SPHERE_NUM_LON;
            float* vptr = &Vertices[VERTEX_SIZE * SINDEX(i,j)];

            vptr[TCOORD_OFFSET + 0] = lon;
            vptr[TCOORD_OFFSET + 1] = 1-lat;

            lat *= PI;
            lon *= 2 * PI;
            real_t sinlat = sin( lat );

            vptr[NORMAL_OFFSET + 0] = vptr[VERTEX_OFFSET + 0] = sinlat * sin( lon );
            vptr[NORMAL_OFFSET + 1] = vptr[VERTEX_OFFSET + 1] = cos( lat ),
            vptr[NORMAL_OFFSET + 2] = vptr[VERTEX_OFFSET + 2] = sinlat * cos( lon );
        }
    }

    for ( int i = 0; i < SPHERE_NUM_LAT; i++ ) {
        for ( int j = 0; j < SPHERE_NUM_LON; j++ ) {
            unsigned int* iptr = &Indices[6 * ( SPHERE_NUM_LON * i + j )];

            unsigned int i00 = SINDEX(i,  j  );
            unsigned int i10 = SINDEX(i+1,j  );
            unsigned int i11 = SINDEX(i+1,j+1);
            unsigned int i01 = SINDEX(i,  j+1);

            iptr[0] = i00;
            iptr[1] = i10;
            iptr[2] = i11;
            iptr[3] = i11;
            iptr[4] = i01;
            iptr[5] = i00;
        }
    }

    initialized = true;
}

Sphere::Sphere()
    : radius(0), material(0) {}

Sphere::~Sphere() {}

void Sphere::render() const
{
    // create geometry if we haven't already
    init_sphere();

    if ( material )
        material->set_gl_state();

    // just scale by radius and draw unit sphere
    glPushMatrix();
    glScaled( radius, radius, radius );
    glInterleavedArrays( GL_T2F_N3F_V3F, VERTEX_SIZE * sizeof Vertices[0], Vertices );
    glDrawElements( GL_TRIANGLES, SPHERE_NUM_INDICES, GL_UNSIGNED_INT, Indices );
    glPopMatrix();

    if ( material )
        material->reset_gl_state();
}

/* - Ammar Husain - */
Intersection* Sphere::hasHit( Ray& r )
{

    Intersection* interParams = new Intersection();
  
    // transform ray for instancing
    Vector4 e_tmp(r.e.x, r.e.y, r.e.z, 1);
    /// invMat is the transformation matrix stored in parent
    Vector4 t_r_e = invMat*e_tmp;
    Vector4 d_tmp(r.d.x, r.d.y, r.d.z, 0);
    Vector4 t_r_d = invMat*d_tmp;

    /// following the naming convention from Shirley text
    Vector3 e = Vector3( t_r_e.x, t_r_e.y, t_r_e.z );
    Vector3 d = Vector3( t_r_d.x, t_r_d.y, t_r_d.z );
    /// for instancing the sphere is placed at (0,0,0)
    Vector3 c = Vector3(0,0,0);

    Vector3 e_sub_c = e-c;

    
    /// compute the discriminant
    real_t discriminant =
        (pow(dot(d, e_sub_c), 2)) -
        (dot(d, d)*(dot(e_sub_c, e_sub_c) - (pow(radius, 2))));

    /// declare variables before the jump labels
    real_t t_1, t_2, t, sqrtDiscriminant;
    
    /// check if there is a real solution
    if (discriminant < 0){
        /// ray does not hit the sphere
        goto final;
    }
    
    sqrtDiscriminant = sqrt(discriminant);

    t_1 = (-dot(d,e_sub_c) + sqrtDiscriminant)/ dot(d,d);
    t_2 = (-dot(d,e_sub_c) - sqrtDiscriminant)/ dot(d,d);

    /// pick the smaller of the two t's
    /// thats the one entering the sphere; the other is exiting
    t = t_1;
    if (t_2 < t_1) 
        t = t_2;

    ///  check for a sane t
    if ((t > interParams->t) || (t < interParams->epsilon_t))
        goto final;

    /// update the intersection structure
    interParams->t = t;

    /// store the uninstanced ray
    interParams->ray = r;
    interParams->instanced_ray = Ray(e, d);
    
  final: return interParams;    
    
}

void Sphere::populateHit( Intersection* hit )
{
    /// sanity check
    if (hit == NULL)
    {
        std::cerr << "Invalid Intersection object passed! WTF!"
                  << std::endl;

        assert(0);
    }

    /// compute position
    hit->int_point.position = hit->ray.e + (hit->ray.d*hit->t);
    /// Shirley Text: unit N = (p-c)/R
    hit->int_point.normal =
        normalize((hit->int_point.position - position)/radius);

    /// -------------!!!!!!!!!    HACK :(   !!!!!!!!!------------- ///
    Vector3 localNormal =
        hit->instanced_ray.e + (hit->instanced_ray.d*hit->t);
    Matrix4 normalMatrix;
    transpose( &normalMatrix, invMat );

    hit->int_point.normal =
        normalize(multiplyVector(normalMatrix,localNormal));
      

    /// compute texture coordinate on sphere
    hit->int_point.tex_coord =
        ComputeSphereTextureCoord(hit->int_point.position);
    
    
    /// populate the material properties
    hit->int_material.ambient = material->ambient;
    hit->int_material.diffuse = material->diffuse;
    hit->int_material.specular = material->specular;
    hit->int_material.refractive_index = material->refractive_index;

    int width, height;
    int pix_x, pix_y;
    material->get_texture_size(&width, &height);
    pix_x = (int) fmod(width*hit->int_point.tex_coord.x, width);
    pix_y = (int) fmod(height*hit->int_point.tex_coord.y, height);
     
    hit->int_material.texture = material->get_texture_pixel(pix_x, pix_y);
    
    return;
}


/** ----------------------------------------------------------------------
 * Computes the texture coordinate for point on sphere
 * Mathematical Details: Shirley 11.2 (2D Texture Mapping)
 * @param hitPosition 
 * 
 * @return 
 ---------------------------------------------------------------------- */
Vector2 Sphere::ComputeSphereTextureCoord(Vector3 hitPosition)
{
    
    Vector2 tex_coord;
    real_t theta = acos((hitPosition.z - position.z)/radius);
    real_t phi =
        atan2(hitPosition.y - position.y, hitPosition.x - position.x);
    if (phi < 0)
        phi = phi + (2*PI);
    
    tex_coord.x = phi/(2*PI);
    tex_coord.y = (PI - theta)/PI;
    return tex_coord;
    
}

} /* _462 */

