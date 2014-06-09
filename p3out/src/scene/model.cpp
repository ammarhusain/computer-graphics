/**
 * @file model.cpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#include "scene/model.hpp"
#include "scene/material.hpp"
#include "application/opengl.hpp"
#include "scene/triangle.hpp"
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>


namespace _462 {

Model::Model() : mesh( 0 ), material( 0 ) { }
Model::~Model() { }

void Model::render() const
{
    if ( !mesh )
        return;
    if ( material )
        material->set_gl_state();
    mesh->render();
    if ( material )
        material->reset_gl_state();
}

/** - Ammar Husain -**/
Intersection* Model::hasHit( Ray& r )
{

    // transform ray for instancing
    Vector4 e_tmp(r.e.x, r.e.y, r.e.z, 1);
    /// invMat is the transformation matrix stored in parent
    Vector4 t_r_e = invMat*e_tmp;
    Vector4 d_tmp(r.d.x, r.d.y, r.d.z, 0);
    Vector4 t_r_d = invMat*d_tmp;

    Ray t_r( Vector3( t_r_e.x, t_r_e.y, t_r_e.z ),
             Vector3( t_r_d.x, t_r_d.y, t_r_d.z ) );

    // instantiate an intersection container
    ModelIntersection* closestMesh = new ModelIntersection;

    // iterate through all the triangles in the mesh                     
    for (unsigned int tri_ctr = 0;
         tri_ctr < mesh->num_triangles();
         tri_ctr++)
    {
        MeshTriangle currTri = mesh->triangles[ tri_ctr ];
        // pull out triangle vertices   
        MeshVertex v_a = mesh->vertices[ currTri.vertices[ 0 ] ];
        MeshVertex v_b = mesh->vertices[ currTri.vertices[ 1 ] ];
        MeshVertex v_c = mesh->vertices[ currTri.vertices[ 2 ] ];

        // set up cramers rule solution to system of equation
        // refer the Shirley Text - 4.4.2 (2nd ed
        // beta column
        double a = v_a.position.x - v_b.position.x;
        double b = v_a.position.y - v_b.position.y;
        double c = v_a.position.z - v_b.position.z;
        // gamma column
        double d = v_a.position.x - v_c.position.x;
        double e = v_a.position.y - v_c.position.y;
        double f = v_a.position.z - v_c.position.z;
        // t column
        double g = t_r.d.x;
        double h = t_r.d.y;
        double i = t_r.d.z;
        // b column
        double j = v_a.position.x - t_r.e.x;
        double k = v_a.position.y - t_r.e.y;
        double l = v_a.position.z - t_r.e.z;

        // common terms
        real_t ei_minus_hf = (e*i) - (h*f);
        real_t gf_minus_di = (g*f) - (d*i);
        real_t dh_minus_eg = (d*h) - (e*g);
        real_t ak_minus_jb = (a*k) - (j*b);
        real_t jc_minus_al = (j*c) - (a*l);
        real_t bl_minus_kc = (b*l) - (k*c);
        // determinant
        real_t det = (a*ei_minus_hf) + (b*gf_minus_di) + (c*dh_minus_eg);
        real_t t =
            -((f*ak_minus_jb) + (e*jc_minus_al) + (d*bl_minus_kc)) / det;
        // check for a sane t
        if ( (t > closestMesh->t) || (t < closestMesh->epsilon_t) )
            continue;
	    
        // proceed to calculate gamma
        real_t gamma =
            ((i*ak_minus_jb) + (h*jc_minus_al) + (g*bl_minus_kc)) / det;
        // check for a sane gamma
        if ( (gamma < 0.0) || (gamma > 1.0) )
            continue;
	
        // proceed to calculate beta
        real_t beta =
            ((j*ei_minus_hf) + (k*gf_minus_di) + (l*dh_minus_eg)) / det;
        // check for a sane beta
        if ( (beta < 0.0) || (beta + gamma > 1) )
            continue;
	
        // update this intersection as closest
        closestMesh->t = t;
        closestMesh->beta = beta;
        closestMesh->gamma = gamma;
        closestMesh->triangle_id = tri_ctr;
    
    }
    
    closestMesh->ray = r;
    closestMesh->instanced_ray = t_r;
    
    //    std::cout << closest->t << "\t";
    return static_cast<Intersection*>( closestMesh );

}


void Model::populateHit( Intersection* hit )
{

    if (hit == NULL) {
        std::cerr <<
            "Invalid Intersection object passed! WTF!!"
                  << std::endl;
        assert(0);
    }

    // cast the incoming pointer to a local intersection object
    ModelIntersection* thisHit =
        static_cast<ModelIntersection*>(hit);
    // compute alpha
    real_t alpha = 1.0 - (thisHit->beta + thisHit->gamma);

    if (mesh == NULL)
        std::cout << "MESH IS NULL!!!" << std::endl;
    
    MeshTriangle tri = mesh->triangles[ thisHit->triangle_id ];

    MeshVertex v_a = mesh->vertices[ tri.vertices[ 0 ] ];
    MeshVertex v_b = mesh->vertices[ tri.vertices[ 1 ] ];
    MeshVertex v_c = mesh->vertices[ tri.vertices[ 2 ] ];

    hit->int_point.position =
        hit->ray.e + (hit->ray.d*hit->t);

    Vector3 localNormal = (alpha*v_a.normal) + 
        (thisHit->beta*v_b.normal) + 
        (thisHit->gamma*v_c.normal);
    
    Matrix4 normalMatrix;
    transpose( &normalMatrix, invMat ); 
    hit->int_point.normal =
        normalize(multiplyVector(normalMatrix,localNormal));

    /// compute the texture coordinate
    hit->int_point.tex_coord = (alpha*v_a.tex_coord) + 
        (thisHit->beta*v_b.tex_coord) + 
        (thisHit->gamma*v_c.tex_coord);
    
    /// store the material details
    hit->int_material.diffuse = material->diffuse;
    
    hit->int_material.ambient = material->ambient;
    
    hit->int_material.specular = material->specular;
    
    hit->int_material.refractive_index = material->refractive_index;

    int width, height;
    int pix_x, pix_y;
    material->get_texture_size(&width, &height);
    pix_x = (int) fmod(width*hit->int_point.tex_coord.x, width);
    pix_y = (int) fmod(height*hit->int_point.tex_coord.y, height);
    
    hit->int_material.texture =
        material->get_texture_pixel(pix_x, pix_y);
/*
    std::cout << "P: " << pix_x << "," << pix_y << " \t" <<width << "," << height << std::endl;
    
    std::cout << hit->int_point.tex_coord.x << "," << hit->int_point.tex_coord.y << " mesh color: " << hit->int_material.texture << std::endl; 
*/
    
    return;

}

} /* _462 */
