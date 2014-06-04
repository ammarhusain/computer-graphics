/*! -------------------------------------------------------------- 
 * @file   triangle.cpp
 *
 * @author Ammar Husain (ahusain@nrec.ri.cmu.edu)
 * @date   06/01/2014
 * ---------------------------------------------------------------- 
 */
/**
 * @file triangle.cpp
 * @brief Function definitions for the Triangle class.
 *
 * @author Eric Butler (edbutler)
 */

#include "scene/triangle.hpp"
#include "application/opengl.hpp"

namespace _462 {

Triangle::Triangle()
{
	vertices[0].material = 0;
	vertices[1].material = 0;
	vertices[2].material = 0;
}

Triangle::~Triangle() { }

void Triangle::render() const
{
	bool materials_nonnull = true;
	for ( int i = 0; i < 3; ++i )
	    materials_nonnull = materials_nonnull && vertices[i].material;

	// this doesn't interpolate materials. Ah well.
	if ( materials_nonnull )
	    vertices[0].material->set_gl_state();

	glBegin(GL_TRIANGLES);

	glNormal3dv( &vertices[0].normal.x );
	glTexCoord2dv( &vertices[0].tex_coord.x );
	glVertex3dv( &vertices[0].position.x );

	glNormal3dv( &vertices[1].normal.x );
	glTexCoord2dv( &vertices[1].tex_coord.x );
	glVertex3dv( &vertices[1].position.x);

	glNormal3dv( &vertices[2].normal.x );
	glTexCoord2dv( &vertices[2].tex_coord.x );
	glVertex3dv( &vertices[2].position.x);

	glEnd();

	if ( materials_nonnull )
	    vertices[0].material->reset_gl_state();
    }

    /* - Ammar Husain - */
Intersection* Triangle::hasHit( Ray& r )
{
	
	/// instantiate an intersection container
	TriangleIntersection* interParams = new TriangleIntersection();
	// transform ray for instancing
	Vector4 e_tmp(r.e.x, r.e.y, r.e.z, 1);

	/// invMat is the transformation matrix stored in parent
	Vector4 t_r_e = invMat*e_tmp;

	Vector4 d_tmp(r.d.x, r.d.y, r.d.z, 0);

	Vector4 t_r_d = invMat*d_tmp;

	Ray t_r( Vector3( t_r_e.x, t_r_e.y, t_r_e.z ),
             Vector3( t_r_d.x, t_r_d.y, t_r_d.z ) );
	  
	/// pull the individual vertices
	Vertex v_a = vertices[0];
	Vertex v_b = vertices[1];
	Vertex v_c = vertices[2];

	// set up cramers rule solution to system of equation
	// refer the Shirley Text - 10.3.2 (2nd ed
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
	real_t det =
	    (a*ei_minus_hf) + (b*gf_minus_di) + (c*dh_minus_eg);

	real_t t =
	    -((f*ak_minus_jb) + (e*jc_minus_al) + (d*bl_minus_kc)) / det;
	
	/// declare beta and gamma
	real_t beta, gamma;
	    
	// check for a sane t
	if ( (t > interParams->t) || (t < interParams->epsilon_t) )
	    goto final;

	// proceed to calculate gamma
	gamma =
	    ((i*ak_minus_jb) + (h*jc_minus_al) + (g*bl_minus_kc)) / det;

	// check for a sane gamma
	if ( (gamma < 0.0) || (gamma > 1.0) )
	    goto final;
	
	// proceed to calculate beta
	beta =
	    ((j*ei_minus_hf) + (k*gf_minus_di) + (l*dh_minus_eg)) / det;

	// check for a sane beta
	if ( (beta < 0.0) || (beta + gamma > 1) )
	    goto final;
	
	
	// update this intersection as closest
	interParams->t = t;
	interParams->beta = beta;
	interParams->gamma = gamma;
    /// save the untransformed ray
    interParams->ray = r;
    
	goto final;

	/// cast the derived object within base object
  final: return static_cast<Intersection*>(interParams);
}

void Triangle::populateHit( Intersection* hit )
{

    /// sanity check
	if (hit == NULL)
	{
	    
	    std::cerr << "Invalid Intersection object passed! WTF!!"
                  << std::endl;
	    assert(0);
	}

	/// cast the intersection object in the local type
	TriangleIntersection* thisHit =
	    static_cast<TriangleIntersection*>(hit);
	/// compute alpha
	real_t alpha = 1.0 - (thisHit->beta + thisHit->gamma);

	/// get the triangle vertices
	Vertex v_a = vertices[0];
	Vertex v_b = vertices[1];
	Vertex v_c = vertices[2];

	/// store the hit coordinate in parent structure
    hit->int_point.position = hit->ray.e + (hit->ray.d*hit->t);
    
	Vector3 localNormal = (alpha*v_a.normal) +
	    (thisHit->beta*v_b.normal) +
	    (thisHit->gamma*v_c.normal);

    //localNormal = normalize(localNormal);

    
	Matrix4 normalMatrix;
	transpose( &normalMatrix, invMat );
	hit->int_point.normal =
        normalize(multiplyVector(normalMatrix,localNormal));

    /// -------------!!!!!!!!!    HACK :(   !!!!!!!!!------------- ///
    //hit->int_point.normal = normalize(localNormal);
    
    
	hit->int_point.tex_coord = (alpha*v_a.tex_coord) +
	    (thisHit->beta*v_b.tex_coord) +
	    (thisHit->gamma*v_c.tex_coord);

	InterpolateMaterials(hit, alpha, thisHit->beta, thisHit->gamma);
	
	return;
}


void
Triangle::InterpolateMaterials (Intersection* hit,
                                real_t alpha,
                                real_t beta,
                                real_t gamma)
{

	/// Get the triangle vertices
	Vertex v_a = vertices[0];
	Vertex v_b = vertices[1];
	Vertex v_c = vertices[2];

	/// interpolate ambient color
	hit->int_material.ambient =
	    (alpha * v_a.material->ambient) +
	    (beta * v_b.material->ambient) +
	    (gamma * v_c.material->ambient);
	/// interpolate diffuse color
	hit->int_material.diffuse =
	    (alpha * v_a.material->diffuse) +
	    (beta * v_b.material->diffuse) +
	    (gamma * v_c.material->diffuse);
	/// interpolate specular color
	hit->int_material.specular =
	    (alpha * v_a.material->specular) +
	    (beta * v_b.material->specular) +
	    (gamma * v_c.material->specular);
	/// interpolate refractive index
	hit->int_material.refractive_index =
	    (alpha * v_a.material->refractive_index) +
	    (beta * v_b.material->refractive_index) +
	    (gamma * v_c.material->refractive_index);

    int width,height;
    int pix_x, pix_y;

    /// texture color for vertex a
    v_a.material->get_texture_size(&width, &height);
    pix_x = (int) fmod(width*hit->int_point.tex_coord.x, width);
    pix_y = (int) fmod(height*hit->int_point.tex_coord.y, height);
    Color3 a_tex =
        v_a.material->get_texture_pixel(pix_x, pix_y);

    /// texture color for vertex b
    v_b.material->get_texture_size(&width, &height);
    pix_x = (int) fmod(width*hit->int_point.tex_coord.x, width);
    pix_y = (int) fmod(height*hit->int_point.tex_coord.y, height);
    Color3 b_tex =
        v_b.material->get_texture_pixel(pix_x, pix_y);

    /// texture color for vertex c
    v_c.material->get_texture_size(&width, &height);
    pix_x = (int) fmod(width*hit->int_point.tex_coord.x, width);
    pix_y = (int) fmod(height*hit->int_point.tex_coord.y, height);
    Color3 c_tex =
        v_c.material->get_texture_pixel(pix_x, pix_y);

    /// interpolate the texture colors
    hit->int_material.texture =
        (alpha * a_tex) +
        (beta * b_tex) +
        (gamma * c_tex);
    
}
    

    
} /* _462 */
