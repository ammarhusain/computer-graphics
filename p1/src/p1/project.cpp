/**
 * @file project.cpp
 * @brief OpenGL project
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "p1/project.hpp"

// use this header to include the OpenGL headers
// DO NOT include gl.h or glu.h directly; it will not compile correctly.
#include "application/opengl.hpp"

// A namespace declaration. All proejct files use this namespace.
// Add this declration (and its closing) to all source/headers you create.
// Note that all #includes should be BEFORE the namespace declaration.
namespace _462 {

// definitions of functions for the OpenglProject class

// constructor, invoked when object is created
OpenglProject::OpenglProject()
{
    // TODO any basic construction or initialization of members
    // Warning: Although members' constructors are automatically called,
    // ints, floats, pointers, and classes with empty contructors all
    // will have uninitialized data!
}

// destructor, invoked when object is destroyed
OpenglProject::~OpenglProject()
{
    // TODO any final cleanup of members
    // Warning: Do not throw exceptions or call virtual functions from deconstructors!
    // They will cause undefined behavior (probably a crash, but perhaps worse).
}

/**
 * Initialize the project, doing any necessary opengl initialization.
 * @param camera An already-initialized camera.
 * @param scene The scene to render.
 * @return true on success, false on error.
 */
bool OpenglProject::initialize( Camera* camera, Scene* scene )
{
  // copy scene
  this->scene = *scene;

  // TODO opengl initialization code and precomputation of mesh/heightmap

  initializePoolMesh();

  // set the camera
  setCamera( camera );


  // add lighting to the scene
  /*
  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[] = { 50.0 };
  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_SMOOTH);

  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
  */
  GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
  GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };

  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
  

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);

  // initialize the watersurface
  createWaterSurfaceMesh();

  return true;
}

/**
 * Function to precompute the normals vectors of the pool and add color
 */
void OpenglProject::initializePoolMesh()
{

  // allocate normal vectors for pool mesh
  scene.mesh.normals = new Vector3 [scene.mesh.num_vertices];

  // compute vertex normals in pool mesh
  computeMeshNormals( scene.mesh );

  //!! iterate through the first 10 vertices and their normals
  for(unsigned int i = 0; i < 100; i++) {
    std::cout << "V: " << scene.mesh.vertices[i] <<std::endl;
    std::cout << "N: " << scene.mesh.normals[i] << std::endl;
  }
  // allocate color vectors - RGBA
  scene.mesh.vertex_color = new Vector4 [scene.mesh.num_vertices];
  // iterate through and add the color - Red
  for(unsigned int v_ctr = 0; v_ctr < scene.mesh.num_vertices; v_ctr++)
    {
      scene.mesh.vertex_color[ v_ctr ].x = 1.0;//Vector4(1.0, 0.0, 0.0, 1.0);
      scene.mesh.vertex_color[ v_ctr ].y = 0.0;
      scene.mesh.vertex_color[ v_ctr ].z = 0.0;
      scene.mesh.vertex_color[ v_ctr ].w = 1.0;
    }

}


/**
 * Compute watersurface mesh and initialize vertex normals.
 */
void OpenglProject::createWaterSurfaceMesh()
{

  // compute and create the watersurface mesh
  // iterate in the x-direction and query for height while meshifying
  // resolution is 100
  double range[2] = {-1.0, 1.0};
  double x = range[0], z = range[0];
  // resolution inclusive of both start and end points
  unsigned int resolution = 6;

  double increment = (double)( (range[1] - range[0]) / (resolution-1) );

  // allocate vertices of mesh
  this->scene.heightmap->mesh.num_vertices = resolution*resolution;
  this->scene.heightmap->mesh.vertices = new Vector3 [this->scene.heightmap->mesh.num_vertices];

  // triagular tesslation of grids: 2*(n-1)^2
  this->scene.heightmap->mesh.num_triangles = 2*( std::pow( resolution-1, 2 ) );
  this->scene.heightmap->mesh.triangles = new Triangle [this->scene.heightmap->mesh.num_triangles];

  unsigned int vertex_ctr = 0;
  unsigned int triangle_ctr =0;
  while( z < range[1] ) 
    {
      x = range[0];
      while( x < range[1] )
	{
	  // create new vertex
	  this->scene.heightmap->mesh.vertices[vertex_ctr].x = x;
	  this->scene.heightmap->mesh.vertices[vertex_ctr].z = z;
	  // query vertex for height
	  this->scene.heightmap->mesh.vertices[vertex_ctr].y = 
	    this->scene.heightmap->compute_height( Vector2(x,z) );

	  // get the indices of the grid where current vertex is the top left
	  // x --- v1
	  // |    / |
	  // |   /  |
	  // |  /   |
	  // v3 --- v2
	  unsigned int v1 = vertex_ctr + 1;
	  unsigned int v2 = vertex_ctr + resolution + 1;
	  unsigned int v3 = vertex_ctr + resolution;

	  // add triangle with neighboring vertices
	  // these new vertices will be queried when they become primary or in the end
	  
	  // triangle 1
	  this->scene.heightmap->mesh.triangles[triangle_ctr].vertices[0] = vertex_ctr;
	  this->scene.heightmap->mesh.triangles[triangle_ctr].vertices[1] = v1;
	  this->scene.heightmap->mesh.triangles[triangle_ctr].vertices[2] = v3;
	  // triangle 2
	  this->scene.heightmap->mesh.triangles[triangle_ctr + 1].vertices[0] = v1;
	  this->scene.heightmap->mesh.triangles[triangle_ctr + 1].vertices[1] = v2;
	  this->scene.heightmap->mesh.triangles[triangle_ctr + 1].vertices[2] = v3;
	  
	  // increment triangle counter
	  triangle_ctr += 2;
	  // increment x
	  x = x + increment;
	  // increment vertex counter
	  vertex_ctr++;

	  // check whether current iteration is the last grid
	  if ( x == range[1] ) {
	    // query the next vertex (v1) in that case
	    this->scene.heightmap->mesh.vertices[vertex_ctr].x = range[1];
	    this->scene.heightmap->mesh.vertices[vertex_ctr].z = z;
	    // query vertex for height
	    this->scene.heightmap->mesh.vertices[vertex_ctr].y = 
	      this->scene.heightmap->compute_height( Vector2(x,z) );

	    // increment vertex counter
	    vertex_ctr++;

	  }


	}

      // increment z
      z = z + increment; 

      // check whether last z-row has been reached
      if ( z == range[1] ) {

	x = range[0];
	// iterate through the x-columns and query all vertices
	while ( x <= range[1] )
	  {
	    this->scene.heightmap->mesh.vertices[vertex_ctr].x = x;
	    this->scene.heightmap->mesh.vertices[vertex_ctr].z = range[1];
	    // query vertex for height
	    this->scene.heightmap->mesh.vertices[vertex_ctr].y = 
	      this->scene.heightmap->compute_height( Vector2(x,z) );

	    // no triangles need to be added: already added for these vertices

	    // increment x
	    x = x + increment;
	    // increment vertex counter
	    vertex_ctr++;
	  }
      }

    }

  // allocate the mesh normals
  scene.heightmap->mesh.normals = new Vector3 [scene.heightmap->mesh.num_vertices];
  // compute normals for watersurface mesh
  computeMeshNormals( scene.heightmap->mesh );

  // allocate color vectors - RGBA  
  scene.heightmap->mesh.vertex_color = new Vector4 [scene.heightmap->mesh.num_vertices];
  // add the color Blue to the mesh
  for(unsigned int v_ctr = 0; v_ctr < scene.heightmap->mesh.num_vertices; v_ctr++)
    {
      scene.heightmap->mesh.vertex_color[ v_ctr ].x = 0.0;//Vector4(1.0, 0.0, 0.0, 1.0);
      scene.heightmap->mesh.vertex_color[ v_ctr ].y = 0.0;
      scene.heightmap->mesh.vertex_color[ v_ctr ].z = 255.0;
      scene.heightmap->mesh.vertex_color[ v_ctr ].w = 1.0;
    }
}


/**
 * Clean up the project. Free any memory, etc.
 */
void OpenglProject::destroy()
{
    // TODO any cleanup code, e.g., freeing memory
}

/**
 * Perform an update step. This happens on a regular interval.
 * @param dt The time difference from the previous frame to the current.
 */
void OpenglProject::update( real_t dt )
{
    // update our heightmap
    scene.heightmap->update( dt );

    // TODO any update code, e.g. commputing heightmap mesh positions and normals
    // iterate through heightmap mesh and update height
    for (unsigned int v_ctr = 0; v_ctr < scene.heightmap->mesh.num_vertices; v_ctr++) 
      {
	Vector2 pos( scene.heightmap->mesh.vertices[v_ctr].x, 
		     scene.heightmap->mesh.vertices[v_ctr].z );

	scene.heightmap->mesh.vertices[v_ctr].y = scene.heightmap->compute_height(pos);
      }

    // recompute watersurface mesh normals
    computeMeshNormals( scene.heightmap->mesh );

}

/**
 * Clear the screen, then render the mesh using the given camera.
 * @param camera The logical camera to use.
 * @see math/camera.hpp
 */
void OpenglProject::render( const Camera* camera )
{
  // TODO render code
  // clear the buffer
  glClear( GL_COLOR_BUFFER_BIT |
  	   GL_DEPTH_BUFFER_BIT );

  // set the camera
  setCamera( camera );


  //renderPool();
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  glColor4f(255.0/255.0, 1.0/255.0, 1.0/255.0, 1.0);
  float mat_specular[] = {0.992157, 0.941176, 0.807843, 1.0};
  float shininess = 100;
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialf(GL_FRONT, GL_SHININESS, shininess);
  renderMesh( scene.mesh, scene.mesh_position );
  glDisable(GL_COLOR_MATERIAL);

  //renderWaterSurface();
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  glColor4f(1.0/255.0, 1.0/255.0, 255.0/255.0, 1.0);
  float mat_specular_2[] = {0.992157, 0.941176, 0.807843, 1.0};
  float shininess_2 = 100;
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular_2);
  glMaterialf(GL_FRONT, GL_SHININESS, shininess_2);
  renderMesh( scene.heightmap->mesh, scene.heightmap_position );
  glDisable(GL_COLOR_MATERIAL);


  if (glGetError() != GL_NO_ERROR)
    std::cout << "We got problems!" <<std::endl;

}


/**
 * Function to render a mesh at a givven position
 */
void OpenglProject::renderMesh( const MeshData& mesh, const PositionData& position )
{
  // apply transforms   
  glEnable( GL_NORMALIZE );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  glTranslated( position.position.x, 
		position.position.y, 
		position.position.z );
  glRotated( position.orientation.w,
	     position.orientation.x,
	     position.orientation.y,
	     position.orientation.z );
  glScaled( position.scale.x,
  	    position.scale.y,
  	    position.scale.z );

  // render the mesh using Vertex Arrays
  // enable functionality
  glEnableClientState( GL_VERTEX_ARRAY );  
  glEnableClientState( GL_COLOR_ARRAY );
  glEnableClientState( GL_NORMAL_ARRAY );

  // pass the vertex pointer
  glVertexPointer( 3,
		   GL_DOUBLE, // real_t is a typedef of double
		   0,
		   mesh.vertices );

  // pass the color pointer
  glColorPointer( 4,
  		 GL_DOUBLE, 
  		 0, 
		 mesh.vertex_color );

  // add the normals to the vertex arrays
  glNormalPointer( GL_DOUBLE,
		   0,
		   mesh.normals );

  // create triangle meshes using the indices specified in triangles
  for (unsigned int i = 0; i < mesh.num_triangles; i++) {
    glDrawElements( GL_TRIANGLES,
		    3, 
		    GL_UNSIGNED_INT, 
		    &mesh.triangles[i].vertices[0] );
  }

  // disable the client states
  glDisableClientState( GL_VERTEX_ARRAY );  
  glDisableClientState( GL_NORMAL_ARRAY );  
  glDisableClientState( GL_COLOR_ARRAY );
}


/**
 * Function to compute and update the normals vectors for vertices in a mesh
 * CAUTION: the normals array must already be allocated to the size num_vertices
 * function will corrupt memory otherwise
 */
void OpenglProject::computeMeshNormals( MeshData& mesh )
{

  //use VertexNormal struct to efficiently average normals for eachh vertex
  VertexNormalStruct* normal_accumulator = new VertexNormalStruct [ mesh.num_vertices ];

  // iterate through triangles and compute the normal vectors
  for(unsigned int i = 0; i < mesh.num_triangles; i ++)  
    {
      // normal vector = (b-a)x(c-a)
      Vector3 a = mesh.vertices[ mesh.triangles[i].vertices[0] ];
      Vector3 b = mesh.vertices[ mesh.triangles[i].vertices[1] ];
      Vector3 c = mesh.vertices[ mesh.triangles[i].vertices[2] ];

      // compute normal
      Vector3 current_normal = normalize( cross( b-a, c-a ) );

      // add the normal vector to the respective vertices
      for(unsigned int triangle_iter = 0; triangle_iter < 3; triangle_iter++) 
	{
	  int vertex_index = mesh.triangles[i].vertices[triangle_iter];
	  normal_accumulator[ vertex_index ].normal += current_normal;
	  normal_accumulator[ vertex_index ].num_adjacent_triangles += 1;
	}
      
    }

  //!! Function assumes that the normals array has already been allocated
  // iterate through all the vertices and average the normals
  for(unsigned int v_ctr = 0; v_ctr < mesh.num_vertices; v_ctr++)
    {
      
      mesh.normals[ v_ctr ] = 
	normalize ( normal_accumulator[ v_ctr ].normal / normal_accumulator[ v_ctr ].num_adjacent_triangles );

    }

  // delete locally allocated memory
  delete normal_accumulator;
}


/**
 * Function to set the camera by applying viewing transforms
 */
void OpenglProject::setCamera( const Camera* camera )
{

  // apply camera transforms
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( camera->get_fov_degrees(),
		  camera->get_aspect_ratio(),
		  camera->get_near_clip(),
		  camera->get_far_clip() );

  gluLookAt( camera->get_position().x,
	     camera->get_position().y,
	     camera->get_position().z,
	     camera->get_direction().x,
	     camera->get_direction().y,
	     camera->get_direction().z,
	     camera->get_up().x,
	     camera->get_up().y,
	     camera->get_up().z );
}


} /* _462 */

// NOTE: The normals calculations might have some bugs. The lighting does not look quite right
