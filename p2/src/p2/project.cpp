/**
 * @file project.cpp
 * @brief Geometry project
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "project.hpp"
#include "application/opengl.hpp"

/*
  A namespace declaration. All project files use this namespace.
  Add this declaration (and its closing) to all source/headers you create.
  Note that all #includes should be BEFORE the namespace declaration.
*/
namespace _462 {

  // definitions of functions for the GeometryProject class

  // constructor, invoked when object is allocated
  GeometryProject::GeometryProject() { }

  // destructor, invoked when object is de-allocated
  GeometryProject::~GeometryProject() { }

  /**
   * Initialize the project, doing any necessary opengl initialization.
   * @param camera An already-initialized camera.
   * @param mesh The mesh to be rendered and subdivided.
   * @param texture_filename The filename of the texture to use with the mesh.
   *  Is null if there is no texture data with the mesh or no texture filename
   *  was passed in the arguments, in which case textures should not be used.
   * @return true on success, false on error.
   */
  bool GeometryProject::initialize( const Camera* camera, const MeshData* mesh, const char* texture_filename )
  {
    this->mesh = *mesh;

    // TODO opengl initialization code

    // position camera
    setCamera( camera );

    // add lighting to the scene
    GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };

    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);

    // a vertex index is mapped to the corresponding triangle indices
    vertexToTriangleMapping( );

    return true;
  }

  /**
   * Clean up the project. Free any memory, etc.
   */
  void GeometryProject::destroy()
  {
    // TODO any cleanup code
  }

  /**
   * Clear the screen, then render the mesh using the given camera.
   * @param camera The logical camera to use.
   * @see scene/camera.hpp
   */
  void GeometryProject::render( const Camera* camera )
  {
    // TODO render code

    // set the camera
    setCamera( camera );

    // render mesh
    renderMesh();
  }

  /**
   * Subdivide the mesh that we are rendering using Loop subdivison.
   */
  void GeometryProject::subdivide()
  {
    // TODO perform a single subdivision.

    subdivisionData subdivision;

    // allocate the split triangle data structure
    subdivision.splitTriangleContainer = new TriangleSplit[ mesh.num_triangles ];
    // allocate the array for new odd vertices formed by splitting edges
    // #edges = 2*#triangles + 1;
    subdivision.num_oddVertices = (2*mesh.num_triangles) + 1;
    subdivision.oddVertices = new Vertex[ subdivision.num_oddVertices ];
    // allocate array to store new triangles
    subdivision.num_newTriangles = 4*mesh.num_triangles;
    subdivision.newTriangles = new Triangle[ subdivision.num_newTriangles ];


    std::cout << "# Orig Triangles: " << mesh.num_triangles <<std::endl;
    std::cout << "# Orig Vertices: " << mesh.num_vertices  <<std::endl;


    initiateLoopSubdivision( subdivision );


    std::cout << "# Triangles: " << mesh.num_triangles <<std::endl;
    std::cout << "# Vertices: " << mesh.num_vertices  <<std::endl;

  }


  /**
   *
   */  
  void GeometryProject::initiateLoopSubdivision( subdivisionData& subdivision )
  {

    // instantiate a counter to keep track off odd vertices
    unsigned int oddVertex_ctr = 0;

    // iterate through the triangles and access pairs of vertices
    for (unsigned int triangle_ind = 0; triangle_ind < mesh.num_triangles; triangle_ind++)
      {
	// pull out pointers to the corresponding triangle structures
	Triangle* rawTriangle = &mesh.triangles[ triangle_ind ];
	TriangleSplit* splitTriangle = &subdivision.splitTriangleContainer[ triangle_ind ];
 
	// iterate through the edge pairs in triangle
	for (unsigned int edge_ind = 0; edge_ind < 3; edge_ind++)
	  {
	    // initialize a variable to hold the new odd vertex that will be created
	    Vertex newOddVertex;

	    // pull first vertex triangles
	    unsigned int ev1 = rawTriangle->vertices[ edge_ind%3 ];

	    // pull second vertex triangles
	    unsigned int ev2 = rawTriangle->vertices[ (edge_ind+1)%3 ];

	    // third vertex not in edge
	    unsigned int v3 = rawTriangle->vertices[ (edge_ind+2)%3 ];

	    // index of the triangle that neighbors current edge 
	    int neighborTriangleInd = findNeighborTriangle( ev1, ev2, triangle_ind );

	    // collect all vertices required to split edge
	    Vertex parentVertices[2];
	    parentVertices[ 0 ] = mesh.vertices[ ev1 ];
	    parentVertices[ 1 ] = mesh.vertices[ ev2 ];

	    // add the new triangle index as neighbor
	    if ( neighborTriangleInd >= 0 ) {

	      // add the triangle index in the neighboring triangles list
	      splitTriangle->neighbor_triangles[ edge_ind ] = neighborTriangleInd;

	      // get pointers to the neighboring triangle structures
	      Triangle* rawTriangleNeighbor = &mesh.triangles[ neighborTriangleInd ];
	      TriangleSplit* splitTriangleNeighbor = &subdivision.splitTriangleContainer[ neighborTriangleInd ];

	      // create a array of support vectors
	      Vertex supportVertices[2];
	      supportVertices[0] = mesh.vertices[ v3 ];

	      // get the non-common vertex from neighboring triangle
	      unsigned int v_ctr = 0;
	      while( v_ctr < 3)
		{
		  unsigned int neighboring_vertex = rawTriangleNeighbor->vertices[ v_ctr ];
		  if ( ( neighboring_vertex != ev1 ) &&
		      ( neighboring_vertex != ev2 ) ) {
		    supportVertices[ 1 ] = mesh.vertices[ neighboring_vertex ];
		    break;
		  }
		  v_ctr ++;
		}

	      // check whether the neighbor has already split this edge
	      if ( splitTriangleNeighbor->isSplit ) {
		// borrow the split vertex
		// if edge order is :0-1, 1-2, 2-0
		// we must find the extra vertex and look at next edge index
		splitTriangle->splitEdgeVertices[ edge_ind ] = splitTriangleNeighbor->splitEdgeVertices[ (v_ctr+1)%3 ];

	      }
	      else {
		// call helper function to split edge and add the newOddVertex to container
		subdivision.oddVertices[ oddVertex_ctr ] = 
		  loopSubdivideOddVertex( parentVertices, supportVertices );
	
		splitTriangle->splitEdgeVertices[ edge_ind ] = oddVertex_ctr;
		oddVertex_ctr++;
	      }
	    }
	    else {
	      // set the on boundary flag for the 2 vertices - used in even vertices
	      mesh.vertices[ ev1 ].onBoundary = true;
	      mesh.vertices[ ev2 ].onBoundary = true;

	      splitTriangle->neighbor_triangles[ edge_ind ] = -1;
	      
	      // call helper function to split edge and add the newOddVertex to container
	      subdivision.oddVertices[ oddVertex_ctr ] = 
		loopSubdivideOddVertex( parentVertices );
	      
	      splitTriangle->splitEdgeVertices[ edge_ind ] = oddVertex_ctr;
	      oddVertex_ctr++;
	    }
	    
	    // WAIT but I am updating triangles now
	    // add the current triangle index in the vertex register
	    // this eliminates the need for future vertex to triangle mapping 
	    // in further subdivision operations
	    //subdivision.oddVertices[ splitTriangle->splitEdgeVertices[ edge_ind ] ].triangles.push_back( triangle_ind );

	    // sanity check on odd vertex counter
	    if ( oddVertex_ctr >=subdivision.num_oddVertices ) {
	      std::cout << "I have not allocated enough odd vertices !!" << std::endl;
	      std::cout << "# odd vertices = " << subdivision.num_oddVertices
			<< " and vertices added = " << oddVertex_ctr
			<< std::endl;

	      assert(0);
	    }
	    // TODO: implement functionality to add new triangles

	  }

	// set the triangle split flag to true
	splitTriangle->isSplit = true;

	// once all edges have been split create 4 new triangles
	divideTriangleByFour( subdivision, triangle_ind );

      }

    // all odd vertices have been computed at this point
    // compute (approximate) the even vertices with the old neighbors
    // no parameters need to be passed as all the neighbors are stored within the vertex
    loopSubdivideEvenVertices( );

    std::cout << "# odd vertices allocated = " << subdivision.num_oddVertices
	      << " and vertices added = " << oddVertex_ctr
	      <<" original # vertices = " << mesh.num_vertices
	      << std::endl;
    

    // modify the mesh with the subdivided triangles and vertices
    delete[] mesh.triangles;
    mesh.num_triangles = subdivision.num_newTriangles;
    // simply point the mesh triangles to the subdivision triangles
    // subdivision triangles will never be deleted
    // mesh must make sure to deallocate this memory
    // saves time for copying over the array
    mesh.triangles = &subdivision.newTriangles[0];

    // merge the subdivision vertices with the original vertices
    Vertex* subdividedVertices = new Vertex[ mesh.num_vertices + oddVertex_ctr ];
    copyElements( mesh.vertices, 
		  mesh.vertices + mesh.num_vertices, 
		  subdividedVertices );
    copyElements( subdivision.oddVertices, 
		  subdivision.oddVertices + oddVertex_ctr, 
		  subdividedVertices + mesh.num_vertices );

    // free the previously allocated container
    delete[] mesh.vertices;
    delete[] subdivision.oddVertices;
    // set the mesh vertices to point to new container
    mesh.num_vertices += oddVertex_ctr;
    mesh.vertices = &subdividedVertices[0];

    // recompute the vertex normals
    computeVertexNormals( );
  }


  /**
   *
   */  
  void GeometryProject::divideTriangleByFour( subdivisionData& subdivision,
					      unsigned int triangle_ind )
  {
    /* Figure to show triangle division
                 0
                / \
               /t-0\
             v0---- v2
             / \t-3/ \
            /t-1\ /t-2\
           1 ---v1---- 2

    */

    // pull out pointers to the corresponding triangle structures
    Triangle* rawTriangle = &mesh.triangles[ triangle_ind ];
    TriangleSplit* splitTriangle = &subdivision.splitTriangleContainer[ triangle_ind ];

    // instantiate container to store triangle
    Triangle newVertexTriangle;
    // instantiate index variable to store new triangles
    unsigned int new_triangle_index;
    // iterate through the vertices in triangle
    // counter starts from 3 in order to go back one step (3-2)%3 and access last elem
    for (unsigned int vertex_ind = 3; vertex_ind < 6; vertex_ind++)
      {
	// the two triangle edges for a vertex i are stored in i and (i-1)%3
	// of splitEdgeVertices -- look at header
	newVertexTriangle.vertices[ 0 ] = rawTriangle->vertices[ vertex_ind%3 ];
	// the new oddVertices will be appended to existing vertices container in mesh
	// therefore offset the index by size of current vertices array
	newVertexTriangle.vertices[ 1 ] = mesh.num_vertices + 
	  splitTriangle->splitEdgeVertices[ vertex_ind%3 ];
	newVertexTriangle.vertices[ 2 ] = mesh.num_vertices + 
	  splitTriangle->splitEdgeVertices[ (vertex_ind-1)%3 ];

	// add the triangle to the new triangles structure
	// offset the index by triangle index plus the vertex index
	// since 4x the # of triangles have been allocated for new triangles
	new_triangle_index = (4*triangle_ind) + (vertex_ind%3);
	subdivision.newTriangles[ new_triangle_index ] = newVertexTriangle;

      }

    // add the fourth and final triangle (t-3 in figure above)
    newVertexTriangle.vertices[ 0 ] = mesh.num_vertices + 
      splitTriangle->splitEdgeVertices[ 0 ];
    newVertexTriangle.vertices[ 1 ] = mesh.num_vertices + 
      splitTriangle->splitEdgeVertices[ 1 ];
    newVertexTriangle.vertices[ 2 ] = mesh.num_vertices + 
      splitTriangle->splitEdgeVertices[ 2 ];
    new_triangle_index = (4*triangle_ind) + 3;
    subdivision.newTriangles[ new_triangle_index ] = newVertexTriangle;

  }


  /**
   * Function uses the loop subdivision to re approximate the even vertices
   * by positioning with respect to its neighbors
   */  
  void GeometryProject::loopSubdivideEvenVertices()
  {
    // iterate through the list of currently active vertices (even)
    for (unsigned int vertex_ctr = 0; vertex_ctr < mesh.num_vertices; vertex_ctr++)
      {
	// pull out the current vertex
	Vertex currentVertex = mesh.vertices[ vertex_ctr ];
	// instantiate a container to store all neighboring vertices
	std::vector<unsigned int> neighborVertices;
	unsigned int num_neighbor_triangles = currentVertex.triangles.size();
	// iterate through the neighvoring triangles for this vertex
	for (unsigned int tri_ctr = 0; tri_ctr < num_neighbor_triangles; tri_ctr++)
	  {
	    unsigned int triangle_index = currentVertex.triangles[ tri_ctr ];
	    Triangle neighborT = mesh.triangles[ triangle_index ];
	    // iterate through the triangle
	    for(unsigned int tri_v_ctr = 0; tri_v_ctr < 3; tri_v_ctr++)
	      {
		  
		// add all the vertices (except me!)
		if(neighborT.vertices[ tri_v_ctr ] != vertex_ctr)
		  neighborVertices.push_back( neighborT.vertices[ tri_v_ctr ] );
	      }
	  }

	// instantiate a new vertex struct
	Vertex evenVertex;
	// call the appropriate re approximation routine 
	if (!currentVertex.onBoundary)
	  evenVertex = computeEvenVertex( currentVertex, neighborVertices );
	else
	  evenVertex = computeEvenBoundaryVertex( currentVertex, neighborVertices );

	// update the vertex in mesh
	mesh.vertices[ vertex_ctr ] = evenVertex;

      }

  }


  /**
   * Function to compute an even vertex given neighbors
   */
  Vertex GeometryProject::computeEvenVertex( const Vertex& currentVertex, 
					     const std::vector<unsigned int>& neighborVertices )
  {
    // compute the parameters
    unsigned int N = neighborVertices.size();
    double beta = (1.0/N)*((5.0/8.0) - pow(((3.0/8.0) +((1/4.0)*cos((2*PI)/N))), 2));

    // recompute and save the vertex
    Vertex evenVertex;
    evenVertex.position = ((1.0 - beta*N)*currentVertex.position);
    // iterate through the neighbor vertices and add their weighted contribution
    for (unsigned int n_ctr = 0; n_ctr < N; n_ctr++)
      {
	Vertex neighbor = mesh.vertices[ neighborVertices.at(n_ctr) ];
	evenVertex.position += (beta*neighbor.position);
      }

    return evenVertex;
  }

  /**
   * Function to compute even boundary vertex given neighbors 
   */
  Vertex GeometryProject::computeEvenBoundaryVertex( const Vertex& currentVertex,
						     const std::vector<unsigned int>& neighborVertices )
  {
    // there should only be two other boundary vertices within the neighbors
    Vertex supportVertices[2];

    // total number of neighbor vertices
    unsigned int N = neighborVertices.size();

    //
    unsigned int boundaryVertex_ctr = 0;
    // iterate through the neighbor vertices and fetch the other boundary vertices
    for (unsigned int n_ctr = 0; n_ctr < N; n_ctr++)
      {
	// if boundary flag is set on this vertex - 
	// we only add other boundary vertices
	if (!mesh.vertices[ neighborVertices.at(n_ctr) ].onBoundary)
	  continue;

	// sanity check on the number of boundary vertices
	if (boundaryVertex_ctr >= 2) {
	  std::cout << "There are more than 2 neighbor boundary vertices!" <<std::endl;
	  assert(0);
	}
	supportVertices[ boundaryVertex_ctr ] = 
	  mesh.vertices[ neighborVertices.at(n_ctr) ];
	
	boundaryVertex_ctr++;

      }

    // recompute and save the vertex
    Vertex evenVertex;
    evenVertex.position = ((3.0/4.0)*currentVertex.position) +
      ((1.0/8.0)*supportVertices[ 0 ].position) +
      ((1.0/8.0)*supportVertices[ 1 ].position);

    return evenVertex;

  }


  /**
   * Function uses the loop subdivision to calculate the odd vertex split
   * the two arrays must be of size 2 if allocated
   */  
  Vertex GeometryProject::loopSubdivideOddVertex( Vertex parentVertices[], 
						  Vertex supportVertices[] )
  {
    Vertex newOddVertex;
    if (supportVertices != NULL) {
      newOddVertex.position = 
	( (3.0/8.0)*parentVertices[0].position ) +
	( (3.0/8.0)*parentVertices[1].position ) +
	( (1.0/8.0)*supportVertices[0].position ) +
	( (1.0/8.0)*supportVertices[1].position );
    }
    else {
      newOddVertex.position = 
	( (1.0/2.0)*parentVertices[0].position ) +
	( (1.0/2.0)*parentVertices[1].position );
      }

    return newOddVertex;
  }


  /**
   *
   */  
  void GeometryProject::vertexToTriangleMapping()
  {

    // Inefficient way : Optimize later
    for (unsigned int vertex_ctr = 0; vertex_ctr < mesh.num_vertices; vertex_ctr++)
      mesh.vertices[ vertex_ctr ].triangles.clear();

    // iterate through the triangles
    for (unsigned int triangle_ind = 0; triangle_ind < mesh.num_triangles; triangle_ind++)
      {
	for (unsigned int ctr_vertex = 0; ctr_vertex < 3; ctr_vertex++)
	  {
	    unsigned int vertex_index = mesh.triangles[triangle_ind].vertices[ctr_vertex];
	    mesh.vertices[ vertex_index ].triangles.push_back( triangle_ind );
	  }
      }

  }


  /**
   * Function to compute the vertex normals by averaging neighboring triangle normals
   */
  void GeometryProject::computeVertexNormals() 
  {

    // remap the vertices to triangles
    // could be incorporated within main loop - optmize later
    vertexToTriangleMapping( );

    // instantiate container to store all triangle normals
    Vector3* triangleNormals = new Vector3[ mesh.num_triangles ];

    // iterate through all the triangles and compute their normals
    for(unsigned int triangle_ctr = 0; triangle_ctr < mesh.num_triangles; triangle_ctr++)
      {
	Vector3 v1 = mesh.vertices[ mesh.triangles[ triangle_ctr ].vertices[0] ].position;
	Vector3 v2 = mesh.vertices[ mesh.triangles[ triangle_ctr ].vertices[1] ].position;
	Vector3 v3 = mesh.vertices[ mesh.triangles[ triangle_ctr ].vertices[2] ].position;
	// compute the normalized normal vector for the triangle
	triangleNormals[ triangle_ctr ] = normalize( cross( v2-v1, v3-v1 ) );

      }

    // iterate through all the vertices and average the normals over their triangles
    for (unsigned int vertex_ctr = 0; vertex_ctr < mesh.num_vertices; vertex_ctr++)
      {
	// pull out the current vertex
	Vertex currentVertex = mesh.vertices[ vertex_ctr ];
	// instantiate and initialize an accumulator container to add normal vectors
	Vector3 normalAccumulator(0.0, 0.0, 0.0);
	// fetch size of neighboring triangles
	unsigned int num_neighbor_triangles = currentVertex.triangles.size();
	// iterate through the neighboring triangles for this vertex
	for (unsigned int tri_ctr = 0; tri_ctr < num_neighbor_triangles; tri_ctr++)
	  {
	    // accumulate the normal vector for these triangles
	    normalAccumulator += triangleNormals[ currentVertex.triangles[ tri_ctr ] ];
	  }
	// average to compute the normal at vertex
	currentVertex.normal = normalize( normalAccumulator/num_neighbor_triangles );
	// update the mesh with the new normal
	mesh.vertices[ vertex_ctr ] = currentVertex;
      }
  }


  /**
   * Function to search for the neighboring triangle 
   * to the edge defined by parameter vertices.
   * Returns -1 is edge is a boundary edge
   */
  int GeometryProject::findNeighborTriangle( unsigned int vertex1, unsigned int vertex2, unsigned int current_triangle )
  {
    int rv = -1;
    std::vector<unsigned int> v1Triangles = mesh.vertices[ vertex1 ].triangles;
    std::vector<unsigned int> v2Triangles = mesh.vertices[ vertex2 ].triangles;

    // iterate through vectors and find common triangles: other than self
    // thankfully the vectors are already sorted
    std::vector<unsigned int>::iterator v1It = v1Triangles.begin();
    std::vector<unsigned int>::iterator v2It = v2Triangles.begin();

    while ( (v1It  != v1Triangles.end()) ||  (v2It != v2Triangles.end()) )
      {

	if (*v1It < *v2It)
	  v1It++;
	else if (*v2It < *v1It)
	  v2It++;
	else { 
	  // make sure they are different than current triangle
	  if (*v1It == current_triangle) {
	    // increment both iterators
	    v1It++;
	    v2It++;
	  }
	  else { // we have a match
	    rv = *v1It;
	    break;
	  }
	}
      } // end while

    return rv;

  }

  /**
   *
   */  
  void GeometryProject::renderMesh()
  {

    // clear the buffer
    glClear( GL_COLOR_BUFFER_BIT |
    	     GL_DEPTH_BUFFER_BIT );
    
    /*
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glColor4f(255.0/255.0, 1.0/255.0, 1.0/255.0, 1.0);
    float mat_specular[] = {0.992157, 0.941176, 0.807843, 1.0};
    float shininess = 100;
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);
    */

    // render the mesh using Vertex Arrays
    // enable functionality 
    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_NORMAL_ARRAY );

    // pass the vertex pointer
    glVertexPointer( 3,
		     GL_DOUBLE, // real_t is a typedef of double       
		     sizeof(Vertex),
		     &mesh.vertices->position );

    // add the normals to the vertex arrays
    glNormalPointer( GL_DOUBLE,
    		     sizeof(Vertex),
    		     &mesh.vertices->normal );
    
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
    //glDisable(GL_COLOR_MATERIAL);

  }



  /**
   * Function to set the camera by applying viewing transforms   
   */
  void GeometryProject::setCamera( const Camera* camera )
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

  /**
   * Function to copy element by element of the vertex array   
   */
  Vertex* GeometryProject::copyElements( Vertex* first, Vertex* last, Vertex* result )
  {
    // iterate through the range of elements specified
    while(first != last)
      {
	// dereference and copy into result container
	*result = *first;
	// increment pointers
	result++; first++;
      }
    return result;

  }

} /* _462 */

