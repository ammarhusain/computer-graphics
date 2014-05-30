/**
 * @file project.hpp
 * @brief Geometry project
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#ifndef _462_GEOMETRY_PROJECT_HPP_
#define _462_GEOMETRY_PROJECT_HPP_

#include "math/camera.hpp"

/* STL includes */
#include <utility>
#include <vector>
#include <math.h>

#define PI 3.14159265

/*
  A namespace declaration. All project files use this namespace.
  Add this declaration (and its closing) to all source/headers you create.
  Note that all #includes should be BEFORE the namespace declaration.
*/
namespace _462 {

  struct Triangle
  {
    // index into vertex list of the 3 vertices of this triangle
    unsigned int vertices[3];

    // overloaded operator to print to screen
    friend std::ostream& operator<<( std::ostream& os, const Triangle& t )
    {
      return os << '(' << t.vertices[0] 
		<< ',' << t.vertices[1] 
		<< ',' << t.vertices[2]
		<< ')' ;
    }

  };
  

  

  struct TriangleSplit
  {
    /**
          0
         / \
        /   \
       1 --- 2
    **/

    // each edge has a neighboring triangle - unless a boundary edge
    int neighbor_triangles[3];

    // index for new vertex of each split edge
    // if triangle has vertices: 0, 1, 2
    // splitEdgeVertices: 0-1, 1-2, 2-0
    unsigned int splitEdgeVertices[3];

    // flag to indicate if the triangle has been split
    bool isSplit;

    // default constructor
    TriangleSplit() : isSplit(false) {}

    // overloaded operator to print to screen
    friend std::ostream& operator<<( std::ostream& os, const TriangleSplit& t )
    {
      return os << "Odd Vertex Indices: "
		<< '(' << t.splitEdgeVertices[0] 
		<< ',' << t.splitEdgeVertices[1] 
		<< ',' << t.splitEdgeVertices[2]
		<< ')' ;
    }
  };

  struct Vertex
  {
    // the position of the vertex
    Vector3 position;
    // the normal of the vertex
    Vector3 normal;
    // the texture coordinate of the vertex
    Vector2 texture_coord;
    // triangles that use this vertex
    std::vector <unsigned int> triangles;

    // flag to indicate if this vertex has an edge on the boundary
    bool onBoundary;
    
    // default constructor
    Vertex() : position(0.0, 0.0, 0.0), normal(0.0, 0.0, 0.0), onBoundary(false) {}

  };

  struct MeshData
  {
    // array of vertices
    Vertex* vertices;
    // size of vertex array
    size_t num_vertices;

    // array of triangles
    Triangle* triangles;
    // size of triangle array
    size_t num_triangles;
  };

  struct subdivisionData
  {
    // array of new odd vertices
    Vertex* oddVertices;
    // number of new odd vertices
    size_t num_oddVertices;

    // array of new triangle: 4x the size of original triangles
    Triangle* newTriangles;
    // size of triangles array with odd vertices
    size_t num_newTriangles;

    // container to store triangle neightbors and split edges
    TriangleSplit* splitTriangleContainer;

  };
 
  class GeometryProject
  {
  public:

    // constructor, invoked when object is created
    GeometryProject();
    // destructor, invoked when object is destroyed
    ~GeometryProject();

    // more detailed specifications for each function are in project.cpp.

    // Initialize the project, loading the mesh from the given filename.
    // Returns true on success.
    bool initialize( const Camera* camera, const MeshData* mesh, const char* texture_filename );
    // Clean up the project, free any memory, etc.
    void destroy();
    // Render the mesh using the given camera.
    void render( const Camera* camera );
    // Subdivide the mesh
    void subdivide();

    // function to render mesh
    void renderMesh();

    // position camera in scene
    void setCamera( const Camera* camera );

    /* subdivide helper function */
    // populates a data structure to map all triangles to its 3 neighbors
    void initiateLoopSubdivision( subdivisionData& subdivision );

    // creates a mapping between vertices to their triangles
    void vertexToTriangleMapping();

    // method to compute the odd vertices
    Vertex loopSubdivideOddVertex( Vertex parentVertices[], 
			     Vertex supportVertices[] = NULL);

    // method to approximate and recompute the even vertices
    void loopSubdivideEvenVertices();

    // method to compute the normal vectors of a mesh
    void computeVertexNormals();

    // method to approximate an even vertex
    Vertex computeEvenVertex( const Vertex& currentVertex,
			      const std::vector<unsigned int>& neighborVertices );

    // method to approximate an even boundary vertex
    Vertex computeEvenBoundaryVertex( const Vertex& currentVertex,
				      const std::vector<unsigned int>& neighborVertices );

    // 
    int findNeighborTriangle( unsigned int vertex1, 
			      unsigned int vertex2,
			      unsigned int current_triangle );

    // function to create 4 new triangles from an existing triangle
    void divideTriangleByFour( subdivisionData& subdivision,
						unsigned int triangle_ind );
      
    Vertex* copyElements( Vertex* first, Vertex* last, Vertex* result );

  private:

    MeshData mesh;

    // TODO add any other private members/functions here.

    // since this has no meaningful assignment/copy, prevent the compiler from
    // automatically generating those functions
    GeometryProject( const GeometryProject& );
    GeometryProject& operator=( const GeometryProject& );
  };

} /* _462 */

#endif /* _462_GEOMETRY_PROJECT_HPP_ */

