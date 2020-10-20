/************************************************************************
     File:        Maze.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for Maze class. Manages the maze.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "Maze.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <FL/Fl.h>
#include <FL/fl_draw.h>

const char Maze::X = 0;
const char Maze::Y = 1;
const char Maze::Z = 2;

const float Maze::BUFFER = 0.1f;


//**********************************************************************
//
// * Constructor for the maze exception
//======================================================================
MazeException::
MazeException(const char *m)
//======================================================================
{
	message = new char[strlen(m) + 4];
	strcpy(message, m);
}


//**********************************************************************
//
// * Constructor to create the default maze
//======================================================================
Maze::
Maze(const int nx, const int ny, const float sx, const float sy)
//======================================================================
{
	// Build the connectivity structure.
	Build_Connectivity(nx, ny, sx, sy);

	// Make edges transparent to create a maze.
	Build_Maze();

	// Set the extents of the maze
	Set_Extents();

	// Default values for the viewer.
	viewer_posn[X] = viewer_posn[Y] = viewer_posn[Z] = 0.0;
	viewer_dir = 0.0;
	viewer_fov = 45.0;

	// Always start on the 0th frame.
	frame_num = 0;
}


//**********************************************************************
//
// * Construtor to read in precreated maze
//======================================================================
Maze::
Maze(const char *filename)
//======================================================================
{
	char    err_string[128];
	FILE    *f;
	int	    i;

	// Open the file
	if ( ! ( f = fopen(filename, "r") ) )
		throw new MazeException("Maze: Couldn't open file");

	// Get the total number of vertices
	if ( fscanf(f, "%d", &num_vertices) != 1 )
		throw new MazeException("Maze: Couldn't read number of vertices");

	// Read in each vertices
	vertices = new Vertex*[num_vertices];
	for ( i = 0 ; i < num_vertices ; i++ ) {
		float x, y;
		if ( fscanf(f, "%g %g", &x, &y) != 2 )	{
			sprintf(err_string, "Maze: Couldn't read vertex number %d", i);
			throw new MazeException(err_string);
		}
		vertices[i] = new Vertex(i, x, y);
	}

	// Get the number of edges
	if ( fscanf(f, "%d", &num_edges) != 1 )
		throw new MazeException("Maze: Couldn't read number of edges");

	// read in all edges
	edges = new Edge*[num_edges];
	for ( i = 0 ; i < num_edges ; i++ ){
		int     vs, ve, cl, cr, o;
		float	r, g, b;
		if ( fscanf(f, "%d %d %d %d %d %g %g %g",
						&vs, &ve, &cl, &cr, &o, &r, &g, &b) != 8) {
			sprintf(err_string, "Maze: Couldn't read edge number %d", i);
			throw new MazeException(err_string);
		}
		edges[i] = new Edge(i, vertices[vs], vertices[ve], r, g, b);
		edges[i]->Add_Cell((Cell*)cl, Edge::LEFT);
		edges[i]->Add_Cell((Cell*)cr, Edge::RIGHT);
		edges[i]->opaque = o ? true : false;
	}

	// Read in the number of cells
	if ( fscanf(f, "%d", &num_cells) != 1 )
		throw new MazeException("Maze: Couldn't read number of cells");


	// Read in all cells
	cells = new Cell*[num_cells];
	for ( i = 0 ; i < num_cells ; i++ )	{
		int epx, epy, emx, emy;
		if ( fscanf(f, "%d %d %d %d", &epx, &epy, &emx, &emy) != 4 ){
			sprintf(err_string, "Maze: Couldn't read cell number %d", i);
			throw new MazeException(err_string);
		}
		cells[i] = new Cell(i, epx >= 0 ? edges[epx] : NULL,
									epy >= 0 ? edges[epy] : NULL,
									emx >= 0 ? edges[emx] : NULL,
									emy >= 0 ? edges[emy] : NULL);
		if ( cells[i]->edges[0] ) {
			if ( cells[i]->edges[0]->neighbors[0] == (Cell*)i )
				cells[i]->edges[0]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[0]->neighbors[1] == (Cell*)i )
				cells[i]->edges[0]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
						  "Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[0]->index);
				throw new MazeException(err_string);
			}
		}

		if ( cells[i]->edges[1] )	{
			if ( cells[i]->edges[1]->neighbors[0] == (Cell*)i )
				cells[i]->edges[1]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[1]->neighbors[1] == (Cell*)i )
				cells[i]->edges[1]->neighbors[1] = cells[i];
			else {
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[1]->index);
				throw new MazeException(err_string);
			}
		}
		if ( cells[i]->edges[2] ) {
			if ( cells[i]->edges[2]->neighbors[0] == (Cell*)i )
				cells[i]->edges[2]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[2]->neighbors[1] == (Cell*)i )
				cells[i]->edges[2]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[2]->index);
				throw new MazeException(err_string);
			}
		}
		if ( cells[i]->edges[3] ) {
			if ( cells[i]->edges[3]->neighbors[0] == (Cell*)i )
				cells[i]->edges[3]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[3]->neighbors[1] == (Cell*)i )
				cells[i]->edges[3]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[3]->index);
				throw new MazeException(err_string);
			}
		}
	}

	if ( fscanf(f, "%g %g %g %g %g",
					 &(viewer_posn[X]), &(viewer_posn[Y]), &(viewer_posn[Z]),
					 &(viewer_dir), &(viewer_fov)) != 5 )
		throw new MazeException("Maze: Error reading view information.");

	// Some edges have no neighbor on one side, so be sure to set their
	// pointers to NULL. (They were set at -1 by the save/load process.)
	for ( i = 0 ; i < num_edges ; i++ )	{
		if ( edges[i]->neighbors[0] == (Cell*)-1 )
			edges[i]->neighbors[0] = NULL;
		if ( edges[i]->neighbors[1] == (Cell*)-1 )
			edges[i]->neighbors[1] = NULL;
	}

	fclose(f);

	Set_Extents();

	// Figure out which cell the viewer is in, starting off by guessing the
	// 0th cell.
	Find_View_Cell(cells[0]);

	frame_num = 0;
}


//**********************************************************************
//
// * Destructor must free all the memory allocated.
//======================================================================
Maze::
~Maze(void)
//======================================================================
{
	int i;

	for ( i = 0 ; i < num_vertices ; i++ )
		delete vertices[i];
	delete[] vertices;

	for ( i = 0 ; i < num_edges ; i++ )
		delete edges[i];
	delete[] edges;

	for ( i = 0 ; i < num_cells ; i++ )
		delete cells[i];
	delete[] cells;
}


//**********************************************************************
//
// * Randomly generate the edge's opaque and transparency for an empty maze
//======================================================================
void Maze::
Build_Connectivity(const int num_x, const int num_y,
                   const float sx, const float sy)
//======================================================================
{
	int	i, j, k;
	int edge_i;

	// Ugly code to allocate all the memory for a new maze and to associate
	// edges with vertices and faces with edges.

	// Allocate and position the vertices.
	num_vertices = ( num_x + 1 ) * ( num_y + 1 );
	vertices = new Vertex*[num_vertices];
	k = 0;
	for ( i = 0 ; i < num_y + 1 ; i++ ) {
		for ( j = 0 ; j < num_x + 1 ; j++ )	{
			vertices[k] = new Vertex(k, j * sx, i * sy);
			k++;
		}
	}

	// Allocate the edges, and associate them with their vertices.
	// Edges in the x direction get the first num_x * ( num_y + 1 ) indices,
	// edges in the y direction get the rest.
	num_edges = (num_x+1)*num_y + (num_y+1)*num_x;
	edges = new Edge*[num_edges];
	k = 0;
	for ( i = 0 ; i < num_y + 1 ; i++ ) {
		int row = i * ( num_x + 1 );
		for ( j = 0 ; j < num_x ; j++ ) {
			int vs = row + j;
			int ve = row + j + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	edge_i = k;
	for ( i = 0 ; i < num_y ; i++ ) {
		int row = i * ( num_x + 1 );
		for ( j = 0 ; j < num_x + 1 ; j++ )	{
			int vs = row + j;
			int ve = row + j + num_x + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	// Allocate the cells and associate them with their edges.
	num_cells = num_x * num_y;
	cells = new Cell*[num_cells];
	k = 0;
	for ( i = 0 ; i < num_y ; i++ ) {
		int row_x = i * ( num_x + 1 );
		int row_y = i * num_x;
		for ( j = 0 ; j < num_x ; j++ )	{
			int px = edge_i + row_x + 1 + j;
			int py = row_y + j + num_x;
			int mx = edge_i + row_x + j;
			int my = row_y + j;
			cells[k] = new Cell(k, edges[px], edges[py], edges[mx], edges[my]);
			edges[px]->Add_Cell(cells[k], Edge::LEFT);
			edges[py]->Add_Cell(cells[k], Edge::RIGHT);
			edges[mx]->Add_Cell(cells[k], Edge::RIGHT);
			edges[my]->Add_Cell(cells[k], Edge::LEFT);
			k++;
		}
	}
}


//**********************************************************************
//
// * Add edges from cell to the set that are available for removal to
//   grow the maze.
//======================================================================
static void
Add_To_Available(Cell *cell, int *available, int &num_available)
//======================================================================
{
	int i, j;

	// Add edges from cell to the set that are available for removal to
	// grow the maze.

	for ( i = 0 ; i < 4 ; i++ ){
		Cell    *neighbor = cell->edges[i]->Neighbor(cell);

		if ( neighbor && ! neighbor->counter )	{
			int candidate = cell->edges[i]->index;
			for ( j = 0 ; j < num_available ; j++ )
				if ( candidate == available[j] ) {
					printf("Breaking early\n");
					break;
			}
			if ( j == num_available )  {
				available[num_available] = candidate;
				num_available++;
			}
		}
	}

	cell->counter = 1;
}


//**********************************************************************
//
// * Grow a maze by removing candidate edges until all the cells are
//   connected. The edges are not actually removed, they are just made
//   transparent.
//======================================================================
void Maze::
Build_Maze()
//======================================================================
{
	Cell    *to_expand;
	int     index;
	int     *available = new int[num_edges];
	int     num_available = 0;
	int	    num_visited;
	int	    i;

	srand(time(NULL));

	// Choose a random starting cell.
	index = (int)floor((rand() / (float)RAND_MAX) * num_cells);
	to_expand = cells[index];
	Add_To_Available(to_expand, available, num_available);
	num_visited = 1;

	// Join cells up by making edges opaque.
	while ( num_visited < num_cells && num_available > 0 ) {
		int ei;

		index = (int)floor((rand() / (float)RAND_MAX) * num_available);
		to_expand = NULL;

		ei = available[index];

		if ( edges[ei]->neighbors[0] && 
			 !edges[ei]->neighbors[0]->counter )
			to_expand = edges[ei]->neighbors[0];
		else if ( edges[ei]->neighbors[1] && 
			 !edges[ei]->neighbors[1]->counter )
			to_expand = edges[ei]->neighbors[1];

		if ( to_expand ) {
			edges[ei]->opaque = false;
			Add_To_Available(to_expand, available, num_available);
			num_visited++;
		}

		available[index] = available[num_available-1];
		num_available--;
	}

	for ( i = 0 ; i < num_cells ; i++ )
		cells[i]->counter = 0;
}


//**********************************************************************
//
// * Go through all the vertices looking for the minimum and maximum
//   extents of the maze.
//======================================================================
void Maze::
Set_Extents(void)
//======================================================================
{
	int i;

	min_xp = vertices[0]->posn[Vertex::X];
	max_xp = vertices[0]->posn[Vertex::X];
	min_yp = vertices[0]->posn[Vertex::Y];
	max_yp = vertices[0]->posn[Vertex::Y];
	for ( i = 1 ; i < num_vertices ; i++ ) {
		if ( vertices[i]->posn[Vertex::X] > max_xp )
			 max_xp = vertices[i]->posn[Vertex::X];
		if ( vertices[i]->posn[Vertex::X] < min_xp )
			 min_xp = vertices[i]->posn[Vertex::X];
		if ( vertices[i]->posn[Vertex::Y] > max_yp )
			 max_yp = vertices[i]->posn[Vertex::Y];
		if ( vertices[i]->posn[Vertex::Y] < min_yp )
			 min_yp = vertices[i]->posn[Vertex::Y];
    }
}


//**********************************************************************
//
// * Figure out which cell the view is in, using seed_cell as an
//   initial guess. This procedure works by repeatedly checking
//   whether the viewpoint is in the current cell. If it is, we're
//   done. If not, Point_In_Cell returns in new_cell the next cell
//   to test. The new cell is the one on the other side of an edge
//   that the point is "outside" (meaning that it might be inside the
//   new cell).
//======================================================================
void Maze::
Find_View_Cell(Cell *seed_cell)
//======================================================================
{
	Cell    *new_cell;

	// 
	while ( ! ( seed_cell->Point_In_Cell(viewer_posn[X], viewer_posn[Y],
													 viewer_posn[Z], new_cell) ) ) {
		if ( new_cell == 0 ) {
			// The viewer is outside the top or bottom of the maze.
			throw new MazeException("Maze: View not in maze\n");
		}

		seed_cell = new_cell;
    }
    
    view_cell = seed_cell;
}


//**********************************************************************
//
// * Move the viewer's position. This method will do collision detection
//   between the viewer's location and the walls of the maze and prevent
//   the viewer from passing through walls.
//======================================================================
void Maze::
Move_View_Posn(const float dx, const float dy, const float dz)
//======================================================================
{
	Cell    *new_cell;
	float   xs, ys, zs, xe, ye, ze;

	// Move the viewer by the given amount. This does collision testing to
	// prevent walking through walls. It also keeps track of which cells the
	// viewer is in.

	// Set up a line segment from the start to end points of the motion.
	xs = viewer_posn[X];
	ys = viewer_posn[Y];
	zs = viewer_posn[Z];
	xe = xs + dx;
	ye = ys + dy;
	ze = zs + dz;

	// Fix the z to keep it in the maze.
	if ( ze > 1.0f - BUFFER )
		ze = 1.0f - BUFFER;
	if ( ze < BUFFER - 1.0f )
		ze = BUFFER - 1.0f;

	// Clip_To_Cell clips the motion segment to the view_cell if the
	// segment intersects an opaque edge. If the segment intersects
	// a transparent edge (through which it can pass), then it clips
	// the motion segment so that it _starts_ at the transparent edge,
	// and it returns the cell the viewer is entering. We keep going
	// until Clip_To_Cell returns NULL, meaning we've done as much of
	// the motion as is possible without passing through walls.
	while ( ( new_cell = view_cell->Clip_To_Cell(xs, ys, xe, ye, BUFFER) ) )
		view_cell = new_cell;

	// The viewer is at the end of the motion segment, which may have
	// been clipped.
	viewer_posn[X] = xe;
	viewer_posn[Y] = ye;
	viewer_posn[Z] = ze;
}

//**********************************************************************
//
// * Set the viewer's location 
//======================================================================
void Maze::
Set_View_Posn(float x, float y, float z)
//======================================================================
{
	// First make sure it's in some cell.
	// This assumes that the maze is rectangular.
	if ( x < min_xp + BUFFER )
		x = min_xp + BUFFER;
	if ( x > max_xp - BUFFER )
		x = max_xp - BUFFER;
	if ( y < min_yp + BUFFER )
		y = min_yp + BUFFER;
	if ( y > max_yp - BUFFER )
		y = max_yp - BUFFER;
	if ( z < -1.0f + BUFFER )
		z = -1.0f + BUFFER;
	if ( z > 1.0f - BUFFER )
		z = 1.0f - BUFFER;

	viewer_posn[X] = x;
	viewer_posn[Y] = y;
	viewer_posn[Z] = z;

	// Figure out which cell we're in.
	Find_View_Cell(cells[0]);
}


//**********************************************************************
//
// * Set the angle in which the viewer is looking.
//======================================================================
void Maze::
Set_View_Dir(const float d)
//======================================================================
{
	viewer_dir = d;
}


//**********************************************************************
//
// * Set the horizontal field of view.
//======================================================================
void Maze::
Set_View_FOV(const float f)
//======================================================================
{
	viewer_fov = f;
}


//**********************************************************************
//
// * Draws the map view of the maze. It is passed the minimum and maximum
//   corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Map(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i;

	// Figure out scaling factors and the effective height of the window.
	scale_x = ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y = ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	// Draw all the opaque edges.
	for ( i = 0 ; i < num_edges ; i++ )
		if ( edges[i]->opaque )	{
			float   x1, y1, x2, y2;

			x1 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
			y1 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
			x2 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
			y2 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

			fl_color((unsigned char)floor(edges[i]->color[0] * 255.0),
					 (unsigned char)floor(edges[i]->color[1] * 255.0),
					 (unsigned char)floor(edges[i]->color[2] * 255.0));
			fl_line_style(FL_SOLID);
			fl_line(min_x + (int)floor((x1 - min_xp) * scale),
					  min_y + height - (int)floor((y1 - min_yp) * scale),
					  min_x + (int)floor((x2 - min_xp) * scale),
					  min_y + height - (int)floor((y2 - min_yp) * scale));
		}
}


//**********************************************************************
//
// * Draws the first-person view of the maze. It is passed the focal distance.
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================
void Maze::
Draw_View(const float focal_dist, float* projection_matrix, float* modelview_matrix)
//======================================================================
{
	frame_num++;

	// glClear(GL_DEPTH_BUFFER_BIT);
	// glEnable(GL_DEPTH_TEST);
	int max = max_xp > max_yp ? max_xp : max_yp;
	LineSeg leftFru(viewer_posn[X], viewer_posn[Y], viewer_posn[X] +  2 * max * (cos(To_Radians(viewer_dir + viewer_fov / 2.0))), viewer_posn[Y] + 2 * max * (sin(To_Radians(viewer_dir + viewer_fov / 2.0))));
	LineSeg rightFru(viewer_posn[X], viewer_posn[Y], viewer_posn[X] + 2 * max * (cos(To_Radians(viewer_dir - viewer_fov / 2.0))), viewer_posn[Y] + 2 * max * (sin(To_Radians(viewer_dir - viewer_fov / 2.0))));

	// cout << "\n\nstart\n";
	// cout << viewer_dir << endl;
	// cout << "("<< leftFru.start[0] << " , " << leftFru.start[1] << ") (" << leftFru.end[0] << " , " << leftFru.end[1] << ")" << endl;
	// cout << "("<< rightFru.start[0] << " , " << rightFru.start[1] << ") (" << rightFru.end[0] << " , " << rightFru.end[1] << ")" << endl;
	// cout << max << endl;
	// cout << leftFru.start[0] << " " << leftFru.start[1] << endl;
	// cout << leftFru.end[0] << " " << leftFru.end[1] << endl;
	// cout << rightFru.end[0] << " " << rightFru.end[1] << endl;

	// for (int i = 0; i < (int)this->num_edges; i++) {
	// 	float edge_start[2] = {
	// 		this->edges[i]->endpoints[Edge::START]->posn[Vertex::X],
	// 		this->edges[i]->endpoints[Edge::START]->posn[Vertex::Y]
	// 	};
	// 	float edge_end[2] = {
	// 		this->edges[i]->endpoints[Edge::END]->posn[Vertex::X],
	// 		this->edges[i]->endpoints[Edge::END]->posn[Vertex::Y]
	// 	};
	// 	float color[3] = { this->edges[i]->color[0], this->edges[i]->color[1], this->edges[i]->color[2] };
	// 	LineSeg edge(edges[i]);

	// 	if (this->edges[i]->opaque) {
	// 		// clip then draw
	// 		float lpoint = leftFru.Cross_Param(edge), rpoint = rightFru.Cross_Param(edge);
	// 		cout << lpoint << " " << rpoint << endl;
	// 		if ((lpoint > 1 || lpoint < 0) && (rpoint > 1 || rpoint < 0)) {
	// 			continue;
	// 		}

	// 		cout << "intersect\n";
			
	// 		Draw_Wall(edge_start, edge_end, color, projection_matrix, modelview_matrix);
	// 	}
	// }

	//###################################################################
	// TODO
	// The rest is up to you!
	//###################################################################
	
	// cout << "start" << endl;
	for (int i = 0; i < num_cells; ++i) {
		cells[i]->visited = false;
	}
	Find_View_Cell(Maze::view_cell);
	Draw_Cell(Maze::view_cell, projection_matrix, modelview_matrix, leftFru, rightFru);
	// cout << "end" << endl;
}

void Maze::Draw_Cell(Cell* cell, float* projection_matrix, float* modelview_matrix, LineSeg& leftFru, LineSeg& rightFru) {
	cout << "draw cell\n";
	cell->visited = true;
	for (int i = 0; i < 4; ++i) {
		if (cell->edges[i]->opaque) {
			float edge_start[2] = {
				cell->edges[i]->endpoints[Edge::START]->posn[Vertex::X],
				cell->edges[i]->endpoints[Edge::START]->posn[Vertex::Y]
			};
			float edge_end[2] = {
				cell->edges[i]->endpoints[Edge::END]->posn[Vertex::X],
				cell->edges[i]->endpoints[Edge::END]->posn[Vertex::Y]
			};

			LineSeg edge(cell->edges[i]);
			float leftIntersectPoint = leftFru.Cross_Param(edge), rightIntersectPoint = rightFru.Cross_Param(edge);
			float leftFruEdgeIntersect[2] = { 0.0 }, rightFruEdgeIntersect[2] = { 0.0 };
	
			bool leftFruInEdge = LineSeg::get_line_intersection(leftFru.start[0], leftFru.start[1], leftFru.end[0], leftFru.end[1],
				edge_start[0], edge_start[1], edge_end[0], edge_end[1], &leftFruEdgeIntersect[0], &leftFruEdgeIntersect[1]);
			bool rightFruInEdge = LineSeg::get_line_intersection(rightFru.start[0], rightFru.start[1], rightFru.end[0], rightFru.end[1],
				edge_start[0], edge_start[1], edge_end[0], edge_end[1], &rightFruEdgeIntersect[0], &rightFruEdgeIntersect[1]);

			// if start is inside frustum, set end with intersection
			// else if end is inside frustum, set start with intersection
			// else if both arenot inside frustim, set them with intersections
			int startAngle = To_Degrees(atan2(edge_start[1] - viewer_posn[Y], edge_start[0] - viewer_posn[X]));
			int endAngle = To_Degrees(atan2(edge_end[1] - viewer_posn[Y], edge_end[0] - viewer_posn[X]));
			int leftAngle = (viewer_dir + viewer_fov / 2.0), rightAngle = (viewer_dir - viewer_fov / 2.0);

			if (startAngle < 0.0) startAngle += 360;
			if (endAngle < 0.0) endAngle += 360;
			if (leftAngle < 0.0) leftAngle += 360;
			if (rightAngle < 0.0) rightAngle += 360;
			startAngle %= 360;
			endAngle %= 360;
			leftAngle %= 360;
			rightAngle %= 360;

			bool startInFru = (leftAngle > rightAngle) ? (rightAngle <= startAngle && startAngle <= leftAngle) : (rightAngle <= startAngle || startAngle <= leftAngle);
			bool endInFru = (leftAngle > rightAngle) ? (rightAngle <= endAngle && endAngle <= leftAngle) : (rightAngle <= endAngle || endAngle <= leftAngle);

			// ///
			// cout << "---\n";
			// cout << "start edge " << edge_start[0] << " , " << edge_start[1] << endl;			
			// cout << "end edge " << edge_end[0] << " , " << edge_end[1] << endl;			
			// cout << "IntersectPointValue " << leftIntersectPoint << " , " << rightIntersectPoint << endl;
			// cout << "FruInEdge " << leftFruInEdge << " , " << rightFruInEdge << endl;
			// cout << "left Intersects " << leftFruEdgeIntersect[0] << " , " << leftFruEdgeIntersect[1] << endl;
			// cout << "right Intersects " << rightFruEdgeIntersect[0] << " , " << rightFruEdgeIntersect[1] << endl;
			// cout << "Start & End angles " <<  startAngle << " " << endAngle << endl;
			// cout << "Left & Right angles " << leftAngle << " " << rightAngle << endl;
			// cout << "InFru " << startInFru << " " << endInFru << endl;
			// cout << "---\n";
			// ///

			if (!leftFruInEdge && !rightFruInEdge) {
				if (startInFru && endInFru) {
					Draw_Wall(edge_start, edge_end, cell->edges[i]->color, projection_matrix, modelview_matrix);
				} else {
					continue;
				}
			}

			if (leftFruInEdge && rightFruInEdge) {
				// if (!startInFru && !endInFru) {
					Draw_Wall(leftFruEdgeIntersect, rightFruEdgeIntersect, cell->edges[i]->color, projection_matrix, modelview_matrix);
				// }
			}
			if (leftFruInEdge && !rightFruInEdge) {
				if (endInFru) {
					Draw_Wall(leftFruEdgeIntersect, edge_end, cell->edges[i]->color, projection_matrix, modelview_matrix);
				} else if (startInFru) {
					Draw_Wall(edge_start, leftFruEdgeIntersect, cell->edges[i]->color, projection_matrix, modelview_matrix);
				} else {
					Draw_Wall(edge_start, edge_end, cell->edges[i]->color, projection_matrix, modelview_matrix);
				}
			}
			if (!leftFruInEdge && rightFruInEdge) {
				if (startInFru) {
					Draw_Wall(edge_start, rightFruEdgeIntersect, cell->edges[i]->color, projection_matrix, modelview_matrix);
				} else if (endInFru) {
					Draw_Wall(rightFruEdgeIntersect, edge_end, cell->edges[i]->color, projection_matrix, modelview_matrix);
				} else {
					Draw_Wall(edge_start, edge_end, cell->edges[i]->color, projection_matrix, modelview_matrix);
				}
			}			
		} else {
			Cell* next = cell->edges[i]->Neighbor(cell);
			if (next == NULL) continue;
			if (!next->visited) {
				Draw_Cell(next, projection_matrix, modelview_matrix, leftFru, rightFru);
			}
		}
	}
}

void Maze::
Draw_Wall(const float start[2], const float end[2], const float color[3], float* projection_matrix, float* modelview_matrix) {
	float edge0[3] = { start[Y], 0.0f, start[X] };
	float edge1[3] = { end[Y], 0.0f, end[X] };

	Matrix4 p(projection_matrix);
	Matrix4 m(modelview_matrix);

	Vector4 v1(start[Y], 1.0f, start[X], 1.0f);
	Vector4 v2(end[Y], 1.0f, end[X], 1.0f);
	Vector4 v3(end[Y], -1.0f, end[X], 1.0f);
	Vector4 v4(start[Y], -1.0f, start[X], 1.0f);
	
	p *= m;
	v1 = p * v1;
	v2 = p * v2;
	v3 = p * v3;
	v4 = p * v4;

	if (v1.w < 0 && v4.w < 0) return;
	if (v2.w < 0 && v3.w < 0) return;

	glBegin(GL_POLYGON);
	glColor3fv(color);

	v1 /= v1.w;
	v2 /= v2.w;
	v3 /= v3.w;
	v4 /= v4.w;

	// cout << v1 << v2 << v3 << v4 << endl;

	glVertex2f(v1.x, v1.y);
	glVertex2f(v2.x, v2.y);
	glVertex2f(v3.x, v3.y);
	glVertex2f(v4.x, v4.y);
	// glVertex3f(v1.x, v1.y, v1.z);
	// glVertex3f(v2.x, v2.y, v2.z);
	// glVertex3f(v3.x, v3.y, v3.z);
	// glVertex3f(v4.x, v4.y, v4.z);
	// glVertex4f(v1.x, v1.y, v1.z, v1.w);
	// glVertex4f(v2.x, v2.y, v2.z, v2.w);
	// glVertex4f(v3.x, v3.y, v3.z, v3.w);
	// glVertex4f(v4.x, v4.y, v4.z, v4.w);
	// glVertex3f(edge0[X], 1.0f, edge0[Z]);
	// glVertex3f(edge1[X], 1.0f, edge1[Z]);
	// glVertex3f(edge1[X], -1.0f, edge1[Z]);
	// glVertex3f(edge0[X], -1.0f, edge0[Z]);
	glEnd();
}

int Maze::getCode(const Vector4& v) const {
	int code = INSIDE;

	if (v.x < x_min) {
		code |= LEFT;
	} else if (v.x > x_max) {
		code |= RIGHT;
	}
	if (v.y < y_min) {
		code |= BOTTOM;
	} else if (v.y > y_max) {
		code |= TOP;
	}
	if (v.z < z_min) {
		code |= BEHIND;
	} else if (v.z > z_max) {
		code |= FRONT;
	}

	return code;
}

void Maze::Clip(Vector4& start, Vector4& end, Vector4& border) {
	int code1 = getCode(start);
	int code2 = getCode(end);

	bool accept = false;
	while(true) {
		if ((code1 == 0) && (code2 == 0)) {
			accept = true;
			break;
		} else if (code1 & code2) {
			break;
		}
		else {
			int code_out;
			double x, y, z;

			if (code1 != 0) {
				code_out = code1;
			} else {
				code_out = code2;
			}

			if (code_out & TOP) {
				x = start.x + (end.x - start.x) * (y_max - start.y) / (end.y - start.y);
				y = y_max;
				z = start.z + (end.z - start.z) * (y_max - start.y) / (end.y - start.y);
			} else if (code_out & BOTTOM) {
				x = start.x + (end.x - start.x) * (y_min - start.y) / (end.y - start.y);
				y = y_min;
				z = start.z + (end.z - start.z) * (y_min - start.y) / (end.y - start.y);
			} else if (code_out & RIGHT) {
				x = x_max;
				y = start.y + (end.y - start.y) * (x_max - start.x) / (end.x - start.x);
				z = start.z + (end.z - start.z) * (x_max - start.x) / (end.x - start.x);
			} else if (code_out & LEFT) {
				x = x_min;
				y = start.y + (end.y - start.y) * (x_min - start.x) / (end.x - start.x);
				z = start.z + (end.z - start.z) * (x_min - start.x) / (end.x - start.x);
			} else if (code_out & BEHIND) {

			} else if (code_out & FRONT) {

			}

			if (code_out == code1) {
				start.x = x;
				start.y = y;
				start.z = z;
				code1 = getCode(start);
			} else {
				end.x = x;
				end.y = y;
				end.z = z;
				code2 = getCode(end);
			}
		}
	}
}


//**********************************************************************
//
// * Draws the frustum on the map view of the maze. It is passed the
//   minimum and maximum corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Frustum(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	  height;
	float   scale_x, scale_y, scale;
	float   view_x, view_y;

	// Draws the view frustum in the map. Sets up all the same viewing
	// parameters as draw().
	scale_x	= ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y	= ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale		= scale_x > scale_y ? scale_y : scale_x;
	height	= (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	view_x = ( viewer_posn[X] - min_xp ) * scale;
	view_y = ( viewer_posn[Y] - min_yp ) * scale;
	fl_line(min_x + (int)floor(view_x + 
			  cos(To_Radians(viewer_dir+viewer_fov / 2.0)) * scale),
			  min_y + height- 
			  (int)floor(view_y + 
							 sin(To_Radians(viewer_dir+viewer_fov / 2.0)) * 
							 scale),
				min_x + (int)floor(view_x),
				min_y + height - (int)floor(view_y));
	fl_line(min_x + (int)floor(view_x + 
										cos(To_Radians(viewer_dir-viewer_fov / 2.0))	* 
										scale),
				min_y + height- 
				(int)floor(view_y + sin(To_Radians(viewer_dir-viewer_fov / 2.0)) *
				scale),
				min_x + (int)floor(view_x),
				min_y + height - (int)floor(view_y));
	}


//**********************************************************************
//
// * Draws the viewer's cell and its neighbors in the map view of the maze.
//   It is passed the minimum and maximum corners of the window in which
//   to draw.
//======================================================================
void Maze::
Draw_Neighbors(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i, j;

	// Draws the view cell and its neighbors in the map. This works
	// by drawing just the neighbor's edges if there is a neighbor,
	// otherwise drawing the edge. Every edge is shared, so drawing the
	// neighbors' edges also draws the view cell's edges.

	scale_x = ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y = ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	for ( i = 0 ; i < 4 ; i++ )   {
		Cell	*neighbor = view_cell->edges[i]->Neighbor(view_cell);

		if ( neighbor ){
			for ( j = 0 ; j < 4 ; j++ ){
				Edge    *e = neighbor->edges[j];

				if ( e->opaque )	{
					float   x1, y1, x2, y2;

					x1 = e->endpoints[Edge::START]->posn[Vertex::X];
					y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
					x2 = e->endpoints[Edge::END]->posn[Vertex::X];
					y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

					fl_color((unsigned char)floor(e->color[0] * 255.0),
							  (unsigned char)floor(e->color[1] * 255.0),
							  (unsigned char)floor(e->color[2] * 255.0));
					fl_line_style(FL_SOLID);
					fl_line( min_x + (int)floor((x1 - min_xp) * scale),
							 min_y + height - (int)floor((y1 - min_yp) * scale),
							 min_x + (int)floor((x2 - min_xp) * scale),
							 min_y + height - (int)floor((y2 - min_yp) * scale));
				}
			}
		}
		else {
			Edge    *e = view_cell->edges[i];

			if ( e->opaque ){
				float   x1, y1, x2, y2;

				x1 = e->endpoints[Edge::START]->posn[Vertex::X];
				y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
				x2 = e->endpoints[Edge::END]->posn[Vertex::X];
				y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

				fl_color((unsigned char)floor(e->color[0] * 255.0),
							 (unsigned char)floor(e->color[1] * 255.0),
							 (unsigned char)floor(e->color[2] * 255.0));
				fl_line_style(FL_SOLID);
				fl_line(min_x + (int)floor((x1 - min_xp) * scale),
							min_y + height - (int)floor((y1 - min_yp) * scale),
							min_x + (int)floor((x2 - min_xp) * scale),
							min_y + height - (int)floor((y2 - min_yp) * scale));
			 }
		}
	}
}


//**********************************************************************
//
// * Save the maze to a file of the given name.
//======================================================================
bool Maze::
Save(const char *filename)
//======================================================================
{
	FILE    *f = fopen(filename, "w");
	int	    i;

	// Dump everything to a file of the given name. Returns false if it
	// couldn't open the file. True otherwise.

	if ( ! f )  {
		return false;
   }

	fprintf(f, "%d\n", num_vertices);
	for ( i = 0 ; i < num_vertices ; i++ )
		fprintf(f, "%g %g\n", vertices[i]->posn[Vertex::X],
			      vertices[i]->posn[Vertex::Y]);

		fprintf(f, "%d\n", num_edges);
	for ( i = 0 ; i < num_edges ; i++ )
	fprintf(f, "%d %d %d %d %d %g %g %g\n",
				edges[i]->endpoints[Edge::START]->index,
				edges[i]->endpoints[Edge::END]->index,
				edges[i]->neighbors[Edge::LEFT] ?
				edges[i]->neighbors[Edge::LEFT]->index : -1,
				edges[i]->neighbors[Edge::RIGHT] ?
				edges[i]->neighbors[Edge::RIGHT]->index : -1,
				edges[i]->opaque ? 1 : 0,
				edges[i]->color[0], edges[i]->color[1], edges[i]->color[2]);

	fprintf(f, "%d\n", num_cells);
	for ( i = 0 ; i < num_cells ; i++ )
		fprintf(f, "%d %d %d %d\n",
					cells[i]->edges[0] ? cells[i]->edges[0]->index : -1,
					cells[i]->edges[1] ? cells[i]->edges[1]->index : -1,
					cells[i]->edges[2] ? cells[i]->edges[2]->index : -1,
					cells[i]->edges[3] ? cells[i]->edges[3]->index : -1);

	   fprintf(f, "%g %g %g %g %g\n",
					viewer_posn[X], viewer_posn[Y], viewer_posn[Z],
					viewer_dir, viewer_fov);

	fclose(f);

	return true;
}