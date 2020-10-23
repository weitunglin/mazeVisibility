/************************************************************************
     File:        LineSeg.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison
						Class header file for LineSeg class.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "LineSeg.h"

//**********************************************************************
//
// * Constructor from an edge
//======================================================================
LineSeg::
LineSeg(Edge *e)
//======================================================================
{
	start[0] = e->endpoints[Edge::START]->posn[Vertex::X];
	start[1] = e->endpoints[Edge::START]->posn[Vertex::Y];

	end[0] = e->endpoints[Edge::END]->posn[Vertex::X];
	end[1] = e->endpoints[Edge::END]->posn[Vertex::Y];
}

//**********************************************************************
//
// * Constructor for specifyng start and end point
//======================================================================
LineSeg::
LineSeg(float xs, float ys, float xe, float ye)
//======================================================================
{
	start[0] = xs;
	start[1] = ys;
	end[0] = xe;
	end[1] = ye;
}

LineSeg::LineSeg(const LineSeg& r) {
    start[0] = r.start[0];
	start[1] = r.start[1];
	end[0] = r.end[0];
	end[1] = r.end[1];
}

//**********************************************************************
//
// * Return the parameter value at which this segment crosses the given
//   segment. This will return parameter values outside the range 0,1
//   THIS FUNCTION IS EXTREMELY USEFUL FOR CLIPPING, but it 
//   DOES NOT tell you whether the edge is "entering" or "leaving".
//   But you can use tests like Edge::Point_Side() to figure that out.
//======================================================================
float LineSeg::
Cross_Param(LineSeg e)
//======================================================================
{
	float   dx1, dy1, dx2, dy2;
	float   denom, s;

	// This computation comes from writing each segment in parametric form,
	// and solving a simulataneous equation to determine the parameter
	// value of the intersection point on this LineSeg.
	dx1 = e.end[0] - e.start[0];
	dy1 = e.end[1] - e.start[1];
	dx2 = end[0] - start[0];
	dy2 = end[1] - start[1];

	if ( ( denom = dx2 * dy1 - dy2 * dx1 ) == 0.0 )
		// Parallel segments.
		return 1.0e20f;

	s = ( e.start[0] - start[0] ) * dy1 - ( e.start[1] - start[1] ) * dx1;

	return s / denom;
}

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
bool LineSeg::get_line_intersection(float p0_x, float p0_y, float p1_x, float p1_y, 
    float p2_x, float p2_y, float p3_x, float p3_y, float *i_x, float *i_y)
{
    float s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    float s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (i_x != 0)
            *i_x = p0_x + (t * s1_x);
        if (i_y != 0)
            *i_y = p0_y + (t * s1_y);
        return true;
    }

    return false; // No collision
}