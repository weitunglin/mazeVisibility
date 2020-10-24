/************************************************************************
     File:        LineSeg.h

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison
						Class header file for LineSeg class.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#ifndef _LINESEG_H_
#define _LINESEG_H_

#include "Edge.h"

class LineSeg {

	public:
		// First constructor takes the endpoints.
		LineSeg(float xs, float ys, float xe, float ye);

		// Second constructor takes an edge. The LineSeg created has the same
		// start and end points as the edge.
		LineSeg(Edge*);

		LineSeg(const LineSeg&);

	public:
		// Return the parameter value at which this segment crosses the given
		// segment. This will return parameter values outside the range 0,1
		// THIS FUNCTION IS EXTREMELY USEFUL FOR CLIPPING, but it 
		// DOES NOT tell you whether the edge is "entering" or "leaving".
		// But you can use tests like Edge::Point_Side() to figure that out.
		float   Cross_Param(LineSeg);
		static bool get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y, 
     	  double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y);

	public:
		double   start[2];	// Starting point, x and y.
		double   end[2];	// Ending point, x and y.
};

#endif

