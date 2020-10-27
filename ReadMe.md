# How To Run?

## Compile the program
> cd src \
> make runmaze \
> make buildmaze

## execute the program

### runmaze
> ./runmaze

## buildmaze
> ./buildmaze

## clean up executables
> make clean

# Technical Document

## Step1 Compute Matries (projection & modelview)
after proofing the mathematics
i simplying use the result from above

## Step2 Clipping
i use the algorithm teacher gave
i have four major parameters
* leftFruInEdge (check if left view frustum intersects the edge, if yes record the intersection point)
* rightFruInEdge (check if right view frustum intersects the edge, if yes record the intersection point)
* startInEdge (check if edge's start point in view frustum, by using radians)
* endInEdge (check if edge's end point in view frustum, by using radians)
and use lots of if to decide how to draw the wall

## Step3 Visibility (Recursion)
same as above, but recursive
main idea is restrict view frustum then pass it to next cell

