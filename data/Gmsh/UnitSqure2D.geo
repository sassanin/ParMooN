// Gmsh project created on Thu Jan 26 2016

//gmsh file to generate 2D model

lc = 0.50;
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1,0,lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};

Physical Line(1000) = {1};  //Straight line
Physical Line(1001) = {2};  //Straight line
Physical Line(1002) = {3};  //Straight line
Physical Line(1003) = {4};  //Straight line

Plane Surface(6) = {5};
Physical Surface(5000) = {6};
