lc1 = 1.5;
zstart = 0;
zend = 1;
xstart = 0;
xend = 1;
ystart = 0;
yend = 1;

Point(1) = {xstart, ystart, zstart, lc1};
Point(2) = {xend, ystart, zstart, lc1};
Point(3) = {xend, yend, zstart, lc1};
Point(4) = {xstart, yend, zstart, lc1};

Point(5) = {xstart, ystart, zend, lc1};
Point(6) = {xend, ystart, zend, lc1} ;
Point(7) = {xend, yend, zend, lc1} ;
Point(8) = {xstart, yend, zend, lc1} ;

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,5};
Line(9)={1,5};
Line(10)={2,6};
Line(11)={3,7};
Line(12)={4,8};

Line Loop(1)={9,-8,-12,4}; //left
Line Loop(2)={5,-10,-1,9}; //bottom
Line Loop(3)={10,6,-11,-2}; //right
Line Loop(4)={-7,-11,3,12}; //top
Line Loop(5)={5,6,7,8}; //front
Line Loop(6)={1,2,3,4}; //back

Plane Surface(1)={1};
Plane Surface(2)={2};
Plane Surface(3)={3};
Plane Surface(4)={4};
Plane Surface(5)={5};
Plane Surface(6)={6};

Surface Loop(1)={1:6};

Physical Surface(1000)={1};
Physical Surface(1001)={2};
Physical Surface(1002)={3};
Physical Surface(1003)={4};
Physical Surface(1004)={5};
Physical Surface(1005)={6};

Volume(1)={1};
Physical Volume(1)={1};

// To generate quadrangles instead of triangles, we can simply add
//Mesh.RecombineAll =1;
//Mesh.Recombine3DAll =1;
