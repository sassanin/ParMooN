lc1 = 0.5;
zstart = -2;
zend = 2;
xstart = -2;
xend = 2;
ystart = -2;
yend = 2;

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

Physical Surface(1000)={1};
Physical Surface(1001)={2};
Physical Surface(1002)={3};
Physical Surface(1003)={4};
Physical Surface(1004)={5};
Physical Surface(1005)={6};

Surface Loop(1)={1:6};
//========================================================
lc1=1;
rad=1;
j=1000;
Point(1+j) = {0,rad,0,lc1}; //top
Point(2+j) = {0,-rad,0,lc1}; //bottom
Point(3+j) = {-rad,0,0,lc1}; //left
Point(4+j) = {rad,0,0,lc1}; //right
Point(5+j) = {0,0,rad,lc1}; //front
Point(6+j) = {0,0,-rad,lc1}; //back
Point(7+j) = {0,0,0,lc1}; //centre

Circle(1+j) = {5+j,7+j,4+j};
Circle(2+j) = {4+j,7+j,6+j};
Circle(3+j) = {6+j,7+j,3+j};
Circle(4+j) = {3+j,7+j,5+j};
Circle(5+j) = {5+j,7+j,2+j};
Circle(6+j) = {2+j,7+j,6+j};
Circle(7+j) = {6+j,7+j,1+j};
Circle(8+j) = {1+j,7+j,5+j};
Circle(9+j) = {4+j,7+j,1+j};
Circle(10+j) = {1+j,7+j,3+j};
Circle(11+j) = {3+j,7+j,2+j};
Circle(12+j) = {2+j,7+j,4+j};

Line Loop(1+j)={1+j,9+j,8+j};
Line Loop(2+j)={2+j,7+j,-(9+j)};
Line Loop(3+j)={3+j,-(10+j),-(7+j)};
Line Loop(4+j)={4+j,-(8+j),10+j};
Line Loop(5+j)={1+j,-(12+j),-(5+j)};
Line Loop(6+j)={2+j,-(6+j),12+j};
Line Loop(7+j)={3+j,11+j,6+j};
Line Loop(8+j)={4+j,5+j,-(11+j)};
Ruled Surface(1+j)={1+j};
Ruled Surface(2+j)={2+j};
Ruled Surface(3+j)={3+j};
Ruled Surface(4+j)={4+j};
Ruled Surface(5+j)={5+j};
Ruled Surface(6+j)={6+j};
Ruled Surface(7+j)={7+j};
Ruled Surface(8+j)={8+j};

Surface Loop(1+j) = {1+j:8+j};
Physical Surface(2000)={1+j:8+j};

Volume(1+j)={1+j};         //sphere volume
Physical Volume(0)={1+j};  //tag for sphere volume

Volume(1)={1,1+j};        //cube volume without the sphere
Physical Volume(1)={1};   //tag for the cube volume without the sphere 
//================================
View "comments" {
  // Add a text string in window coordinates, 10 pixels from the left
  // and 10 pixels from the bottom:
  T2(10, -10, 0){ "Copyright (C) NMSC" };
  // Add another text string in window coordinates, 10 pixels from the
  // left and 15 pixels from the top, using the StrCat() function to
  // concatenate a string with the current date:
  T2(10, 15, 0){ StrCat("File created on ", Today) };
};
