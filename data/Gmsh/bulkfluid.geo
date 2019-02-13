
lc1 = 3;
lcspline = 3;
xmax = 20;
ymax = 6;
zmax = 8;
// xshift = -0.2;
// yshift = 0.8;
wavh = 0.8;

Point(1) = {0,0,0,lc1};
Point(2) = {xmax,0,0,lc1};
For t In {xmax:0:-0.1}
Point(newp) = {t,ymax+wavh*Sin(2*3.141592653589793*t/xmax),0,lcspline};
EndFor
pend1=newp-1;
Line(1) = {1,2};
Line(2) = {2,3};
Spline(3) = {3:pend1};
Line(4) = {pend1,1};
Line Loop(1) = {1:4};
Plane Surface(1) = {1};
//fluback
//===================================================================
p1=newp;
Point(p1) = {0,0,zmax,lc1};
Point(newp) = {xmax,0,zmax,lc1};
For t In {xmax:0:-0.1}
Point(newp) = {t,ymax+wavh*Sin(2*3.141592653589793*t/xmax),zmax,lcspline};
EndFor
pend2 = newp-1;
l1=newl;
Line(l1) = {p1,p1+1};
Line(newl) = {p1+1,p1+2};
Spline(newl) = {p1+2:pend2};
l2=newl;
Line(l2) = {pend2,p1};
Line Loop(2) = {l1:l2};
Plane Surface(2) = {2};
//flufront
//==================================================================
Line(newl) = {1,pend1+1};
Line(newl) = {2,pend1+2};
Line(newl) = {3,pend1+3};
Line(newl) = {pend1,pend2};
//===============================
Line Loop(3) = {9,-8,-12,4};
Plane Surface(3) = {3};  //fluleft
Line Loop(4) = {9,5,-10,-1};
Plane Surface(4) = {4}; //flubottom
Line Loop(5) = {10,6,-11,-2};
Plane Surface(5) = {5};  //fluright
Line Loop(6) = {-7,-11,3,12};
Ruled Surface(6) = {6};  //flutop
//=================================
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};
//================================
Physical Surface(1000)={1};
Physical Surface(1001)={2};
Physical Surface(1002)={3};
Physical Surface(1003)={4};
Physical Surface(1004)={5};
Physical Surface(1005)={6};
Physical Volume(0) = {1};


