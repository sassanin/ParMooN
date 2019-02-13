// Gmsh project created on Thu Jun 18 14:38:14 2015

//Modified: 22-July-2015
//gmsh file to generate 2D model



lc = 0.50;
crad = 8.50;

//  calculate the angle for the circle, needed in PRM file
//     cout<< " rad  "<< atan2(2.2, (0.289641-crad)) << " Angle  "<<(180/Pi)*atan2(2.2, (0.289641-crad))<<endl;
//     cout<< " x  "<<crad+ crad*cos(atan2(2.2, (0.289641-crad))) << " y  "<<crad*sin(atan2(2.2, (0.289641-crad)))<<endl;
//  for crad = 8.50 ==> rad  2.87979 Angle  165


Point(1) = {0, 0, 0, lc};
Point(2) = {crad,0,0,lc};
Point(3) = {15.00,  0.00, 0, lc};
Point(4) = {15.00,  1.25, 0, lc};
Point(5) = {12.90,  1.25, 0, lc};
Point(6) = {12.00,  1.25, 0, lc};
Point(7) = { 6.30,  1.50, 0, lc};
Point(8) = { 6.30,  3.50, 0, lc};
Point(9) = {10.70,  4.50, 0, lc};
Point(10) = {10.70, 10.0, 0, lc};
Point(11) = {0.00, 10.0, 0, lc};
Point(12) = {0.00, 5.70, 0, lc};
Point(13) = {2.05, 5.70, 0, lc};
Point(14) = {2.05, 3.70, 0, lc};
Point(15) = {1.50, 3.70, 0, lc};
Point(16) = {0.50, 2.50, -0, lc};
Point(17) = {0, 2.50, 0, lc};
Point(18) = {0, 2.20, 0, lc};
Point(19) = {0.289641, 2.20, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};

Circle(19) = {1, 2, 19};

//add plane surface using lines
Line Loop(20) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, -19};


Physical Line(1000) = {1}; //Straight line
Physical Line(1001) = {2}; //Straight line
Physical Line(1002) = {3}; //Straight line
Physical Line(1003) = {4}; //Straight line
Physical Line(1004) = {5}; //Straight line
Physical Line(1005) = {6}; //Straight line
Physical Line(1006) = {7}; //Straight line
Physical Line(1007) = {8}; //Straight line
Physical Line(1008) = {9}; //Straight line
Physical Line(1009) = {10}; //Straight line
Physical Line(1010) = {11}; //Straight line
Physical Line(1011) = {12}; //Straight line
Physical Line(1012) = {13}; //Straight line
Physical Line(1013) = {14}; //Straight line
Physical Line(1014) = {15}; //Straight line
Physical Line(1015) = {16}; //Straight line
Physical Line(1016) = {17}; //Straight line
Physical Line(1017) = {18}; //Straight line

Physical Line(2000) = {19};   //Curved line

Plane Surface(21) = {20};
Physical Surface(5000) = {21};

