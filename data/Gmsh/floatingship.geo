
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
fluback=1;
Line Loop(fluback) = {1:4};
Plane Surface(fluback) = {fluback};
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
//==================================================================
Line(newl) = {1,pend1+1};
Line(newl) = {2,pend1+2};
Line(newl) = {3,pend1+3};
Line(newl) = {pend1,pend2};
//===============================
Line Loop(3) = {9,-8,-12,4};
Plane Surface(3) = {3};
Line Loop(4) = {9,5,-10,-1};
Plane Surface(4) = {4};
Line Loop(5) = {10,6,-11,-2};
Plane Surface(5) = {5};
Line Loop(6) = {-7,-11,3,12};
// Ruled Surface(6) = {6};
// //=================================
ybottom = 0.7*ymax;
ytop = 1.2*ymax;
xleft = 0.3*xmax;
xright = 0.7*xmax;
zback = 0.3*zmax;
zfront = 0.7*zmax;

ps1 = newp;
Point(ps1) = {xleft,ybottom,zback,lc1};
Point(ps1+1) = {xright,ybottom,zback,lc1};

For t In {xright:xleft:-0.1}
Point(newp) = {t,ymax+wavh*Sin(2*3.141592653589793*t/xmax),zback,lcspline};
EndFor
psend1 = newp-1;

sl1 = newl;
Line(sl1) = {ps1,ps1+1};
Line(sl1+1) = {ps1+1,ps1+2};
Spline(sl1+2) = {ps1+2:psend1};
Line(sl1+3) = {psend1,ps1};

Line Loop(7) = {sl1:sl1+3};
Plane Surface(7) = {7};
//=================================
ps2 = newp;
Point(ps2) = {xleft,ybottom,zfront,lc1};
Point(ps2+1) = {xright,ybottom,zfront,lc1};

For t In {xright:xleft:-0.1}
Point(newp) = {t,ymax+wavh*Sin(2*3.141592653589793*t/xmax),zfront,lcspline};
EndFor
psend2 = newp-1;

sl2 = newl;
Line(sl2) = {ps2,ps2+1};
Line(sl2+1) = {ps2+1,ps2+2};
Spline(sl2+2) = {ps2+2:psend2};
Line(sl2+3) = {psend2,ps2};

Line Loop(8) = {sl2:sl2+3};
Plane Surface(8) = {8};
// //==============================================
sline = newl;
Line(sline) = {ps1,ps2};
Line(sline+1) = {ps1+1,ps2+1};
Line(sline+2) = {ps1+2,ps2+2};
Line(sline+3) = {psend1,psend2};

shback = newll;
Line Loop(shback) = {sline,-(sl2+3),-(sline+3),sl1+3};
Plane Surface(shback) = {shback};
shbottom = newll;
Line Loop(shbottom) = {sline,sl2,-(sline+1),-sl1};
Plane Surface(shbottom) = {shbottom};
shfront = newll;
Line Loop(shfront) = {sline+1,sl2+1,-(sline+2),-(sl1+1)};
Plane Surface(shfront) = {shfront};

fluidcrossh = newll;
Line Loop(fluidcrossh) = {sl1+2,sline+3,-(sl2+2),-(sline+2)};
Ruled Surface(6) = {6,fluidcrossh};

Surface Loop(1) = {26,8,27,7,25,6,3,2,5,1,4};
Volume(1) = {1} ;
//Physical Volume(250) = {1};
//=============================================
ybottom = 0.7*ymax;
ytop = 1.2*ymax;
xleft = 0.3*xmax;
xright = 0.7*xmax;
zback = 0.3*zmax;
zfront = 0.7*zmax;
ymaxship = ymax+2;

stopp1 = newp;
Point(stopp1) = {xleft,ymaxship,zback,lc1};
Point(stopp1+1) = {xleft,ymaxship,zfront,lc1};
Point(stopp1+2) = {xright,ymaxship,zfront,lc1};
Point(stopp1+3) = {xright,ymaxship,zback,lc1};
topl = newl;
Line(topl) = {stopp1,stopp1+1};
Line(topl+1) = {stopp1+1,stopp1+2};
Line(topl+2) = {stopp1+2,stopp1+3};
Line(topl+3) = {stopp1+3,stopp1};
corners = newl;
Line(corners) = {stopp1,psend1};
Line(corners+1) = {stopp1+1,psend2};
Line(corners+2) = {stopp1+2,ps2+2};
Line(corners+3) = {stopp1+3,ps1+2};

shiptop = newll;
Line Loop(shiptop) = {topl:topl+3};
Plane Surface(shiptop) = {shiptop};
shipback=newll;
Line Loop(shipback) = {-topl,corners,sline+3,-(corners+1)};
Plane Surface(shipback) = {shipback};
shipforward = newll;
Line Loop(shipforward) = {corners+1,-(sl2+2),-(corners+2),-(topl+1)};
Plane Surface(shipforward) = {shipforward};
shipfront=newll;
Line Loop(shipfront) = {topl+2,corners+3,sline+2,-(corners+2)};
Plane Surface(shipfront) = {shipfront};
shipbehind=newll;
Line Loop(shipbehind) = {topl+3,corners,-(sl1+2),-(corners+3)};
Plane Surface(shipbehind) = {shipbehind};
 
Surface Loop(2) = {shiptop,shipback,shipforward,shipfront,shipbehind,7,shback,8,shfront,shbottom};
Volume(2) = {2};
//Physical Volume(260) = {2};


Physical Surface(1000)={1};
Physical Surface(1001)={2};
Physical Surface(1002)={3};
Physical Surface(1003)={4};
Physical Surface(1004)={5};
Physical Surface(1005)={6};
Physical Surface(1006)={7};
Physical Surface(1007)={8};
Physical Surface(1008)={shback};
Physical Surface(1009)={shbottom};
Physical Surface(1010)={shfront};
Physical Surface(1011)={shipbehind};
Physical Surface(1012)={shipback};
Physical Surface(1013)={shipforward};
Physical Surface(1014)={shipfront};
Physical Surface(1015)={shiptop};

Physical Volume(0) = {1};
Physical Volume(1) = {2};


// Color Red{Surface{21};}
// Field[1] = Attractor;
// Field[1].NodesList = {5};
// Field[1].NNodesByEdge = 100;
// Field[1].EdgesList = {10};
// 
// Field[2] = Threshold;
// Field[2].IField = 1;
// Field[2].LcMin = lc /15;
// Field[2].LcMax = lc;
// Field[2].DistMin = 0.05;
// Field[2].DistMax = 0.5;
// 
// Field[3] = Min;
// Field[3].FieldsList = {2};
// Background Field = 3;
// 
// Color Grey50{ Surface{ 21 }; }
// Color Red{ Surface{ 23 }; }
// Color Red{ Line{ 1:14 }; }
// Color Yellow{ Line{ 15:20 }; }

