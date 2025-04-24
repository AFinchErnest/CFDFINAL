//Mesh file for a single building as a rectangular prisim
SetFactory("OpenCASCADE");

//Constants
L = 20; //long length of surround
W = 10; //short length of surroundings
H = 10; //height of surroundings

a = 2; //side length of building
b = 6; //height of building

h1 = 1.0; //outer mesh size
h2 = 0.5*h1; //inner mesh size

//Surroundings Lines and Points
Point(1) = {0, 0, 0, h1};
Point(2) = {W, 0, 0, h1};
Point(3) = {W, L, 0, h1};
Point(4) = {0, L, 0, h1};
Point(5) = {0, 0, H, h1};
Point(6) = {W, 0, H, h1};
Point(7) = {W, L, H, h1};
Point(8) = {0, L, H, h1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

//Building Lines and Points
Point(9) = {(W-a)/2,(W-a)/2,0,h2};
Point(10) = {(W-a)/2,(W+a)/2,0,h2};
Point(11) = {(W+a)/2,(W+a)/2,0,h2};
Point(12) = {(W+a)/2,(W-a)/2,0,h2};
Point(13) = {(W-a)/2,(W-a)/2,b,h2};
Point(14) = {(W-a)/2,(W+a)/2,b,h2};
Point(15) = {(W+a)/2,(W+a)/2,b,h2};
Point(16) = {(W+a)/2,(W-a)/2,b,h2};

Line(13) = {9,10};
Line(14) = {10,11};
Line(15) = {11,12};
Line(16) = {12,9};
Line(17) = {13,14};
Line(18) = {14,15};
Line(19) = {15,16};
Line(20) = {16,13};
Line(21) = {9,13};
Line(22) = {10,14};
Line(23) = {11,15};
Line(24) = {12,16};

//City South
Curve Loop(1) = {9, 5, -10, -1};
Plane Surface(1) = {1};

Physical Surface(25) = {1};

//City North
Curve Loop(2) = {11, 7, -12, -3};
Plane Surface(2) = {2};

Physical Surface(26) = {2};

//City East
Curve Loop(3) = {4, 9, -8, -12};
Plane Surface(3) = {3};

Physical Surface(27) = {3};

//City West
Curve Loop(4) = {6, -11, -2, 10};
Plane Surface(4) = {4};

Physical Surface(28) = {4};

//City Sky
Curve Loop(5) = {8, 5, 6, 7};
Plane Surface(5) = {5};

Physical Surface(29) = {5};

//City Ground
Curve Loop(6) = {2, 3, 4, 1};
Curve Loop(7) = {15, 16, 13, 14}; //floor of building
Plane Surface(6) = {6, 7};

Physical Surface(30) = {6};

//Building 1
Curve Loop(8) = {16, 21, -20, -24}; //South
Plane Surface(7) = {8}; //South
Curve Loop(9) = {15, 24, -19, -23}; // East
Plane Surface(8) = {9}; //East
Curve Loop(10) = {14, 23, -18, -22}; //North
Plane Surface(9) = {10}; //North
Curve Loop(11) = {13, 22, -17, -21}; //West
Plane Surface(10) = {11};//West
Curve Loop(12) = {17, 18, 19, 20}; //roof
Plane Surface(11) = {12}; //roof

Physical Surface(31) = {7};
Physical Surface(32) = { 8};
Physical Surface(33) = { 9};
Physical Surface(34) = { 10};
Physical Surface(35) = { 11};

//City Domain
Surface Loop(1) = {4, 5, 3, 6, 2, 1, 9, 8, 7, 10, 11};
Volume(1) = {1};
Physical Volume(1) = {1};
