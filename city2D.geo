SetFactory("OpenCASCADE");

// -------------------
// CONSTANTS
// -------------------
Lx = 350; Ly = 150; // Domain size
h = 5; h1 = 0.5*h;  // Mesh size

// -------------------
// DOMAIN
// -------------------
Point(1) = {-Lx/3.5, -40, 0, h1};
Point(2) = {Lx, -40, 0, h1};
Point(3) = {Lx, Ly, 0, h1};
Point(4) = {-Lx/3.5, Ly, 0, h1};

Line(1) = {1, 2};  // Bottom
Line(2) = {2, 3};  // Right
Line(3) = {3, 4};  // Top
Line(4) = {4, 1};  // Left

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1}; // Main domain

// -------------------
// BUILDING 1 FOOTPRINT
// -------------------
Point(5) = {10, 10, 0, h1};
Point(6) = {10, 55, 0, h1};
Point(7) = {40, 55, 0, h1};
Point(8) = {40, 30, 0, h1};
Point(9) = {50, 30, 0, h1};
Point(10) = {50, 10, 0, h1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 5};

Curve Loop(2) = {5, 6, 7, 8, 9, 10};
Plane Surface(2) = {2};

// -------------------
// BUILDING 2 FOOTPRINT
// -------------------
Point(11) = {10, 62, 0, h1};
Point(12) = {10, 107, 0, h1};
Point(13) = {101, 107, 0, h1};
Point(14) = {101, 62, 0, h1};

Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};

Curve Loop(3) = {11, 12, 13, 14};
Plane Surface(3) = {3};

// -------------------
// BUILDING 3 FOOTPRINT
// -------------------
Point(15) = {56, 10, 0, h1};
Point(16) = {56, 55, 0, h1};
Point(17) = {101, 55, 0, h1};
Point(18) = {101, 10, 0, h1};

Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 15};

Curve Loop(4) = {15, 16, 17, 18};
Plane Surface(4) = {4};

// -------------------
// BUILDING 4 FOOTPRINT
// -------------------
Point(19) = {126, 69, 0, h1};
Point(20) = {126, 107, 0, h1};
Point(21) = {216, 107, 0, h1};
Point(22) = {216, 69, 0, h1};

Line(19) = {19, 20};
Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 19};

Curve Loop(5) = {19, 20, 21, 22};
Plane Surface(5) = {5};

// -------------------
// BUILDINGS 5 & 6 FOOTPRINT (shared face)
// -------------------
Point(23) = {126, 10, 0, h1};
Point(24) = {126, 28, 0, h1};
Point(25) = {146, 28, 0, h1};
Point(26) = {146, 43, 0, h1};
Point(27) = {156, 43, 0, h1};
Point(28) = {156, 10, 0, h1};
Point(29) = {156, 50, 0, h1};
Point(30) = {216, 50, 0, h1};
Point(31) = {216, 10, 0, h1};

Line(23) = {23, 24};
Line(24) = {24, 25};
Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 29};
Line(28) = {29, 30};
Line(29) = {30, 31};
Line(30) = {31, 28};
Line(31) = {28, 23};

Curve Loop(6) = {23, 24, 25, 26, 27, 28, 29, 30, 31};
Plane Surface(6) = {6};

// -------------------
// BOOLEAN SUBTRACTION
// -------------------
BooleanDifference{ Surface{1}; Delete; }{ Surface{2,3,4,5,6}; Delete; }

// -------------------
// PHYSICAL GROUPS
// -------------------
Physical Surface(1) = {1}; // Fluid domain w/ building removed 

// Physical boundaries for BCs
Physical Line(101) = {35};              // Outlet (South)
Physical Line(102) = {33};              // East
Physical Line(103) = {32};              // Inlet (North)
Physical Line(104) = {34};              // West
Physical Line(200) = {5,6,7,8,9,10, 11,12,13,14, 15,16,17,18, 19,20,21,22, 23,24,25,26,27,28,29,30,31};

// Mesh control
Mesh.CharacteristicLengthMax = h1;
