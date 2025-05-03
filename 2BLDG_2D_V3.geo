SetFactory("OpenCASCADE");

//Constants
a = 40; b = 100; c =25; //Building Size params
delta = 10; //distnace between buildings
L = delta+2*b+80; W = delta+2*a+50; W_offset = a+20; H = c+15; //Domain Size Params

h1 = a/4;
h2 = 0.25*h1;
//--------------------------
// CITY DOMAIN DEFINITION
//--------------------------

// Outer Domain Points 
Point(1) = {-L*.25, -L*.25, 0, h1};
Point(2) = {W, -L*.25, 0, h1};
Point(3) = {W, L*1.25, 0, h1};
Point(4) = {-L*.25, L*1.25, 0, h1};

// Outer domain Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


// outer city surfaces
Curve Loop(21) = {1, 2, 3, 4};         
Plane Surface(31) = {21}; // Bottom

//--------------------------
// BUILDINGS DEFFINITIONS
//--------------------------
//Building 1
x0 = (W_offset - a)/2; x1 = (W_offset + a)/2;
y0 = (L - delta)/2 - b; y1 = (L - delta)/2;

Point(101) = {x0, y0, 0, h2};
Point(102) = {x1, y0, 0, h2};
Point(103) = {x1, y1, 0, h2};
Point(104) = {x0, y1, 0, h2};

Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 101};


// Surfaces
Line Loop(201) = {101, 102, 103, 104};     
Plane Surface(301) = {201}; // building 1


//Building 2
x0 = (W_offset - a)/2; x1 = (W_offset + a)/2;
y0 = (L + delta)/2 + b; y1 = (L + delta)/2;

Point(109) = {x0, y0, 0, h2};
Point(110) = {x1, y0, 0, h2};
Point(111) = {x1, y1, 0, h2};
Point(112) = {x0, y1, 0, h2};

Line(149) = {109, 112};
Line(150) = {112, 111};
Line(151) = {111, 110};
Line(152) = {110, 109};

// Surfaces
Line Loop(207) = {149,150,151,152};
Plane Surface(307) = {207}; //building 2


// ------------------------------
// Subtract buildings from domain
// ------------------------------
BooleanDifference(100) = { Surface{31}; Delete; } { Surface{301, 307}; Delete; };


// ----------------------------
// Physical Groups 
// ----------------------------
Physical Surface(1) = {100}; // domain with building removed 

//Physical Surface(2) = {301}; //building 1

//Physical Surface(3) = {307}; //building 2

// Physical boundaries for BCs
Physical Line(101) = {3};              // Outlet (South)
Physical Line(102) = {4};              // East
Physical Line(103) = {2};              // Inlet (North)
Physical Line(104) = {1};              // West
Physical Line(200) = {5, 6, 7, 8, 9, 10, 11, 12}; // Buildings (shared group)


// Physical boundaries for BCs
//Physical Line(101) = {155};              // Outlet (South)
//Physical Line(102) = {156};              // East
//Physical Line(103) = {154};              // Inlet (North)
//Physical Line(104) = {153};              // West
//Physical Line(200) = {101,102,103,104,149,150,151,152}; // Buildings (shared group)

// Mesh control
Mesh.CharacteristicLengthMax = h1;