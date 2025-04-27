SetFactory("OpenCASCADE");

L = 20; W = 10; H = 10;
a = 2; b = 6;

h1 = 1.0;
h2 = 0.5*h1;

// Outer Domain Points 
Point(1) = {0, 0, 0, h1};
Point(2) = {W, 0, 0, h1};
Point(3) = {W, L, 0, h1};
Point(4) = {0, L, 0, h1};
Point(5) = {0, 0, H, h1};
Point(6) = {W, 0, H, h1};
Point(7) = {W, L, H, h1};
Point(8) = {0, L, H, h1};

// Outer city Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// oter city surfaces
Curve Loop(21) = {1, 2, 3, 4};         
Plane Surface(31) = {21}; // Bottom
Curve Loop(22) = {5, 6, 7, 8};         
Plane Surface(32) = {22}; // Top
Curve Loop(23) = {1, 10, -5, -9};      
Plane Surface(33) = {23}; // South
Curve Loop(24) = {2, 11, -6, -10};     
Plane Surface(34) = {24}; // East
Curve Loop(25) = {3, 12, -7, -11};    
Plane Surface(35) = {25}; // North
Curve Loop(26) = {4, 9, -8, -12};      
Plane Surface(36) = {26}; // West

Surface Loop(40) = {31, 32, 33, 34, 35, 36};
Volume(41) = {40};

// Building
x0 = (W - a)/2; x1 = (W + a)/2;
y0 = (L - a)/2; y1 = (L + a)/2;

Point(101) = {x0, y0, 0, h2};
Point(102) = {x1, y0, 0, h2};
Point(103) = {x1, y1, 0, h2};
Point(104) = {x0, y1, 0, h2};
Point(105) = {x0, y0, b, h2};
Point(106) = {x1, y0, b, h2};
Point(107) = {x1, y1, b, h2};
Point(108) = {x0, y1, b, h2};

Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 101};

Line(105) = {105, 106};
Line(106) = {106, 107};
Line(107) = {107, 108};
Line(108) = {108, 105};

Line(109) = {101, 105};
Line(110) = {102, 106};
Line(111) = {103, 107};
Line(112) = {104, 108};

// Building Surfaces
Line Loop(201) = {101, 102, 103, 104};     
Plane Surface(301) = {201}; // Bottom
Curve Loop(202) = {105, 106, 107, 108};     
Plane Surface(302) = {202}; // Top
Curve Loop(203) = {101, 110, -105, -109};   
Plane Surface(303) = {203}; // South
Curve Loop(204) = {102, 111, -106, -110};   
Plane Surface(304) = {204}; // East
Curve Loop(205) = {103, 112, -107, -111};   
Plane Surface(305) = {205}; // North
Curve Loop(206) = {104, 109, -108, -112};   
Plane Surface(306) = {206}; // West

Surface Loop(310) = {301, 302, 303, 304, 305, 306};
Volume(311) = {310};

// subtract volumes
BooleanDifference(99)={ Volume{41}; Delete; } { Volume{311}; Delete; };

// Physical Groups 
// Domain volume


// Domain faces
//Physical Surface(51) = {35}; // North
//Physical Surface(52) = {33}; // South
//Physical Surface(53) = {36}; // West
//Physical Surface(54) = {34}; // East
//Physical Surface(55) = {32}; // Top
//Physical Surface(56) = {31}; // Ground
//Physical Surface(57) = {303, 304, 305, 306, 302}; // Building sides + roof

// Mesh control
Mesh.CharacteristicLengthMax = h1;
//+
Physical Surface(100) = {5}; //WEST
Physical Surface(101) = {3} ;//South
Physical Surface(102) = {4} ;//North
Physical Surface(103) = {2} ;//EAST
Physical Surface(104) = {10};  //SKY
Physical Surface(105) = {1};  //Ground


Physical Surface(106)={6,9,7,8,11};
//+

// Domain volume
Physical Volume(202) = {99};
