SetFactory("OpenCASCADE");
// -------------------
//CONSTANTS
// -------------------
H = 20; Lx = 226; Ly = 137; // Surrounding Dimensions
c1 = 12; c2 = 15; c3 = 12; c4 = 15; c5 = 6; c6 = 12; //Building Heights
h= 5; h1 = 0.5*h; //Mesh Size


// -------------------
//Surroundings
// -------------------
Point(1) = {0,0,0,h1};
Point(2) = {Lx,0,0,h1};
Point(3) = {0,Ly,0,h1};
Point(4) = {Lx,Ly,0,h1};
Point(37) = {0,0,H,h1};
Point(38) = {Lx,0,H,h1};
Point(39) = {0,Ly,H,h1};
Point(40) = {Lx,Ly,H,h1};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Line(5) = {37, 38};
Line(6) = {38, 40};
Line(7) = {40, 39};
Line(8) = {39, 37};

Line(9) = {1, 37};
Line(10) = {2, 38};
Line(11) = {4, 40};
Line(12) = {3, 39};

Curve Loop(1) = {8, 5, 6, 7}; // SKY
Plane Surface(1) = {1};
Curve Loop(2) = {1, 2, 3, 4}; //GROUND
Plane Surface(2) = {2};
Curve Loop(3) = {9, 5, -10, -1};
Plane Surface(3) = {3};
Curve Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};
Curve Loop(5) = {3, 12, -7, -11};
Plane Surface(5) = {5};
Curve Loop(6) = {4, 9, -8, -12};
Plane Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

// -------------------
//Building 1
// -------------------
Point(5) = {10,10,0,h1};
Point(6) = {10,55,0,h1};
Point(9) = {40,30,0,h1};
Point(10) = {40,55,0,h1};
Point(11) = {50,10,0,h1};
Point(12) = {50,30,0,h1};

Point(32) = {10,10,c1,h1};
Point(33) = {10,55,c1,h1};
Point(34) = {40,55,c1,h1};
Point(35) = {50,10,c1,h1};
Point(36) = {50,30,c1,h1};
Point(41) = {40,30,c1,h1};

Line(13) = {11, 5};
Line(14) = {5, 6};
Line(15) = {6, 10};
Line(16) = {10, 9};
Line(17) = {9, 12};
Line(18) = {12, 11};

Line(19) = {35, 32};
Line(20) = {32, 33};
Line(21) = {33, 34};
Line(22) = {34, 41};
Line(23) = {41, 36};
Line(24) = {36, 35};

Line(25) = {5, 32};
Line(26) = {6, 33};
Line(27) = {10, 34};
Line(28) = {9, 41};
Line(29) = {36, 12};
Line(30) = {35, 11};

Curve Loop(7) = {13, 25, -19, 30};
Plane Surface(7) = {7};
Curve Loop(8) = {25, 20, -26, -14};
Plane Surface(8) = {8};
Curve Loop(9) = {15, 27, -21, -26};
Plane Surface(9) = {9};
Curve Loop(10) = {16, 28, -22, -27};
Plane Surface(10) = {10};
Curve Loop(11) = {23, 29, -17, 28};
Plane Surface(11) = {11};
Curve Loop(12) = {29, 18, -30, -24};
Plane Surface(12) = {12};
Curve Loop(13) = {20, 21, 22, 23, 24, 19};
Plane Surface(13) = {13};
Curve Loop(14) = {15, 16, 17, 18, 13, 14};
Plane Surface(14) = {14};

Surface Loop(2) = {7,8,9,10,11,12,13,14};
Volume(2) = {2};

// -------------------
//Building 2
// -------------------
Point(7) = {10,62,0,h1};
Point(8) = {10,107,0,h1};
Point(17) = {101,62,0,h1};
Point(18) = {101,107,0,h1};

Point(42) = {10,62,c2,h1};
Point(43) = {10,107,c2,h1};
Point(44) = {101,62,c2,h1};
Point(45) = {101,107,c2,h1};

Line(31) = {42, 43};
Line(32) = {43, 45};
Line(33) = {45, 44};
Line(34) = {44, 42};
Line(35) = {7, 17};
Line(36) = {17, 18};
Line(37) = {18, 8};
Line(38) = {8, 7};
Line(39) = {17, 44};
Line(40) = {18, 45};
Line(41) = {8, 43};
Line(42) = {7, 42};

Curve Loop(15) = {31, 32, 33, 34};
Plane Surface(15) = {15};
Curve Loop(16) = {35, 36, 37, 38};
Plane Surface(16) = {16};
Curve Loop(17) = {42, -34, -39, -35};
Plane Surface(17) = {17};
Curve Loop(18) = {38, 42, 31, -41};
Plane Surface(18) = {18};
Curve Loop(19) = {41, 32, -40, 37};
Plane Surface(19) = {19};
Curve Loop(20) = {36, 40, 33, -39};
Plane Surface(20) = {20};

Surface Loop(3) = {15,16,17,18,19,20};
Volume(3) = {3};

// -------------------
//Building 3
// -------------------
Point(13) = {56,10,0,h1};
Point(14) = {56,55,0,h1};
Point(15) = {101,10,0,h1};
Point(16) = {101,55,0,h1};

Point(46) = {56,10,c3,h1};
Point(47) = {56,55,c3,h1};
Point(48) = {101,10,c3,h1};
Point(49) = {101,55,c3,h1};

Line(43) = {47, 49};
Line(44) = {49, 48};
Line(45) = {48, 46};
Line(46) = {46, 47};
Line(47) = {14, 16};
Line(48) = {16, 15};
Line(49) = {15, 13};
Line(50) = {13, 14};
Line(51) = {14, 47};
Line(52) = {16, 49};
Line(53) = {15, 48};
Line(54) = {13, 46};

Curve Loop(21) = {45, 46, 43, 44};
Plane Surface(21) = {21};
Curve Loop(22) = {48, 49, 50, 47};
Plane Surface(22) = {22};
Curve Loop(23) = {53, 45, -54, -49};
Plane Surface(23) = {23};
Curve Loop(24) = {48, 53, -44, -52};
Plane Surface(24) = {24};
Curve Loop(25) = {47, 52, -43, -51};
Plane Surface(25) = {25};
Curve Loop(26) = {50, 51, -46, -54};
Plane Surface(26) = {26};

Surface Loop(4) = {21,22,23,24,25,26};
Volume(4) = {4};

// -------------------
//Building 4
// -------------------
Point(21) = {126,69,0,h1};
Point(22) = {126,107,0,h1};
Point(30) = {216,69,0,h1};
Point(31) = {216,107,0,h1};
Point(50) = {126,69,c4,h1};
Point(51) = {126,107,c4,h1};
Point(52) = {216,69,c4,h1};
Point(53) = {216,107,c4,h1};

Line(55) = {21, 50};
Line(56) = {22, 51};
Line(57) = {31, 53};
Line(58) = {30, 52};
Line(59) = {52, 50};
Line(60) = {50, 51};
Line(61) = {51, 53};
Line(62) = {53, 52};
Line(63) = {30, 31};
Line(64) = {31, 22};
Line(65) = {22, 21};
Line(66) = {21, 30};

Curve Loop(27) = {66, 58, 59, -55};
Plane Surface(27) = {27};
Curve Loop(28) = {55, 60, -56, 65};
Plane Surface(28) = {28};
Curve Loop(29) = {64, 56, 61, -57};
Plane Surface(29) = {29};
Curve Loop(30) = {62, -58, 63, 57};
Plane Surface(30) = {30};
Curve Loop(31) = {59, 60, 61, 62};
Plane Surface(31) = {31};
Curve Loop(32) = {66, 63, 64, 65};
Plane Surface(32) = {32};

Surface Loop(5) = {27,28,29,30,31,32};
Volume(5) = {5};

// -------------------
//Buldings 5 and 6
// -------------------

//NOTE: Since building 6 and 5 share a face 
//      they need to be one single Volume.

//Building 6 Points
Point(26) = {156,10,0,h1};
Point(27) = {156,50,0,h1};
Point(28) = {216,10,0,h1};
Point(29) = {216,50,0,h1};
Point(60) = {156,10,c6,h1};
Point(61) = {156,50,c6,h1};
Point(62) = {216,10,c6,h1};
Point(63) = {216,50,c6,h1};

//Building 5 Points
Point(19) = {126,10,0,h1};
Point(20) = {126,28,0,h1};
Point(23) = {146,28,0,h1};
Point(24) = {146,43,0,h1};
Point(25) = {156,43,0,h1};
Point(54) = {126,10,c5,h1};
Point(55) = {126,28,c5,h1};
Point(56) = {146,28,c5,h1};
Point(57) = {146,43,c5,h1};
Point(58) = {156,10,c5,h1};
Point(59) = {156,43,c5,h1};

//Buildings 5&6 Lines
Line(67) = {19, 26};
Line(69) = {25, 24};
Line(70) = {24, 23};
Line(71) = {23, 20};
Line(72) = {20, 19};
Line(73) = {19, 54};
Line(74) = {20, 55};
Line(75) = {23, 56};
Line(76) = {24, 57};
Line(77) = {25, 59};
Line(79) = {55, 56};
Line(80) = {56, 57};
Line(81) = {57, 59};
Line(82) = {58, 59};
Line(83) = {55, 54};
Line(84) = {54, 58};
Line(85) = {26, 28};
Line(86) = {28, 29};
Line(87) = {29, 27};
Line(88) = {27, 25};
Line(89) = {27, 61};
Line(90) = {61, 60};
Line(91) = {60, 58};
Line(92) = {61, 63};
Line(93) = {63, 62};
Line(94) = {62, 60};
Line(95) = {29, 63};
Line(96) = {28, 62};

Curve Loop(33) = {85, 96, 94, 91, -84, -73, 67};
Plane Surface(33) = {33};
Curve Loop(34) = {96, -93, -95, -86};
Plane Surface(34) = {34};
Curve Loop(35) = {87, 89, 92, -95};
Plane Surface(35) = {35};
Curve Loop(36) = {88, 77, -82, -91, -90, -89};
Plane Surface(36) = {36};
Curve Loop(37) = {77, -81, -76, -69};
Plane Surface(37) = {37};
Curve Loop(38) = {76, -80, -75, -70};
Plane Surface(38) = {38};
Curve Loop(39) = {75, -79, -74, -71};
Plane Surface(39) = {39};
Curve Loop(40) = {72, 73, -83, -74};
Plane Surface(40) = {40};
Curve Loop(41) = {81, -82, -84, -83, 79, 80};
Plane Surface(41) = {41};
Curve Loop(42) = {92, 93, 94, -90};
Plane Surface(42) = {42};
Curve Loop(43) = {67, 85, 86, 87, 88, 69, 70, 71, 72};
Plane Surface(43) = {43};

Surface Loop(6) = {33,34,35,36,37,38,39,40,41,42,43};
Volume(6) = {6};

// -------------------
// SUBTRACTING VOLUMES
// -------------------
BooleanDifference(7)={ Volume{1}; Delete; } { Volume{2}; Delete; };
BooleanDifference(8)={ Volume{7}; Delete; } { Volume{3}; Delete; };
BooleanDifference(9)={ Volume{8}; Delete; } { Volume{4}; Delete; };
BooleanDifference(10)={ Volume{9}; Delete; } { Volume{5}; Delete; };
BooleanDifference(11)={ Volume{10}; Delete; } { Volume{6}; Delete; };

// -------------------
// PHYSICAL IDS
// -------------------
// Domain faces
Physical Surface(100) = {2}; //WEST
Physical Surface(102) = {5} ;//South
Physical Surface(101) = {3} ;//North
Physical Surface(103) = {4} ;//EAST
Physical Surface(104) = {1};  //SKY
Physical Surface(105) = {6};  //Ground
Physical Surface(106)={7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38}; //Building1

// Domain volume
Physical Volume(202) = {11};

// Mesh control
Mesh.CharacteristicLengthMax = h1;
