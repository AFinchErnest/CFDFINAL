SetFactory("OpenCASCADE");
// ----------------------------
//Constants
// ------------------------
a = 40; b = 100; c =25; //Building Size params
delta_y = 10; delta_x = 10; //distnace between buildings
L = delta_y+2*b+80; W = delta_x+2*a+50; W_offset = a+20; H = c+15; //Domain Size Params

h1 = a/10;
h2 = 0.5*h1;

//--------------------------
// CITY DOMAIN DEFINITION
//--------------------------

// Outer Domain Points 
Point(1) = {0, 0, 0, h1};
Point(2) = {W, 0, 0, h1};
Point(3) = {W, L, 0, h1};
Point(4) = {0, L, 0, h1};
Point(5) = {0, 0, H, h1};
Point(6) = {W, 0, H, h1};
Point(7) = {W, L, H, h1};
Point(8) = {0, L, H, h1};

// Outer domain Lines
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

// outer city surfaces
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

//--------------------------
// BUILDINGS DEFFINITIONS
//--------------------------
//Building 1
x0 = (W_offset - a)/2; x1 = (W_offset + a)/2;
y0 = (L - delta_y)/2 - b; y1 = (L - delta_y)/2;

Point(101) = {x0, y0, 0, h2};
Point(102) = {x1, y0, 0, h2};
Point(103) = {x1, y1, 0, h2};
Point(104) = {x0, y1, 0, h2};
Point(105) = {x0, y0, c, h2};
Point(106) = {x1, y0, c, h2};
Point(107) = {x1, y1, c, h2};
Point(108) = {x0, y1, c, h2};

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

// Surfaces
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

//Building 2
x0 = (W_offset - a)/2; x1 = (W_offset + a)/2;
y0 = (L + delta_y)/2 + b; y1 = (L + delta_y)/2;

Point(109) = {x0, y0, 0, h2};
Point(110) = {x1, y0, 0, h2};
Point(111) = {x1, y1, 0, h2};
Point(112) = {x0, y1, 0, h2};
Point(113) = {x0, y0, c, h2};
Point(114) = {x1, y0, c, h2};
Point(115) = {x1, y1, c, h2};
Point(116) = {x0, y1, c, h2};

Line(149) = {109, 112};
Line(150) = {112, 111};
Line(151) = {111, 110};
Line(152) = {110, 109};

Line(153) = {115, 116};
Line(154) = {116, 113};
Line(155) = {113, 114};
Line(156) = {114, 115};

Line(157) = {115, 111};
Line(158) = {116, 112};
Line(159) = {113, 109};
Line(160) = {110, 114};

// Surfaces
Line Loop(207) = {149,150,151,152};
Plane Surface(307) = {207}; //bottom
Line Loop(208) = {153,154,155,156};
Plane Surface(308) = {208}; //Top
Line Loop(209) = {149,159,154,158};
Plane Surface(309) = {209}; //South
Line Loop(213) = {157,150,158,153};
Plane Surface(310) = {213}; //East
Line Loop(211) = {151,160,156,157};
Plane Surface(311) = {211}; //North
Line Loop(212) = {155,159,152,160};
Plane Surface(312) = {212}; //West

Surface Loop(315) = {307,308,309,310,311,312};
Volume(315) = {315};

//Building 3
x0 = (W_offset - a)/2 + a + delta_x; x1 = (W_offset - a)/2 +2*a + delta_x;
y0 = (L + delta_y)/2 + b; y1 = (L + delta_y)/2;

Point(117) = {x0, y0, 0, h2};
Point(118) = {x1, y0, 0, h2};
Point(119) = {x1, y1, 0, h2};
Point(120) = {x0, y1, 0, h2};
Point(121) = {x0, y0, c, h2};
Point(122) = {x1, y0, c, h2};
Point(123) = {x1, y1, c, h2};
Point(124) = {x0, y1, c, h2};

Line(37) = {122, 123};
Line(38) = {123, 124};
Line(39) = {124, 121};
Line(40) = {121, 122};

Line(41) = {118, 119};
Line(42) = {119, 120};
Line(43) = {120, 117};
Line(44) = {117, 118};

Line(45) = {118, 122};
Line(46) = {117, 121};
Line(47) = {119, 123};
Line(48) = {120, 124};

// Surfaces
Line Loop(216) = {41,42,43,44};
Plane Surface(316) = {216}; //bottom
Line Loop(217) = {37,38,39,40};
Plane Surface(317) = {217}; //Top
Line Loop(218) = {39,48,43,46};
Plane Surface(318) = {218}; //South
Line Loop(228) = {40,46,44,45};
Plane Surface(328) = {228}; //East
Line Loop(240) = {37,45,41,47};
Plane Surface(340) = {240}; //North
Line Loop(221) = {47,42,48,38};
Plane Surface(321) = {221}; //West

Surface Loop(316) = {316,317,318,328,340,321};
Volume(316) = {316};



//Building 4
x0 = (W_offset - a)/2 + a + delta_x; x1 = (W_offset - a)/2 +2*a + delta_x;
y0 = (L - delta_y)/2 - b; y1 = (L - delta_y)/2;

Point(125) = {x0, y0, 0, h2};
Point(126) = {x1, y0, 0, h2};
Point(127) = {x1, y1, 0, h2};
Point(128) = {x0, y1, 0, h2};
Point(129) = {x0, y0, c, h2};
Point(130) = {x1, y0, c, h2};
Point(131) = {x1, y1, c, h2};
Point(132) = {x0, y1, c, h2};

Line(85) = {131, 130};
Line(86) = {130, 129};
Line(87) = {129, 132};
Line(88) = {132, 131};

Line(89) = {127, 126};
Line(90) = {126, 125};
Line(91) = {125, 128};
Line(92) = {128, 127};

Line(93) = {130, 126};
Line(94) = {125, 129};
Line(95) = {131, 127};
Line(96) = {132, 128};

// Surfaces

Line Loop(230) = {89,90,91,92};
Plane Surface(330) = {230}; //bottom
Line Loop(223) = {85,86,87,88};
Plane Surface(323) = {223}; //Top
Line Loop(224) = {96,91,94,87};
Plane Surface(324) = {224}; //South
Line Loop(225) = {93,86,94,90};
Plane Surface(325) = {225}; //East
Line Loop(231) = {93,89,95,85};
Plane Surface(331) = {231}; //North
Line Loop(227) = {95,92,96,88};
Plane Surface(327) = {227}; //West

Surface Loop(317) = {330,323,324,325,331,327};
Volume(317) = {317};


// ----------------------------
// subtract volumes
// ----------------------------
BooleanDifference(99)={ Volume{41}; Delete; } { Volume{311}; Delete; };
BooleanDifference(100)={Volume{99}; Delete; } {Volume{315}; Delete;};
BooleanDifference(101)={Volume{100}; Delete; } {Volume{316}; Delete;};
BooleanDifference(102)={Volume{101}; Delete; } {Volume{317}; Delete;};

// ----------------------------
// Mesh control
// ----------------------------
Mesh.CharacteristicLengthMax = h1;

// ----------------------------
// Physical Groups 
// ----------------------------

// Domain faces
Physical Surface(100) = {2}; //WEST
Physical Surface(102) = {3} ;//South
Physical Surface(101) = {4} ;//North
Physical Surface(103) = {5} ;//EAST
Physical Surface(104) = {22};  //SKY
Physical Surface(105) = {1};  //Ground
Physical Surface(106)={7,6,8,23,9}; //Building1
Physical Surface(107) = {25,17,16,14,15}; //Building2
Physical Surface(108) = {24,11,13,10,12}; //Building3
Physical Surface(109) = {21,26,20,18,19}; //Building4

// Domain volume
Physical Volume(202) = {102};



