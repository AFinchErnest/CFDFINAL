
SetFactory("OpenCASCADE");
// Box dimensions
Lx = 1.0; // length in x
Ly = 1.0; // length in y
Lz = 0.5; // height in z

// Mesh resolution
res = 0.1;

// Define box
Box(1) = {0, 0, 0, Lx, Ly, Lz};

// Physical volume
Physical Volume(1) = {1};

// Define surfaces
// Face naming convention: -Z = bottom, +Z = top, +X = right, -X = left, +Y = back, -Y = front

Physical Surface(1) = {1};  // z = 0
Physical Surface(2)    = {6};  // z = Lz
Physical Surface(3)   = {3};  // x = 0
Physical Surface(4)  = {5};  // x = Lx
Physical Surface(5)  = {2};  // y = 0
Physical Surface(6)   = {4};  // y = Ly

Mesh.CharacteristicLengthMax = res;