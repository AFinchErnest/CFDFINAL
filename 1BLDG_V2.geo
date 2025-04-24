SetFactory("OpenCASCADE");

L = 20; W = 10; H = 10;
a = 2; b = 6;

h1 = 1.0;
h2 = 0.5*h1;

// Outer block
Box(1) = {0, 0, 0, W, L, H};

// Building
Box(2) = {(W - a)/2, (L - a)/2, 0, a, a, b};

// Subtract building from domain
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// Assign volume
Physical Volume("AirDomain") = {3};

// Extract and tag surfaces
inlet[] = Surface In BoundingBox(0, 0, 0, 0, L, H);
Physical Surface("Inlet") = {inlet[]};

outlet[] = Surface In BoundingBox(W, 0, 0, W, L, H);
Physical Surface("Outlet") = {outlet[]};

walls[] = Surface In BoundingBox(0, 0, 0, W, 0, H);
walls[] += Surface In BoundingBox(0, L, 0, W, L, H);
walls[] += Surface In BoundingBox(0, 0, H, W, L, H);
Physical Surface("Walls") = {walls[]};

ground[] = Surface In BoundingBox(0, 0, 0, W, L, 0);
Physical Surface("Ground") = {ground[]};

building[] = Surface In BoundingBox((W - a)/2, (L - a)/2, 0, (W + a)/2, (L + a)/2, b);
Physical Surface("Building") = {building[]};

// Mesh control
Mesh.CharacteristicLengthMax = h1;
