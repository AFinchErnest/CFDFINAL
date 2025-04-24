SetFactory("OpenCASCADE");

L = 20; W = 10; H = 10;
a = 2; b = 6;

h1 = 1.0;
h2 = 0.5*h1;

// City Domain
Box(1) = {0, 0, 0, W, L, H};

// Build 1
Box(2) = {(W - a)/2, (L - a)/2, 0, a, a, b};

// Subtract building from domain
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// Assign sim domain
Physical Volume(1) = {3};

// Extract and tag sides of domain
north[] = Surface In BoundingBox(0, 0, 0, 0, L, H);
Physical Surface(1) = {north[]};

south[] = Surface In BoundingBox(W, 0, 0, W, L, H);
Physical Surface(2) = {south[]};

west[] = Surface In BoundingBox(0, 0, 0, W, 0, H);
Physical Surface(3) = {west[]};
east[] += Surface In BoundingBox(0, L, 0, W, L, H);
Physical Surface(4) = {east[]};
walls[] += Surface In BoundingBox(0, 0, H, W, L, H);
Physical Surface(5) = {sky[]};

ground[] = Surface In BoundingBox(0, 0, 0, W, L, 0);
Physical Surface(6) = {ground[]};

building[] = Surface In BoundingBox((W - a)/2, (L - a)/2, 0, (W + a)/2, (L + a)/2, b);
Physical Surface(7) = {building[]};

// Mesh control
Mesh.CharacteristicLengthMax = h1;
