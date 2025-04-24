h=.00;
\\+
Point(newp) = {1.0000,0.0013,0.0,h};
Point(newp) = {0.9500,0.0114,0.0,h};
Point(newp) = {0.9000,0.0208,0.0,h};
Point(newp) = {0.8000,0.0375,0.0,h};
Point(newp) = {0.7000,0.0518,0.0,h};
Point(newp) = {0.6000,0.0636,0.0,h};
Point(newp) = {0.5000,0.0724,0.0,h};
Point(newp) = {0.4000,0.0780,0.0,h};
Point(newp) = {0.3000,0.0788,0.0,h};
Point(newp) = {0.2500,0.0767,0.0,h};
Point(newp) = {0.2000,0.0726,0.0,h};
Point(newp) = {0.1500,0.0661,0.0,h};
Point(newp) = {0.1000,0.0563,0.0,h};
Point(newp) = {0.0750,0.0496,0.0,h};
Point(newp) = {0.0500,0.0413,0.0,h};
Point(newp) = {0.0250,0.0299,0.0,h};
Point(newp) = {0.0125,0.0215,0.0,h};
Point(newp) = {0.0000,0.0000,0.0,h};
Point(newp) = {0.0125,-0.0165,0.0,h};
Point(newp) = {0.0250,-0.0227,0.0,h};
Point(newp) = {0.0500,-0.0301,0.0,h};
Point(newp) = {0.0750,-0.0346,0.0,h};
Point(newp) = {0.1000,-0.0375,0.0,h};
Point(newp) = {0.1500,-0.0410,0.0,h};
Point(newp) = {0.2000,-0.0423,0.0,h};
Point(newp) = {0.2500,-0.0422,0.0,h};
Point(newp) = {0.3000,-0.0412,0.0,h};
Point(newp) = {0.4000,-0.0380,0.0,h};
Point(newp) = {0.5000,-0.0334,0.0,h};
Point(newp) = {0.6000,-0.0276,0.0,h};
Point(newp) = {0.7000,-0.0214,0.0,h};
Point(newp) = {0.8000,-0.0150,0.0,h};
Point(newp) = {0.9000,-0.0082,0.0,h};
Point(newp) = {0.9500,-0.0048,0.0,h};
Point(newp) = {1.0000,-0.0013,0.0,h};
Point(newp) = {1.5000,-0.1028,0.0,h};
Point(newp) = {1.5000,.1393,0.0,h};
Point(newp) = {-.50000,-0.1028,0.0,h};
Point(newp) = {-.50000,.1393,0.0,h};//+
Spline(1) = {34, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
//+
Line(2) = {38, 36};
//+
Line(3) = {36, 35};
//+
Line(4) = {35, 37};
//+
Line(5) = {38, 37};
//+
Curve Loop(1) = {2, 3, 4, -5};
//+
Curve Loop(2) = {1};
//+
Surface(1) = {1, 2};

//+
Field[2] = Distance;
//+
Field[2].CurvesList = {1};
//+
Field[2].NumPointsPerCurve = 100;
//+
Field[3] = Threshold;
//+
Field[3].InField = 1;
//+
Field[3].DistMax = 1;
//+
Field[3].DistMin = 0;
//+
Field[3].SizeMin = 0.01;
//+
Field[3].SizeMax = 0.1;
//+
Field[3].InField = 2;
//+
Background Field = 3;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
//+


Field[1] = Distance;
//+
Field[1].PointsList = {14, 15, 16, 17, 18};
//+
Field[1].PointsList = {14, 15, 16, 17, 18, 33, 34, 1};
//+
Field[4] = Threshold;
//+
Field[4].DistMax = 0.2;
//+
Field[4].DistMin = 0;
//+
Field[4].SizeMax = 0.1;
//+
Field[4].SizeMin = 0.005;
//+
Field[4].InField = 1;
//+
Field[4].DistMax = 0.05;
//+
Field[5] = Min;
//+
Field[5].FieldsList = {4, 3};
//+
Background Field = 5;
//+
Field[3].SizeMin = 0.05;
//+
Field[3].SizeMax = 0.5;
//+
Field[4].SizeMax = 0.05;
//+
Field[4].DistMax = 0.1;
//+
Field[3].SizeMax = 0.1;
//+
Field[4].StopAtDistMax = 1;
//+
Field[3].DistMax = 0.55;
//+
Field[4].SizeMax = 0.01;
//+
Field[3].SizeMax = 0.05;
//+
Field[3].SizeMin = 0.01;

