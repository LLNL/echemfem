La = 5.;
L = 2.;
Lb = 5.;
h = 1.;
lc = 0.25*h;
lc1 = 0.08*h;
Point(1) = {0, 0, 0, lc};
Point(2) = {La, 0, 0, lc1};
Point(3) = {La+L, 0, 0, lc1};
Point(4) = {La+L+Lb, 0, 0, lc};
Point(5) = {La+L+Lb, h, 0, lc};
Point(6) = {La+L, h, 0, lc1};
Point(7) = {La, h, 0, lc1};
Point(8) = {0, h, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Curve Loop(9) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {9};

Physical Curve("Inlet", 10) = {8};
Physical Curve("Outlet", 11) = {4};
Physical Curve("Anode", 12) = {6};
Physical Curve("Cathode", 13) = {2};
Physical Curve("Wall", 14) = {7, 1, 5, 3};
Physical Surface("Channel", 2) = {1};


Recombine Surface{1};
