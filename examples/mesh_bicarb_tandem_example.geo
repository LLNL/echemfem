La = 0.00025;
L = 0.001;
Lb = 0.00025;
h = 6e-4;
lc = L/7.; 
delta = 1e-8; //boundary layer 1
delta2 = 1e-6;// refined region around points
delta2b = 2e-4; //start of electrode
delta3 = 0.4*h; //boundary layer 2
delta3b = 0.01*h; //transition between layers
lc2 = lc;

Point(1) = {0, 0, 0, lc};
Point(2) = {La-delta2, 0, 0, lc2};
Point(3) = {La, 0, 0, lc2};
Point(4) = {La+0.1*delta2b, 0, 0, lc2};
Point(5) = {La+delta2b, 0, 0, lc2};
Point(6) = {La+L-delta2, 0, 0, lc2};
Point(7) = {La+L, 0, 0, lc2};
Point(8) = {La+L+delta2, 0, 0, lc2};
Point(9) = {La+L+Lb, 0, 0, lc};

Point(10) = {0, delta, 0, lc2};
Point(11) = {La-delta2, delta, 0, lc2};
Point(12) = {La, delta, 0, lc2};
Point(13) = {La+0.1*delta2b, delta, 0, lc2};
Point(14) = {La+delta2b, delta, 0, lc2};
Point(15) = {La+L-delta2, delta, 0, lc2};
Point(16) = {La+L, delta, 0, lc2};
Point(17) = {La+L+delta2, delta, 0, lc2};
Point(18) = {La+L+Lb, delta, 0, lc2};

Point(19) = {0, h, 0, lc};
Point(20) = {La+L+Lb, h, 0, lc};

Point(21) = {La + delta2b, delta3b, 0, lc2};
Point(22) = {La + L - delta2, delta3b, 0, lc2};

Point(23) = {La + delta2b, delta3, 0, lc2};
Point(24) = {La + L - delta2, delta3, 0, lc2};

//horizontal lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};

Line(9) = {10,11};
Line(10) = {11,12};
Line(11) = {12,13};
Line(12) = {13,14};
Line(13) = {14,15};
Line(14) = {15,16};
Line(15) = {16,17};
Line(16) = {17,18};

Line(17) = {19,20};

//vertical lines
Line(18) = {1,10};
Line(19) = {2,11};
Line(20) = {3,12};
Line(21) = {4,13};
Line(22) = {5,14};
Line(23) = {6,15};
Line(24) = {7,16};
Line(25) = {8,17};
Line(26) = {9,18};

Line(27) = {10,19};
Line(28) = {18,20};

// extra lines
Line(29) = {21,22};
Line(30) = {14,21};
Line(31) = {15,22};

Line(32) = {23,24};
Line(33) = {21,23};
Line(34) = {22,24};

Curve Loop(26) = {1,19,-9,-18};
Curve Loop(27) = {2,20,-10,-19};
Curve Loop(28) = {3,21,-11,-20};
Curve Loop(29) = {4,22,-12,-21};
Curve Loop(30) = {5,23,-13,-22};
Curve Loop(31) = {6,24,-14,-23};
Curve Loop(32) = {7,25,-15,-24};
Curve Loop(33) = {8,26,-16,-25};

Curve Loop(34) = {9,10,11,12,30,33,32,-34,-31,14,15,16,28,-17,-27};
Curve Loop(35) = {13,31,-29,-30};
Curve Loop(36) = {29,34,-32,-33};

Plane Surface(1) = {26};
Plane Surface(2) = {27};
Plane Surface(3) = {28};
Plane Surface(4) = {29};
Plane Surface(5) = {30};
Plane Surface(6) = {31};
Plane Surface(7) = {32};
Plane Surface(8) = {33};
Plane Surface(11) = {34};
Plane Surface(9) = {35};
Plane Surface(10) = {36};

Physical Curve("Inlet", 1) = {18,27};
Physical Curve("Outlet", 2) = {26,28};
Physical Curve("Bulk", 3) = {17};
Physical Curve("Cathode", 4) = {3,4,5,6};
Physical Curve("Wall", 5) = {1,2,7,8};
Physical Surface("Channel", 1) = {1,2,3,4,5,6,7,8,9,10,11};


//y-dir
Transfinite Curve{18,19,20,21,22,23,24,25,26} = 26 Using Progression 1.15;


//x-dir electrode
Transfinite Curve{2,10} = 3 Using Progression 1/1.4;
Transfinite Curve{3,11} = 5 Using Progression 1.35;
Transfinite Curve{4,12} = 3 Using Progression 1.0;
Transfinite Curve{6,14} = 3 Using Progression 1/1.4;
Transfinite Curve{7,15} = 3 Using Progression 1.4;
Transfinite Curve{5,13,29,32} = 59 Using Bump 0.4;

//x-dir inlet transition
Transfinite Curve{1,9} = 51 Using Progression .9;

//x-dir outlet transition
Transfinite Curve{8,16} = 71 Using Progression 1./.95;

//y-dir diffusion bdlayer
Transfinite Curve{30,31} = 7 Using Progression 1.4;
Transfinite Curve{33,34} = 21 Using Progression 1;

//boundary layer
Transfinite Surface{1} = {1,2,11,10};
Transfinite Surface{2} = {2,3,12,11};
Transfinite Surface{3} = {3,4,13,12};
Transfinite Surface{4} = {4,5,14,13};
Transfinite Surface{5} = {5,6,15,14};
Transfinite Surface{6} = {6,7,16,15};
Transfinite Surface{7} = {7,8,17,16};
Transfinite Surface{8} = {8,9,18,17};
Transfinite Surface{9} = {14,15,22,21};
Transfinite Surface{10} = {21,22,24,23};

//coarse quads

// Electrode start and end
Field[1] = Distance;
// Field[1].PointsList = {10,11,12,13,14,15};
Field[1].PointsList = {12,16};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 8e-11;
Field[2].SizeMax = lc;
Field[2].DistMin = 2e-10;
Field[2].DistMax = 0.8*h;

Field[3] = Distance;
Field[3].CurvesList = {13,14,15,16};
Field[3].NumPointsPerCurve = 100;

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 0.01/100.;
Field[4].SizeMax = lc;
Field[4].DistMin = 0.5*h;
Field[4].DistMax = h;
//Field[4].Sigmoid = 1;

Field[5] = Distance;
Field[5].CurvesList = {16};
Field[5].NumPointsPerCurve = 100;

Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = 0.35*lc;
Field[6].SizeMax = lc;
Field[6].DistMin = delta3;
Field[6].DistMax = 0.5*h;

Field[7] = Distance;
Field[7].CurvesList = {10};
Field[7].NumPointsPerCurve = 100;

Field[8] = Threshold;
Field[8].InField = 7;
Field[8].SizeMin = 0.8*delta2b/22;
Field[8].SizeMax = lc;
Field[8].DistMin = delta3*0.1;
Field[8].DistMax = 0.5*h;

Field[9] = Min;
Field[9].FieldsList = {2,4,6, 8}; 
Background Field = 9;

Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
Recombine Surface{1,2,3,4,5,6,7,8,9,10,11};
