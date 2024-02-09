boundary_layer=5e-3;

size_boundary=0.1;
size_bas=boundary_layer/8;
size_bottom=boundary_layer/40;
height=0.3;
nb_square=5;
Amplitude=0.05;
discretisation_pts=80;
//+
Point(1) = {1, 0, 0, size_bottom};
//+
Point(2) = {1, height, 0, size_boundary};
//+
Point(3) = {0, height, 0, size_boundary};
//+
Point(4) = {0, 0, 0, size_bottom};
Point(5) = {0,boundary_layer,0,size_bas};
Point(6) = {1,boundary_layer,0,size_bas};

k=5000;
For i In {1:nb_square}

    Point(k)   = {(1/4)/nb_square+(i-1)/nb_square+boundary_layer,boundary_layer,0,size_bas};
    If(i>1)
        Line(k-1)={k-1,k};
    EndIf
    Point(k+1) = {(1/4)/nb_square+(i-1)/nb_square+boundary_layer,boundary_layer-Amplitude,0,size_bas};
    Point(k+2) = {(3/4)/nb_square+(i-1)/nb_square-boundary_layer,boundary_layer-Amplitude,0,size_bas};
    Point(k+3) = {(3/4)/nb_square+(i-1)/nb_square-boundary_layer,boundary_layer,0,size_bas};
    Line(k)={k,k+1};
    Line(k+1)={k+1,k+2};
    Line(k+2)={k+2,k+3};
    k=k+4;
EndFor
Line(k-1)={k-1,6};
Line(4999)={5,5000};
max=k-1;


//+
Line(1) = {1, 6};
//+
Line(2) = {6, 2};
//+
Line(3) = {2, 3};
Line(4)={3,5};
Line(5)={5,4};
//+
k=7;


linep=7;
linec=1000;
list={5,6};

For i In {1:nb_square}
    Point(k) = {(1/4)/nb_square+(i-1)/nb_square,0,0,size_bottom};
    If(i>1)
        Line(linec)={k-1,k};
        list={list[],linec};
    EndIf
    Point(k+1) = {(1/4)/nb_square+(i-1)/nb_square,-Amplitude,0,size_bottom};
    Point(k+2) = {(3/4)/nb_square+(i-1)/nb_square,-Amplitude,0,size_bottom};
    Point(k+3) = {(3/4)/nb_square+(i-1)/nb_square,0,0,size_bottom};
    Line(linec+1)={k,k+1};
    Line(linep+1)={k+1,k+2};
    Line(linec+2)={k+2,k+3};
    list={list[],linec+1,linep+1,linec+2};
    k=k+4;
    linep=linep+1;
    linec=linec+3;
EndFor
Line(6)={4,7};
Line(k-1)={k-1,1};
//+
Physical Curve("Cathode",11)={1001:linec-1};
Physical Curve("Pore",10)={8:linep};
Physical Curve("VerEdgeInlet", 12) = {4:5};
Physical Curve("VerEdgeOutlet", 14) = {1,2};
//+
Physical Curve("HorEdgeTop",13)={3};

Curve Loop(1)={2,3,4,4999:max};
//Curve Loop(2)={1,-max:-4999,5,6:k-1};
NumPoints = #list[];

Printf("The Array Contents are") ;
For index In {0:NumPoints-1}
    Printf("%g ",list[index]) ;
EndFor
Curve Loop(2)={1,-max:-4999,list[],k-1};
Plane Surface(1)={1};
Plane Surface(2)={2};
Physical Surface(14)={1,2};

Mesh 2;

