w = 0.5;
l = 12*w;
t = 0.6*w;
d  = 0.25*w;

s1 = w;
s2 = 3.5*w;
s  = 0.25*w;

Point(1) = {0, 0, 0, d};
Point(2) = {s1 - s, 0, 0, d};
Point(3) = {s1, 0, 0, d};
Point(4) = {s2 - s, 0, 0, d};
Point(5) = {s2, 0, 0, d};
Point(6) = {0.5*(l-t), 0, 0, d};
Point(7) = {0.5*(l+t), 0, 0, d};
Point(8) = {l - s2, 0, 0, d};
Point(9) = {l - s2 + s, 0, 0, d};
Point(10) = {l - s1, 0, 0, d};
Point(11) = {l - s1 + s, 0, 0, d};
Point(12) = {l, 0, 0, d};
Point(13) = {l, w, 0, d};
Point(14) = {l - s1 + s, w, 0, d};
Point(15) = {l - s1, w, 0, d};
Point(16) = {l - s2 + s, w, 0, d};
Point(17) = {l - s2, w, 0, d};
Point(18) = {0.5*(l+t), w, 0, d};
Point(19) = {0.5*(l-t), w, 0, d};
Point(20) = {s2, w, 0, d};
Point(21) = {s2 - s, w, 0, d};
Point(22) = {s1, w, 0, d};
Point(23) = {s1 - s, w, 0, d};
Point(24) = {0, w, 0, d};


Line(1) = {24, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 17};
Line(18) = {17, 18};
Line(19) = {18, 19};
Line(20) = {19, 20};
Line(21) = {20, 21};
Line(22) = {21, 22};
Line(23) = {22, 23};
Line(24) = {23, 24};
Line(25) = {2, 23};
Line(26) = {3, 22};
Line(27) = {4, 21};
Line(28) = {5, 20};
Line(29) = {6, 19};
Line(30) = {7, 18};
Line(31) = {8, 17};
Line(32) = {9, 16};
Line(33) = {10, 15};
Line(34) = {11, 14};
// ======================================================
// Surface
// ======================================================

Curve Loop(1) = {1, 2, 25, 24};
Curve Loop(2) = {-25, 3, 26, 23};
Curve Loop(3) = {4, 27, 22, -26};
Curve Loop(4) = {5, 28, 21, -27};
Curve Loop(5) = {6, 29, 20, -28};
Curve Loop(6) = {7, 30, 19, -29};
Curve Loop(7) = {8, 31, 18, -30};
Curve Loop(8) = {9, 32, 17, -31};
Curve Loop(9) = {10, 33, 16, -32};
Curve Loop(10) = {11, 34, 15, -33};
Curve Loop(11) = {12, 13, 14, -34};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11}; 

Physical Surface(1) = {1,2,3,4,5,7,8,9,10,11};
Physical Surface(2) = {6};

Transfinite Curve {1,25,26,27,28,29,30,31,32,33,34,13}=11;
Transfinite Curve {2,24,12,14} = 4;
Transfinite Curve {3,23,11,15} = 4;
Transfinite Curve {5,21,9,17} = 4;
Transfinite Curve {-6,20,-18,8} = 30 Using Progression 1.05;
Transfinite Curve {7,19} = 9;

Transfinite Surface {1} = {24, 1, 2, 23};
Transfinite Surface {2} = {23, 2, 3, 22};
Transfinite Surface {3} = {22, 3, 4, 21};
Transfinite Surface {4} = {21, 4, 5, 20};
Transfinite Surface {5} = {20, 5, 6, 19};
Transfinite Surface {6} = {19, 6, 7, 18};
Transfinite Surface {7} = {18, 7, 8, 17};
Transfinite Surface {8} = {17, 8, 9, 16};
Transfinite Surface {9} = {16, 9, 10, 15};
Transfinite Surface {10} = {15, 10, 11, 14};
Transfinite Surface {11} = {14, 11, 12, 13};
Recombine Surface {1,2,3,4,5,6,7,8,9,10,11};

Physical Curve(1) = {3};
Physical Curve(2) = {11};
Physical Curve(3) = {17};
Physical Curve(4) = {21};


// ======================================================
// Output
// ======================================================

Mesh.MshFileVersion = 2.2;
