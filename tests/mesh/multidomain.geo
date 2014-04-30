dx = 1;
dy = 1;
dz = 3;

R1 = 4.0;
R2 = 6.0;

lc = 1.5;

//lc_i   = 0.5;
//lc_s   = 0.5;

pio2=Pi/2;

//////////////////////////
// BAR MAGNET
//////////////////////////

Point(11) = { dx,  dy,  dz, lc};
Point(12) = { dx,  dy, -dz, lc};
Point(13) = { dx, -dy,  dz, lc};
Point(14) = { dx, -dy, -dz, lc};
Point(15) = {-dx,  dy,  dz, lc};
Point(16) = {-dx,  dy, -dz, lc};
Point(17) = {-dx, -dy,  dz, lc};
Point(18) = {-dx, -dy, -dz, lc};

Line(101) = {11, 12};
Line(102) = {12, 14};
Line(103) = {14, 13};
Line(104) = {13, 11};

Line(105) = {15, 16};
Line(106) = {16, 18};
Line(107) = {18, 17};
Line(108) = {17, 15};

Line(109) = {11, 15};
Line(110) = {12, 16};
Line(111) = {13, 17};
Line(112) = {14, 18};

Line Loop(11)     = {101,102,103,104};
Line Loop(12)     = {105,106,107,108};
Line Loop(13)     = {101,110,-105,-109};
Line Loop(14)     = {102,112,-106,-110};
Line Loop(15)     = {103,111,-107,-112};
Line Loop(16)     = {104,109,-108,-111};

Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};

Surface Loop(1)  = {11, 12, 13, 14, 15, 16};
Physical Surface(2)  = {11, 12, 13, 14, 15, 16};

//////////////////////////
// SHELL
//////////////////////////

// create inner 1/8 shell
Point(41) = {0, 0, 0, lc};
Point(42) = {-R1, 0, 0, lc};
Point(43) = {0, R1, 0, lc};
Point(46) = {0, 0, R1, lc};

Circle(21) = {42, 41, 43};
Circle(25)= {46,41,42};
Circle(26)= {46,41,43};

Line Loop(10) = {21, -26, 25} ;
Ruled Surface (60) = {10};

// create outer 1/8 shell
t0[] = {60};

// create remaining 7/8 inner shells
t1[] = Rotate {{0,0,1},{0,0,0},pio2}   {Duplicata{Surface{t0[]};}};
t2[] = Rotate {{0,0,1},{0,0,0},pio2*2} {Duplicata{Surface{t0[]};}};
t3[] = Rotate {{0,0,1},{0,0,0},pio2*3} {Duplicata{Surface{t0[]};}};
t4[] = Rotate {{0,1,0},{0,0,0},-pio2}  {Duplicata{Surface{t0[]};}};
t5[] = Rotate {{0,0,1},{0,0,0},pio2}   {Duplicata{Surface{t4[]};}};
t6[] = Rotate {{0,0,1},{0,0,0},pio2*2} {Duplicata{Surface{t4[]};}};
t7[] = Rotate {{0,0,1},{0,0,0},pio2*3} {Duplicata{Surface{t4[]};}};

// create entire inner and outer shell
Surface Loop(100)={t0[0],t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
//Physical Surface(1) = {t0[0],t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
Physical Surface(3) = {t0[0],t1[0],t2[0],t3[0]};
Physical Surface(4) = {t7[0],t4[0],t5[0],t6[0]};

//////////////////////////
// VOLUMES
//////////////////////////

Volume(1000)={1};       // Bar magnet
Volume(2000)={100,1};   // inner sphere

// Physical Volumes
Physical Volume(1) = {1000};
Physical Volume(2) = {2000};
