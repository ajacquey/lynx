L = 6;
r = 1;
lc_i = 5.0e-02;
lc_m = 5.0e-01;

Point(1) = {-L/2., -L/2., 0.0, lc_m};
Point(2) = {L/2., -L/2., 0.0, lc_m};
Point(3) = {L/2., L/2., 0.0, lc_m};
Point(4) = {-L/2., L/2., 0.0, lc_m};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(5) = {0.0, 0.0, 0.0, lc_m};
Point(6) = {r, 0.0, 0.0, lc_i};
Point(7) = {0.0, r, 0.0, lc_i};
Point(8) = {-r, 0.0, 0.0, lc_i};
Point(9) = {0.0, -r, 0.0, lc_i};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1,2};
Plane Surface(2) = {2};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};