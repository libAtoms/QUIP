% AX_3D_CYLINDER_MESH = 24; %
AX_3D_CYLINDER_MESH = 12;
theta = 2*pi/AX_3D_CYLINDER_MESH;
for i = 1:AX_3D_CYLINDER_MESH,
 fprintf(1,'%f,%f, ',cos(i*theta),sin(i*theta));
end
