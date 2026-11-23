// Sphere de test de regression (2026-11-20) : diametre 17 cm, zone fine a l'oreille
// DROITE (axe -y) : hmin dans un rayon r1 autour de (0,-a,0), transition lineaire
// jusqu'a hmax a r2, hmax ailleurs. Unites : metres. Normales vers le fluide (gmsh).
//   gmsh sphere_test.geo -2 -format msh2 -o sphere_test.msh
SetFactory("OpenCASCADE");
DefineConstant[ a = 0.085, hmin = 0.002, hmax = 0.008, r1 = 0.010, r2 = 0.040 ];
Sphere(1) = {0, 0, 0, a};
Physical Surface(1) = Surface{:};
Field[1] = MathEval;
Field[1].F = Sprintf("Min(%g, Max(%g, %g + (Sqrt(x^2+(y+%g)^2+z^2) - %g)*(%g - %g)/(%g - %g)))", hmax, hmin, hmin, a, r1, hmax, hmin, r2, r1);
Background Field = 1;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.Algorithm = 6;
Mesh.RandomSeed = 1;
