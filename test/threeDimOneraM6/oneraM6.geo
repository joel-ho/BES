////////////////////////////////////////////////////////////////////////////////
// oneraM6.geo
//	Creator: Ho Mun Onn, Joel
//
// Onera M6 wing generated from STEP file found on: 
// https://turbmodels.larc.nasa.gov/onerawingnumerics_val.html
////////////////////////////////////////////////////////////////////////////////

Merge "AileM6_with_sharp_TE.stp";
SetFactory("OpenCASCADE");
Dilate {{0, 0, 0}, {0.001, 0.001, 0.001}} {
  Volume{1}; 
}

rootChord = 805.9/1000;
tipChord = 0.56*rootChord;
span = 1196.3/1000;
xTip = 690.6841270315492/1000;
boundingRadius = 10*rootChord;

globalMaxMeshSize = 1.25*rootChord;
wingSurfMeshSize = rootChord/40;
wingEdgeMeshSize = rootChord/135;
wingSurfExternalMeshSize = rootChord;

wingEdgeAttractorRadius = 0.08*rootChord;
wingEdgeAttractorSurfaceMeshSize = wingEdgeMeshSize;
wingEdgeAttractorBoundaryMeshSize = wingSurfMeshSize;

wingTipAttractorRadius = 0.1*rootChord;
wingTipAttractorSurfaceMeshSize = wingEdgeMeshSize;
wingTipAttractorBoundaryMeshSize = wingSurfMeshSize;

wingSurfAttractorRadius = 2.75*rootChord;
wingSurfAttractorSurfaceMeshSize = wingSurfMeshSize;
wingSurfAttractorBoundaryMeshSize = wingSurfExternalMeshSize;

wingtipTEAttractorRadius = 0.15*rootChord;
wingtipTEAttractorMeshSize = wingEdgeMeshSize/3;
wingtipTEAttractorBoundaryMeshSize = wingSurfMeshSize;

// BooleanDifference { object } { tool }
Cylinder(2) = {0, 0, 0, 0, boundingRadius, 0, boundingRadius, 2*Pi};
BooleanDifference { Volume{2}; Delete; } { Volume{1}; Delete; } // Volume{2} is the new boolean-ed entity

Field[1] = Box;
Field[1].XMin = -(boundingRadius+1);
Field[1].XMax = boundingRadius+1;
Field[1].YMin = 0;
Field[1].YMax = boundingRadius+1;
Field[1].ZMin = -(boundingRadius+1);
Field[1].ZMax = boundingRadius;
Field[1].VIn = globalMaxMeshSize;
Field[1].VOut = boundingRadius;

Field[2] = Attractor;
Field[2].EdgesList = {25, 27};
Field[2].NNodesByEdge = 500;

Field[3] = MathEval;
Field[3].F = Sprintf("(%g-%g)*(F2/%g)^2 + %g", wingEdgeAttractorBoundaryMeshSize, wingEdgeAttractorSurfaceMeshSize, wingEdgeAttractorRadius, wingEdgeAttractorSurfaceMeshSize);

Field[4] = Attractor;
Field[4].EdgesList = {8, 10, 3};
Field[4].NNodesByEdge = 500;

Field[5] = MathEval;
Field[5].F = Sprintf("(%g-%g)*(F4/%g)^2 + %g", wingTipAttractorBoundaryMeshSize, wingTipAttractorSurfaceMeshSize, wingTipAttractorRadius, wingTipAttractorSurfaceMeshSize);

Field[6] = Attractor;
Field[6].FacesList = {11, 12, 13, 14};

Field[7] = MathEval;
Field[7].F = Sprintf("(%g-%g)*(F6/%g)^2 + %g", wingSurfAttractorBoundaryMeshSize, wingSurfAttractorSurfaceMeshSize, wingSurfAttractorRadius, wingSurfAttractorSurfaceMeshSize);

Field[8] = Attractor;
Field[8].NodesList = {2};

Field[9] = MathEval;
Field[9].F = Sprintf("(%g-%g)*(F8/%g)^2 + %g", wingtipTEAttractorBoundaryMeshSize, wingtipTEAttractorMeshSize, wingtipTEAttractorRadius, wingtipTEAttractorMeshSize);

Field[100] = Min;
Field[100].FieldsList = {1, 3, 5, 7, 9};
Background Field = 100;

Physical Surface ("FREESTREAM") = {8, 9};
Physical Surface ("SLIP_WALL") = {1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14};

Physical Volume ("FLUID") = {2};

Mesh.Algorithm = 6; // 1=MeshAdapt, 5=Delaunay, 6=Frontal
Mesh.CharacteristicLengthExtendFromBoundary = 0;
