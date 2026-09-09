---
oneliner: C++ demo turning paths on a triangle mesh into geodesics by iteratively flipping edges
tags: [geometry-processing, geodesics, triangle-mesh, edge-flips, computational-geometry, c++, geometry-central, polyscope, gui, siggraph]
stack: [C++, CMake, geometry-central, Polyscope]
generated: 2026-09-09
commit: 48dfb1b
placeholder: false
---
Reference implementation and GUI for the SIGGRAPH Asia 2020 paper "You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges" (Sharp & Crane). Loads a mesh and an initial path (via Dijkstra, manual point picking, or file loaders), then runs the FlipOut algorithm to straighten it into a locally-shortest geodesic, with a Polyscope-based UI for interactive input and visualization. Working demo, built on the geometry-central library.
