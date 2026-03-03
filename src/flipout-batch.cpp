// flipout-batch.cpp - Batch geodesic computation using FlipOut and MMP algorithms
// For comparison with Iterative Funnel Algorithm
//
// Usage: flipout-batch mesh.obj pairs.json output.json
//
// Input pairs.json format:  [{"from": 0, "to": 7}, ...]
// Output format: {"flipout_time_ms": N, "mmp_time_ms": M, "results": [{"from": 0, "to": 7, "distance": 1.234, "mmp_distance": 1.230, "waypoints": [...]}, ...]}

#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <nlohmann/json.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

using namespace geometrycentral;
using namespace geometrycentral::surface;
using json = nlohmann::json;

void printUsage(const char* progName) {
    std::cerr << "Usage: " << progName << " <mesh.obj> <pairs.json> <output.json>\n";
    std::cerr << "\n";
    std::cerr << "Input pairs.json format:\n";
    std::cerr << "  [{\"from\": 0, \"to\": 7}, {\"from\": 10, \"to\": 20}, ...]\n";
    std::cerr << "\n";
    std::cerr << "Output format:\n";
    std::cerr << "  {\"flipout_time_ms\": N, \"mmp_time_ms\": M, \"results\": [{\"from\": 0, \"to\": 7, \"distance\": 1.234, \"mmp_distance\": 1.230, \"waypoints\": [...]}, ...]}\n";
}

int main(int argc, char** argv) {
    if (argc != 4) {
        printUsage(argv[0]);
        return 1;
    }

    std::string meshFile = argv[1];
    std::string pairsFile = argv[2];
    std::string outputFile = argv[3];

    // Load mesh
    std::cout << "Loading mesh: " << meshFile << std::endl;
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;

    try {
        std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshFile);
    } catch (const std::exception& e) {
        std::cerr << "Error loading mesh: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "  Vertices: " << mesh->nVertices() << std::endl;
    std::cout << "  Faces: " << mesh->nFaces() << std::endl;

    // Load pairs
    std::cout << "Loading pairs: " << pairsFile << std::endl;
    json pairs;
    {
        std::ifstream inStream(pairsFile);
        if (!inStream) {
            std::cerr << "Error: could not open " << pairsFile << std::endl;
            return 1;
        }
        try {
            inStream >> pairs;
        } catch (const std::exception& e) {
            std::cerr << "Error parsing JSON: " << e.what() << std::endl;
            return 1;
        }
    }

    if (!pairs.is_array()) {
        std::cerr << "Error: pairs.json must be an array" << std::endl;
        return 1;
    }

    std::cout << "  Pairs to process: " << pairs.size() << std::endl;

    // Validate all pairs first
    std::vector<std::pair<size_t, size_t>> validPairs;
    for (size_t i = 0; i < pairs.size(); i++) {
        const auto& pair = pairs[i];
        if (pair.find("from") == pair.end() || pair.find("to") == pair.end()) {
            std::cerr << "Warning: pair " << i << " missing 'from' or 'to' field, skipping" << std::endl;
            continue;
        }
        size_t fromIdx = pair["from"].get<size_t>();
        size_t toIdx = pair["to"].get<size_t>();
        if (fromIdx >= mesh->nVertices() || toIdx >= mesh->nVertices()) {
            std::cerr << "Warning: pair " << i << " has invalid vertex index, skipping" << std::endl;
            continue;
        }
        validPairs.push_back({fromIdx, toIdx});
    }

    std::cout << "  Valid pairs: " << validPairs.size() << std::endl;

    // ============================================================
    // PHASE 1: FlipOut computation (timed separately)
    // ============================================================
    std::cout << "\n=== Phase 1: FlipOut Geodesics ===" << std::endl;

    struct FlipOutResult {
        size_t from, to;
        double distance;
        json waypoints;
        int cornerCount;
        bool success;
    };
    std::vector<FlipOutResult> flipoutResults(validPairs.size());

    auto flipoutStart = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < validPairs.size(); i++) {
        size_t fromIdx = validPairs[i].first;
        size_t toIdx = validPairs[i].second;
        Vertex vStart = mesh->vertex(fromIdx);
        Vertex vEnd = mesh->vertex(toIdx);

        flipoutResults[i].from = fromIdx;
        flipoutResults[i].to = toIdx;
        flipoutResults[i].cornerCount = 0;
        flipoutResults[i].success = false;

        try {
            // Construct path from Dijkstra
            auto edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, vStart, vEnd);

            if (edgeNetwork == nullptr) {
                std::cerr << "Warning: FlipOut could not construct path for " << fromIdx << " -> " << toIdx << std::endl;
                continue;
            }

            edgeNetwork->posGeom = geometry.get();

            // Iteratively shorten to geodesic
            edgeNetwork->iterativeShorten();

            // Get the 3D polyline
            auto polylines = edgeNetwork->getPathPolyline3D();

            // Also get SurfacePoint version to count actual vertex corners
            auto surfacePointPaths = edgeNetwork->getPathPolyline();
            int cornerCount = 0;
            for (const auto& path : surfacePointPaths) {
                for (const auto& sp : path) {
                    if (sp.type == SurfacePointType::Vertex) {
                        cornerCount++;
                    }
                }
            }

            // Compute total distance and collect waypoints
            double distance = 0.0;
            json waypoints = json::array();

            for (const auto& path : polylines) {
                for (size_t j = 0; j < path.size(); j++) {
                    const Vector3& pt = path[j];
                    waypoints.push_back({pt.x, pt.y, pt.z});
                    if (j > 0) {
                        distance += (path[j] - path[j-1]).norm();
                    }
                }
            }

            flipoutResults[i].distance = distance;
            flipoutResults[i].waypoints = waypoints;
            flipoutResults[i].cornerCount = cornerCount;
            flipoutResults[i].success = true;

        } catch (const std::exception& e) {
            std::cerr << "Warning: FlipOut exception for pair " << i << ": " << e.what() << std::endl;
        }

        // Progress indicator
        if ((i + 1) % 100 == 0 || i + 1 == validPairs.size()) {
            std::cout << "  FlipOut: " << (i + 1) << "/" << validPairs.size() << std::endl;
        }
    }

    auto flipoutEnd = std::chrono::high_resolution_clock::now();
    auto flipoutDuration = std::chrono::duration_cast<std::chrono::milliseconds>(flipoutEnd - flipoutStart);

    size_t flipoutSuccess = 0;
    for (const auto& r : flipoutResults) if (r.success) flipoutSuccess++;
    std::cout << "  FlipOut completed: " << flipoutSuccess << " successful, " << flipoutDuration.count() << " ms" << std::endl;

    // ============================================================
    // PHASE 2: MMP (Exact) computation (timed separately)
    // ============================================================
    std::cout << "\n=== Phase 2: MMP Exact Geodesics ===" << std::endl;

    std::vector<double> mmpDistances(validPairs.size(), -1.0);

    // Require edge lengths for exact algorithm
    geometry->requireEdgeLengths();

    auto mmpStart = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < validPairs.size(); i++) {
        size_t fromIdx = validPairs[i].first;
        size_t toIdx = validPairs[i].second;
        Vertex vStart = mesh->vertex(fromIdx);
        Vertex vEnd = mesh->vertex(toIdx);

        try {
            // Create exact geodesic algorithm and propagate from source
            GeodesicAlgorithmExact exactAlgo(*mesh, *geometry);
            exactAlgo.propagate(vStart, GEODESIC_INF, {vEnd});

            // Get exact distance
            double exactDistance = exactAlgo.getDistance(vEnd);
            mmpDistances[i] = exactDistance;

        } catch (const std::exception& e) {
            std::cerr << "Warning: MMP exception for pair " << i << ": " << e.what() << std::endl;
        }

        // Progress indicator
        if ((i + 1) % 100 == 0 || i + 1 == validPairs.size()) {
            std::cout << "  MMP: " << (i + 1) << "/" << validPairs.size() << std::endl;
        }
    }

    auto mmpEnd = std::chrono::high_resolution_clock::now();
    auto mmpDuration = std::chrono::duration_cast<std::chrono::milliseconds>(mmpEnd - mmpStart);

    size_t mmpSuccess = 0;
    for (double d : mmpDistances) if (d >= 0) mmpSuccess++;
    std::cout << "  MMP completed: " << mmpSuccess << " successful, " << mmpDuration.count() << " ms" << std::endl;

    // ============================================================
    // Combine results
    // ============================================================
    json results = json::array();
    for (size_t i = 0; i < validPairs.size(); i++) {
        if (!flipoutResults[i].success) continue;

        json result = {
            {"from", flipoutResults[i].from},
            {"to", flipoutResults[i].to},
            {"distance", flipoutResults[i].distance},
            {"mmp_distance", mmpDistances[i]},
            {"corner_count", flipoutResults[i].cornerCount},
            {"waypoints", flipoutResults[i].waypoints}
        };
        results.push_back(result);
    }

    // Write output
    std::cout << "\nWriting results: " << outputFile << std::endl;
    {
        json output;
        output["flipout_time_ms"] = flipoutDuration.count();
        output["mmp_time_ms"] = mmpDuration.count();
        output["results"] = results;

        std::ofstream outStream(outputFile);
        if (!outStream) {
            std::cerr << "Error: could not open " << outputFile << " for writing" << std::endl;
            return 1;
        }
        outStream << output.dump(2);
    }

    // Summary
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "  Total pairs: " << validPairs.size() << std::endl;
    std::cout << "  FlipOut: " << flipoutSuccess << " successful, " << flipoutDuration.count() << " ms" << std::endl;
    std::cout << "  MMP:     " << mmpSuccess << " successful, " << mmpDuration.count() << " ms" << std::endl;

    // Quick accuracy check
    double maxDiff = 0, sumDiff = 0;
    size_t nCompared = 0;
    for (size_t i = 0; i < validPairs.size(); i++) {
        if (flipoutResults[i].success && mmpDistances[i] >= 0) {
            double diff = flipoutResults[i].distance - mmpDistances[i];
            double relDiff = std::abs(diff) / mmpDistances[i];
            maxDiff = std::max(maxDiff, relDiff);
            sumDiff += relDiff;
            nCompared++;
        }
    }
    if (nCompared > 0) {
        std::cout << "\n  FlipOut vs MMP accuracy:" << std::endl;
        std::cout << "    Compared: " << nCompared << " paths" << std::endl;
        std::cout << "    Mean relative diff: " << (sumDiff / nCompared * 100) << "%" << std::endl;
        std::cout << "    Max relative diff:  " << (maxDiff * 100) << "%" << std::endl;
    }

    return 0;
}
