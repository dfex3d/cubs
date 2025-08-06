// Terrain Generation Web Worker
// This worker handles the computationally expensive marching cubes algorithm
// to prevent blocking the main thread during terrain sculpting

// Import marching cubes lookup tables
importScripts('marching_tables.js');

console.log('Terrain worker loaded, edgeTable length:', edgeTable.length);

// Terrain generation constants (must match main thread)
const GRID_SIZE = 64;
const GRID_EXTENT = 30; // planetRadius * 1.5, will be updated from main thread
const ISO_LEVEL = 0.0;

let VOXEL_SIZE = (GRID_EXTENT * 2) / GRID_SIZE;

// Voxel indexing function
function voxelIndex(x, y, z) {
  const stride = GRID_SIZE + 1;
  return x + y * stride + z * stride * stride;
}

// Main terrain mesh generation function (moved from main thread)
function generateTerrainMesh(voxelData, gridSize, gridExtent, isoLevel) {
  const positions = [];
  const normals = [];
  const uvs = [];
  const colors = [];
  
  // Update local constants
  const voxelSize = (gridExtent * 2) / gridSize;
  
  // Marching cubes vertex and edge definitions
  const vertexOffset = [
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
  ];
  
  const edgeConnection = [
    [0,1],[1,2],[2,3],[3,0],
    [4,5],[5,6],[6,7],[7,4],
    [0,4],[1,5],[2,6],[3,7]
  ];
  
  const edgeVerts = new Array(12);
  const cubeVerts = new Array(8);
  const cubeValues = new Array(8);

  // Process each cube in the voxel grid
  for (let z = 0; z < gridSize; z++) {
    for (let y = 0; y < gridSize; y++) {
      for (let x = 0; x < gridSize; x++) {
        let cubeIndex = 0;
        
        // Sample voxel values at cube corners
        for (let i = 0; i < 8; i++) {
          const vx = x + vertexOffset[i][0];
          const vy = y + vertexOffset[i][1];
          const vz = z + vertexOffset[i][2];
          const value = voxelData[voxelIndex(vx, vy, vz)];
          cubeValues[i] = value;
          
          // Convert to world coordinates
          const wx = vx * voxelSize - gridExtent;
          const wy = vy * voxelSize - gridExtent;
          const wz = vz * voxelSize - gridExtent;
          cubeVerts[i] = { x: wx, y: wy, z: wz };
          
          // Build cube index for marching cubes lookup
          if (value < isoLevel) cubeIndex |= 1 << i;
        }
        
        // Get edge intersections from lookup table
        const edges = edgeTable[cubeIndex];
        if (!edges) continue;
        
        // Calculate edge intersection points
        for (let i = 0; i < 12; i++) {
          if (edges & (1 << i)) {
            const [a0, b0] = edgeConnection[i];
            const valp1 = cubeValues[a0];
            const valp2 = cubeValues[b0];
            const p1 = cubeVerts[a0];
            const p2 = cubeVerts[b0];
            
            // Linear interpolation to find surface intersection
            const mu = (isoLevel - valp1) / (valp2 - valp1);
            edgeVerts[i] = {
              x: p1.x + mu * (p2.x - p1.x),
              y: p1.y + mu * (p2.y - p1.y),
              z: p1.z + mu * (p2.z - p1.z)
            };
          }
        }
        
        // Generate triangles from lookup table
        const tri = triTable[cubeIndex];
        for (let i = 0; i < 16 && tri[i] !== -1; i += 3) {
          const a = edgeVerts[tri[i]];
          const b = edgeVerts[tri[i+1]];
          const c = edgeVerts[tri[i+2]];
          
          // Add triangle vertices (reversed winding order)
          positions.push(c.x, c.y, c.z, b.x, b.y, b.z, a.x, a.y, a.z);

          // Calculate face normal
          const cb = { x: b.x - c.x, y: b.y - c.y, z: b.z - c.z };
          const ca = { x: a.x - c.x, y: a.y - c.y, z: a.z - c.z };

          // Cross product for normal
          const normal = {
            x: cb.y * ca.z - cb.z * ca.y,
            y: cb.z * ca.x - cb.x * ca.z,
            z: cb.x * ca.y - cb.y * ca.x
          };

          // Normalize
          const length = Math.sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
          if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
          }

          // Add normal for each vertex of the triangle
          for (let j = 0; j < 3; j++) {
            normals.push(normal.x, normal.y, normal.z);
          }

          // Generate UV coordinates based on world position
          // Use spherical mapping for better texture distribution
          const vertices = [c, b, a];
          for (let j = 0; j < 3; j++) {
            const vertex = vertices[j];
            const radius = Math.sqrt(vertex.x * vertex.x + vertex.y * vertex.y + vertex.z * vertex.z);
            const theta = Math.atan2(vertex.z, vertex.x);
            const phi = Math.acos(vertex.y / radius);

            const u = (theta + Math.PI) / (2 * Math.PI);
            const v = phi / Math.PI;

            uvs.push(u, v);
          }

          // Generate random color variation for each triangle
          // Create subtle variation based on cube position for consistency
          const cubeHash = (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
          const randomSeed = (cubeHash % 1000) / 1000.0;

          // Small color variation (±5% brightness)
          const colorVariation = 0.9 + (randomSeed * 0.2);

          // Apply same color variation to all vertices of this triangle
          for (let j = 0; j < 3; j++) {
            colors.push(colorVariation, colorVariation, colorVariation);
          }
        }
      }
    }
  }
  
  return {
    positions: new Float32Array(positions),
    normals: new Float32Array(normals),
    uvs: new Float32Array(uvs),
    colors: new Float32Array(colors)
  };
}

// Worker message handler
self.onmessage = function(e) {
  const { type, data } = e.data;
  
  switch (type) {
    case 'generateTerrain':
      const { voxelData, gridSize, gridExtent, isoLevel } = data;

      console.log('Worker: Starting terrain generation, voxel data length:', voxelData.length);

      // Generate mesh data
      const meshData = generateTerrainMesh(voxelData, gridSize, gridExtent, isoLevel);

      console.log('Worker: Terrain generation complete, positions:', meshData.positions.length, 'normals:', meshData.normals.length);

      // Send result back to main thread
      self.postMessage({
        type: 'terrainGenerated',
        data: meshData
      });
      break;
      
    default:
      console.warn('Unknown message type:', type);
  }
};

// Signal that worker is ready
self.postMessage({ type: 'workerReady' });
