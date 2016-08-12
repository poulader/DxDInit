#include "EPGeometryGenerator.h"
#include <vector>
#include <Windows.h>
#include <d3d11.h>
#include <xnamath.h>
#include "MathHelper.h"

using namespace EPGeometry;

//TODO: Add overload allowing axis to be specified
//width = width in cells along x axis, depth = depth in cells along z axis
//m = number of vertices along a line in the z axis, n = number of vertices in a line along x axis
//meshdata = stores results
HRESULT EPGeometryGenerator::CreateGrid(float width, float depth, UINT m, UINT n, MeshData& meshData)
{

	if (width <= 0 || depth <= 0)
		return -1;
	else if (m < 2 || n < 2)
		return -2;

	//we want to make a grid of width cells by depth cells, width = span along x axis, depth = span along z axis, in local space

	//In order to calculate the vertex coordinates, we work in units of cells rather than vertices to make it simpler.

	//apperently n = x axis, m = z axis

	//total count of vertices
	UINT vertexCount = m*n;
	//total count of triangles to cover grid
	UINT faceCount = 2 * (m - 1) * (n - 1);

	//length of each cell in the x direction, in local space units
	float deltaX = width / (n-1);

	//length of e ach cell in the z direction, in local space units
	float deltaZ = depth / (m-1);

	//half of the width, i.e. distance in x axis from the origin to a corner
	float halfWidth = width / 2;

	//half of the depth, i.e. distance in z axis from the origin to a corner
	float halfDepth = depth / 2;

	//used for texturing, looks like its the percentage each cell takes in each direction
	//wtf, it really looks like he mixed up his x and y axises, may need to change later
	float du = 1.0f / (n - 1);
	float dv = 1.0f / (m - 1);

	//We will have m * n vertices
	meshData.Vertices.resize(vertexCount);

	UINT facePerQuad = 2;

	UINT indicesPerQuad = facePerQuad * 3;

	//We will have 2 faces per quad, for a total of 2*3*(m-1)*(n-1) indices
	meshData.Indices.resize(faceCount * 3);

	//generate vertices, start at +z, -x, move in the +x direction until end is hit, then lower z by 1. Repeat until finished.

	/*
	
	
	.	.	.	.   |	* -------->
	.	.	.	.	z	* -------->
	.	.	.	.	|   * -------->
	.	.	.	.
	.	.	.	.
	_____ x _____
	
	*/

	for (UINT i = 0; i < m; i++)
	{

		float zCoord = halfDepth - i*deltaZ;

		for (UINT j = 0; j < n; j++)
		{
			float xCoord = -halfWidth + j*deltaX;

			meshData.Vertices[i*n + j].Position = XMFLOAT3(xCoord, 0, zCoord);

			//Ignore for now, used for lighting
			meshData.Vertices[i*n + j].Normal = XMFLOAT3(0.0f, 1.0f, 0.0f);
			meshData.Vertices[i*n + j].TangentU = XMFLOAT3(1.0f, 0.0f, 0.0f);

			//Ignore for now, used for texturing
			meshData.Vertices[i*n + j].TexC.x = j*du;
			meshData.Vertices[i*n + j].TexC.y = i*dv;


			//Cant do indices until all the vertices are defined
			//well actually we can, if we lag one behind. But that would introduce a branch into the inner loop,
			//probably faster to just run entire loop again without branch.


		}

	}

	//iterate over each cell. There are (m-1)*(n-1) cells.

	//current row
	for (UINT i = 0; i < m - 1; i++)
	{

		//current column
		for (UINT j = 0; j < n - 1; j++)
		{
			UINT baseIndexIntoArray = (i*(n-1) + j) * 6;
			meshData.Indices[baseIndexIntoArray++] = i*n + j;	
			meshData.Indices[baseIndexIntoArray++] = i*n + j + 1;
			meshData.Indices[baseIndexIntoArray++] = (i+1)*n + j;
			meshData.Indices[baseIndexIntoArray++] = (i+1)*n + j;
			meshData.Indices[baseIndexIntoArray++] = i*n + j + 1;
			meshData.Indices[baseIndexIntoArray] =	(i+1)*n + j + 1;

		}

	}


	return 0;


};


//create a cylinder with specified top and bottom radius, height, slices, stackCount
HRESULT EPGeometryGenerator::CreateCylinder(UINT sliceCount, UINT stackCount, float height, float topRadius, float bottomRadius, MeshData& meshData)
{
	//We need at least four slices or we cannot make triangles using this method
	//we need at least one stackcount or we cannot make a cylinder
	//In this case, height >= 0
	//top and bottom radius >= 0
	if (sliceCount < 3 || stackCount < 1 || height <= 0 || topRadius <= 0 || bottomRadius <= 0)
	{
		return -1;
	}

	meshData.Vertices.clear();
	meshData.Indices.clear();
	
	//Number of rings
	UINT nRings = stackCount + 1;

	//deltaRadius (from bottom going up) = (topRadius - bottomRadius) / stackCount
	float deltaRadius = (topRadius - bottomRadius) / (float)stackCount;

	//height difference per ring = height / stacjCOunt
	float deltaHeight = height / (float)stackCount;

	float baseHeight = height / -2.0f;

	float deltaAnglePerSlice = (2 * XM_PI) / sliceCount;

	//using lots of variables so its easier to debug/trace
	for (UINT i = 0; i < nRings; ++i)
	{

		float ringRadius = bottomRadius + i*deltaRadius;
		float ringHeight = baseHeight + i*deltaHeight;

		//duplicate vertex for first and last position, but texture coords differ
		for (UINT j = 0; j <= sliceCount; ++j)
		{
			Vertex vertex;
			//we have an angle and a radius, convert from polar coords in the current plane, height is constant per ring

			float c = cosf(j*deltaAnglePerSlice);
			float s = sinf(j*deltaAnglePerSlice);
			float y = ringHeight;

			vertex.Position = XMFLOAT3(ringRadius * c, y, ringRadius*s);

			vertex.TexC.x = (float)j / sliceCount;
			vertex.TexC.y = 1.0f - (float)i / stackCount;

			//Parameterize cylinder, where v goes in same direction as v tex-coord, bitangent goes in same direction as v tex coord
			//let r0 be the bottom radius and let r1 be the top radius
			// y(v) = h - hv for v in [0,1]
			// r(v) = r1 + (r0 - r1)*v
			//
			// x(t, v) = r(v)*cos(t)
			// y(t, v) = h - hv
			// z(t, v) = r(v)*sin(t)

			// dx/dt = -r(v)*sin(t)
			// dy/dt = 0
			// dz/dt = r(v)*cos(t)

			// dx/dv = (r0 - r1) *cos(t)
			// dy/dv = -h
			// dz/dv = (r0 -r1)*sin(t)

			//TangentU use unit length
			
			vertex.TangentU = XMFLOAT3(-s, 0.0f, c);

			float dr = bottomRadius - topRadius;

			XMFLOAT3 bitangent(dr*c, -y, dr*s);

			XMVECTOR T = XMLoadFloat3(&vertex.TangentU);
			XMVECTOR B = XMLoadFloat3(&bitangent);
			XMVECTOR N = XMVector3Normalize(XMVector3Cross(T, B));
			XMStoreFloat3(&vertex.Normal, N);


			meshData.Vertices.push_back(vertex);
		}

	}

	//Generate indices, iterate over each slice, index as j = slice, i = ring

	//Note that we need to add one to ringVertexCount as first and last vertex per ring is duplicated with different texcoords
	UINT ringVertexCount = sliceCount + 1;

	for (UINT i = 0; i < stackCount; ++i)
	{

		for (UINT j = 0; j < sliceCount; j++)
		{
			//i*ringVertexCount = index of first vertex on bottom ring of stack, j = offset from the index 0 based ring
			//Note that the final pass in each slice does index to the "extra" vertex, which makes sense as that's the only way
			//it would be fed to the GPU. As the position is the same, it does not alter the faces. However, the texture coords & maybe normal, bitangent, are different.
			meshData.Indices.push_back(i*ringVertexCount + j);
			meshData.Indices.push_back((i + 1)*ringVertexCount + j);
			meshData.Indices.push_back((i + 1)*ringVertexCount + j + 1);

			meshData.Indices.push_back(i*ringVertexCount + j);
			meshData.Indices.push_back((i + 1)*ringVertexCount + j + 1);
			meshData.Indices.push_back(i*ringVertexCount + j + 1);
		}

	}

	//now we need to built the caps, each function will insert an extra vertex at the end of the vertex buffer for the "center".
	BuildCylinderTopCap(bottomRadius, topRadius, height, sliceCount, meshData);
	BuildCylinderBottomCap(bottomRadius, topRadius, height, sliceCount, meshData);

	return 0;

}


//private functions, dont bother validating params for now, add later

HRESULT EPGeometryGenerator::BuildCylinderTopCap(float bottomRadius, float topRadius, float height, UINT sliceCount, MeshData& meshData)
{

	UINT baseIndex = (UINT)meshData.Vertices.size();

	float y = 0.5f*height;
	float s = 0, c = 0;

	float deltaTheta = 2 * XM_PI / sliceCount;

	for (UINT i = 0; i <= sliceCount; ++i)
	{
		Vertex vertex;
		s = sinf(deltaTheta*i);
		c = cosf(deltaTheta*i);

		float x, z;
		x = topRadius*c;
		z = topRadius*s;

		//scale down by the height to try and make top cap texture coord area proportional to base.
		float u = x / height + 0.5f;
		float v = z / height + 0.5f;

		meshData.Vertices.push_back(Vertex( z, y, z, 
											0.0f, 1.0f, 0.0f, 
											1.0f, 0.0f, 0.0f, 
											u, v));

	}

	//cap center vertex
	meshData.Vertices.push_back(Vertex( 0.0f, y, 0.0f,
										0.0f, 1.0f, 0.0f,
										1.0f, 0.0f, 0.0f,
										0.5f, 0.5f));

	//index of center vertex
	UINT centerIndex = (UINT)meshData.Vertices.size() - 1;

	for (UINT i = 0; i < sliceCount; ++i)
	{
		meshData.Indices.push_back(centerIndex);
		meshData.Indices.push_back(baseIndex + i +1);
		meshData.Indices.push_back(baseIndex + i);
	}

	return 0;

}

//private functions, dont bother validating params for now, add later

HRESULT EPGeometryGenerator::BuildCylinderBottomCap(float bottomRadius, float topRadius, float height, UINT sliceCount, MeshData& meshData)
{
	UINT baseIndex = (UINT)meshData.Vertices.size();

	float y = -0.5f*height;
	float s = 0, c = 0;

	float deltaTheta = 2 * XM_PI / sliceCount;

	for (UINT i = 0; i <= sliceCount; ++i)
	{
		Vertex vertex;
		s = sinf(deltaTheta*i);
		c = cosf(deltaTheta*i);

		float x, z;
		x = topRadius*c;
		z = topRadius*s;

		//scale down by the height to try and make top cap texture coord area proportional to base.
		float u = x / height + 0.5f;
		float v = z / height + 0.5f;

		meshData.Vertices.push_back(Vertex(z, y, z,
			0.0f, -1.0f, 0.0f,
			1.0f, 0.0f, 0.0f,
			u, v));

	}

	//cap center vertex
	meshData.Vertices.push_back(Vertex(0.0f, y, 0.0f,
		0.0f, -1.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
		0.5f, 0.5f));

	//index of center vertex
	UINT centerIndex = (UINT)meshData.Vertices.size() - 1;

	for (UINT i = 0; i < sliceCount; ++i)
	{
		meshData.Indices.push_back(centerIndex);
		meshData.Indices.push_back(baseIndex + i);
		meshData.Indices.push_back(baseIndex + i + 1);
	}

	return 0;

}

HRESULT EPGeometryGenerator::CreateSphere(UINT slices, UINT stackCount, float radius, MeshData& meshData)
{

	meshData.Vertices.clear();
	meshData.Indices.clear();

	//Compute vertices starting at top pole and moving down the stacks.
	
	// Poles: note that there will be texture coordinate distortion as there is
	// not a unique point on the texture map to assign to the pole when mapping
	// a rectangular texture onto a sphere.

	Vertex topVertex(0.0f, radius, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	Vertex bottomVertex(0.0f, -radius, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);

	meshData.Vertices.push_back(topVertex);

	//vertical
	float phiStep = XM_PI / stackCount;
	//sweep in xz plane
	float thetaStep = 2 * XM_PI / slices;

	//dont count poles as rings, ringcount then equals stackcount - 1, and we start at first ring not counting north pole, end at last ring not counting south pole
	for (UINT i = 1; i < stackCount ; ++i)
	{

		float ringPhi = i*phiStep;

		//just like cylinder, we repeat the first and last vertex position, but normals and textcoords are different
		for (UINT j = 0; j <= slices; ++j)
		{

			float sliceTheta = j*thetaStep;

			Vertex vertex;

			//coords are based on spherical to cartesian
			vertex.Position.x = sinf(ringPhi)*radius*cosf(sliceTheta);
			vertex.Position.y = cosf(ringPhi)*radius;
			vertex.Position.z = sinf(ringPhi)*radius*sinf(sliceTheta);

			//partial derivative of position w/ respect to theta (P sub(theta), i.e. the vector representing the rate at which position is changing at our theta)
			vertex.TangentU.x = sinf(ringPhi)*radius*-sinf(sliceTheta);
			vertex.TangentU.y = 0;
			vertex.TangentU.z = sinf(ringPhi)*radius*cosf(sliceTheta);

			//use SIMD registers to compute unit tangent and store
			XMVECTOR T = XMLoadFloat3(&vertex.TangentU);
			XMStoreFloat3(&vertex.TangentU, XMVector3Normalize(T));

			//use SIMD registers to compute unit normal (inherently the unit position vector in a sphere) and store
			XMVECTOR N = XMLoadFloat3(&vertex.Position);
			XMStoreFloat3(&vertex.Normal, XMVector3Normalize(N));

			vertex.TexC.x = sliceTheta / XM_2PI;
			vertex.TexC.y = ringPhi / XM_PI;

			meshData.Vertices.push_back(vertex);

		}


	}

	//push bottom pole
	meshData.Vertices.push_back(bottomVertex);

	//compute indices for top stack. The pole was pushed before everything else.

	//so lets think... what will the indices look like
	/*
		
		vertices ....   top pole      0
						ring 1
						ring ...
						ring stackcount - 2
						ring stackcount -1
						bottom pole   vertices.size()-1

		indices = stackcount
	*/
	UINT topPoleIndex = 0, bottomPoleIndex = meshData.Vertices.size() - 1;

	//seems like the best way to do this is manually do the first and last stacks, loop for the middle stacks.
	//I'm going to do it this way and then see how Luna did it. Better to ahve three separate loops instead of one giant loop with
	//possible conditional branches in the body.

	//first stack
	//do the first ring indices, we can reach 1 beyond slices-1 because there are duplicate vertices at start and end
	//not really sure where to start, there are three different points on each triangle we can start with to make clockwise winding order sets of vertices. Don't think it matters
	//through, because texture coords, tangent, normal, binormal, are embedded in each vertex structure
	UINT startRingBase = 1;

	//TO FIX: I fucked up the winding order calculations both here and in cylinder

	//yeh this is the same as Luna he just goes from i = 1 to i <= slices
	for (UINT i = 0; i < slices; i++)
	{
		//push top pole index
		meshData.Indices.push_back(topPoleIndex);
		//push the vertex back which would be equiv to cos(thetastep * i) (note that we skip the final duplicate pos vertex)
		meshData.Indices.push_back(startRingBase + i + 1);
		//push the vertex back which would be equiv to cos(thetastep*i + 1) (which will include the final duplicate pos vertex)
		meshData.Indices.push_back(startRingBase +i);
	}


	//all the other stacks - minus final stack which includes south pole
	//so looking at the stacks, we need to start at the top ring of a slice, have access to the next stack index
	//for the bottom ring of the slice. Basically (stackCount - 2) + 1 rings, starting at the offset into
	//the first ring vertex indices. Starting at i = 0, which points to the first set of ring vertices in the vertex array, we will
	//in the body of the loop, end up acessing (i*stackCount+1) to account for duplicate pos vertex of the next ring down to form the triangle.
	//so if (stackCount -2) +1 rings, i starts at 0, goes through stackCount -2 so we will be able, in final pass, to reference stackCount -2 +1.

	//yep same idea as luna, he just uses offset + scale * index
	for (UINT i = 0; i < stackCount - 2; i++)
	{

		//index of the first vertex of the top ring of the slice (remember we have sliceCount +1 vertices per ring)
		UINT startRingVertexIndex = startRingBase + i*(slices+1);

		//index of the first vertex of the bottom ring of the slice
		UINT stopRingVertexIndex = startRingVertexIndex + (slices + 1);

		//each ring - again note the final duplicate vertex pos (with differing tangent & texture coords) is included in the vertices of the final triangle.
		// A	B		*	B
		// C	*		C	D   
		//where B is the first point formed by the end of a ray between the pole and the cos(thetaStep*j) for the current ring
		for (UINT j = 0; j < slices; j++)
		{
			//First triance, ABC
			meshData.Indices.push_back(startRingVertexIndex + j );
			meshData.Indices.push_back(startRingVertexIndex + j + 1);
			meshData.Indices.push_back(stopRingVertexIndex + j);

			//second triangle, CBD
			meshData.Indices.push_back(stopRingVertexIndex + j);
			meshData.Indices.push_back(startRingVertexIndex + j + 1);
			meshData.Indices.push_back(stopRingVertexIndex + j + 1);
		}

	}


	//OK, that should be it for the middle stacks
	//now lets do the bottom stack

	//same idea as luna, again he uses the var which equals slices +1

	//first vertex of the ring forming the top of the last stack
	UINT bottomStackTopRingStartIndex = bottomPoleIndex - (slices + 1);

	for (UINT i = 0; i < slices; ++i)
	{
		meshData.Indices.push_back(bottomStackTopRingStartIndex + i);
		meshData.Indices.push_back(bottomStackTopRingStartIndex + i + 1);
		meshData.Indices.push_back(bottomPoleIndex);
	}

	return 0;

}