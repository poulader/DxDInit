#include "EPGeometryGenerator.h"
#include <vector>
#include <Windows.h>
#include <d3d11.h>
#include <xnamath.h>

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