#pragma once

#include <vector>
#include <Windows.h>
#include <d3d11.h>
#include <xnamath.h>

//put this in a different namespace than Luna's
namespace EPGeometry
{

	//Class to generate various meshes of geometry
	class EPGeometryGenerator
	{
	public:

		typedef struct Vertex
		{
			Vertex() {};

			Vertex(const XMFLOAT3& p,
				const XMFLOAT3& n,
				const XMFLOAT3& t,
				const XMFLOAT2& uv)
				: Position(p), Normal(n), TangentU(t), TexC(uv)
			{

			};

			Vertex(float px, float py, float pz,
				float nx, float ny, float nz,
				float tx, float ty, float tz,
				float u, float v)
				: Position(px, py, pz), Normal(nx, ny, nz), TangentU(tx, ty, tz), TexC(u, v)
			{

			};


			XMFLOAT3 Position;
			XMFLOAT3 Normal;
			XMFLOAT3 TangentU;
			XMFLOAT2 TexC;

		} Vertex, *pVertex;

		typedef struct MeshData
		{

			std::vector<Vertex> Vertices;
			std::vector<UINT> Indices;

		} MeshData, *pMeshData;

		static HRESULT CreateGrid(float width, float depth, UINT m, UINT n, MeshData& meshData);



	};
}