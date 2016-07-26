#include <Windows.h>
#include <d3d11.h>
#include "DxAppBase.h"
#include "d3dUtil.h"
#include "MathHelper.h"
#include <stdint.h>
#include <tchar.h>
#include "d3dx11effect.h"
#include <xnamath.h>
#include "HillsDemo.h"
#include <iostream>
#include "EPGeometryGenerator.h"

HillsDemo::HillsDemo(HINSTANCE hWnd)
	: DxAppBase(hWnd), pHillsVertexBuffer(NULL), pHillsIndexBuffer(NULL), pFX(NULL), pTech(NULL),
	pfxWorldViewProj(NULL), pInputLayout(NULL), mTheta(1.5f * MathHelper::Pi), mPhi(0.25f * MathHelper::Pi), mRadius(5.0f), mGridIndexCount(0), mGridVertexCount(0)
{

	strMainWindowCaption = _T("Hills Demo");

	mLastMousePos.x = mLastMousePos.y = 0;

	XMMATRIX I = XMMatrixIdentity();

	//set our world, view, proj matrices to identity
	XMStoreFloat4x4(&mWorld, I);
	XMStoreFloat4x4(&mView, I);
	XMStoreFloat4x4(&mProj, I);
}

HillsDemo::~HillsDemo()
{
	ReleaseCOM(pHillsVertexBuffer);
	ReleaseCOM(pHillsIndexBuffer);
	ReleaseCOM(pFX);
	ReleaseCOM(pInputLayout);
}


bool HillsDemo::InitApp()
{

	if (!DxAppBase::InitApp())
	{
		return false;
	}

	BuildGeometryBuffers();
	BuildFX();
	BuildVertexLayout();

	return true;

}


bool HillsDemo::OnResizeHandler()
{

	DxAppBase::OnResizeHandler();

	//Window resized, update aspect ratio, recompute the projection matrix

	XMMATRIX p = XMMatrixPerspectiveFovLH(0.25f*MathHelper::Pi, CurAspectRatio(), 1.0f, 1000.0f);
	XMStoreFloat4x4(&mProj, p);
	return true;
}

void HillsDemo::ProcSceneUpdate(float dt)
{

	//Convert spherical coords to radial
	// (radius * sin(mPhi) = hypotenus of triangle in xz plane from mTheta)
	float x = mRadius*sinf(mPhi)*cosf(mTheta);
	float z = mRadius*sinf(mPhi)*sinf(mTheta);
	float y = mRadius*cosf(mPhi);

	//build view matrix
	XMVECTOR pos = XMVectorSet(x, y, z, 1.0f);
	XMVECTOR target = XMVectorZero();
	XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);

	XMMATRIX V = XMMatrixLookAtLH(pos, target, up);
	XMStoreFloat4x4(&mView, V);
}



void HillsDemo::ProcSceneDraw()
{

	//lock manager
	_dxMgr.LockMgr();

	//get current device context

	ID3D11DeviceContext *pContext = _dxMgr.CurrentDeviceContext();

	//Clear back buffer blue
	pContext->ClearRenderTargetView(_dxMgr.CurrentRenderTargetView(), reinterpret_cast<const float*>(&Colors::Blue));

	//clear depth and stencil buffer
	pContext->ClearDepthStencilView(_dxMgr.CurrentDepthStencilView(), D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0f, 0);

	//set the input layout
	pContext->IASetInputLayout(pInputLayout);

	//set topology
	pContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	UINT stride = sizeof(HillsVertex);
	UINT offset = 0;

	//set vertex buffer
	pContext->IASetVertexBuffers(0, 1, &pHillsVertexBuffer, &stride, &offset);

	//set index buffer
	pContext->IASetIndexBuffer(pHillsIndexBuffer, DXGI_FORMAT_R32_UINT, 0);


	//set constants
	XMMATRIX world = XMLoadFloat4x4(&mWorld);
	XMMATRIX view = XMLoadFloat4x4(&mView);
	XMMATRIX proj = XMLoadFloat4x4(&mProj);

	XMMATRIX worldViewProj = world*view*proj;

	//update the constant buffer wvp matrix for vertex shader
	pfxWorldViewProj->SetMatrix(reinterpret_cast<float*>(&worldViewProj));

	//get a handle to the technique description
	D3DX11_TECHNIQUE_DESC techDesc;
	pTech->GetDesc(&techDesc);

	for (UINT p = 0; p < techDesc.Passes; ++p)
	{

		//get pass by index... update constant buffer in video memory, bind shaders to pipeline
		pTech->GetPassByIndex(p)->Apply(0, pContext);

		//Draw the 36 indicies into the vertex buffer to the back buffer
		pContext->DrawIndexed(mGridIndexCount, 0, 0);

	}

	//flip the back buffer and... front buffer?
	HR(_dxMgr.CurrentSwapChain()->Present(0, 0));

	//unlock mgr
	_dxMgr.UnlockMgr();

}


void HillsDemo::HandleMouseDown(WPARAM bState, int x, int y)
{

	mLastMousePos.x = x;
	mLastMousePos.y = y;

	SetCapture(this->handleMainWindow);
}

void HillsDemo::HandleMouseUp(WPARAM bState, int x, int y)
{
	ReleaseCapture();
}


void HillsDemo::HandleMouseMove(WPARAM bState, int x, int y)
{

	if ((bState & MK_LBUTTON) != 0)
	{

		//make each pixel correspond to a quarter of a degree
		float dx = XMConvertToRadians(0.25f*static_cast<float>(x - mLastMousePos.x));
		float dy = XMConvertToRadians(0.25f*static_cast<float>(y - mLastMousePos.y));

		//Update angles based on input to orbit camera around box
		mTheta += dx;
		mPhi += dy;

		//restrict angle mPhi
		mPhi = MathHelper::Clamp(mPhi, 0.1f, MathHelper::Pi - 0.1f);
	}
	else if ((bState & MK_RBUTTON) != 0)
	{

		//and lets change this so we can move in direction we are facing using y axis
		float dx = 0.010f * static_cast<float>(x - mLastMousePos.x);
		float dy = 0.010f * static_cast<float>(y - mLastMousePos.y);

		//update camera radius based on input
		mRadius += dx -dy;

		//restrict the radius
		mRadius = MathHelper::Clamp(mRadius, 50.0f, 500.0f);
	}

	mLastMousePos.x = x;
	mLastMousePos.y = y;

}

float HillsDemo::GetHeight(float x, float z) const
{
	return 0.3f*(z*sinf(0.1f*x) + x*cosf(0.1f*z));
}

void HillsDemo::BuildGeometryBuffers()
{

	EPGeometry::EPGeometryGenerator::MeshData grid;

	HRESULT gridResult = EPGeometry::EPGeometryGenerator::CreateGrid(160.0f, 160.0f, 50, 50, grid);

	mGridIndexCount = grid.Indices.size();
	mGridVertexCount = grid.Vertices.size();

	std::vector<HillsVertex> vertices(mGridVertexCount);

	for (size_t i = 0; i < mGridVertexCount; i++)
	{

		XMFLOAT3 p = grid.Vertices[i].Position;

		p.y = GetHeight(p.x, p.z);

		vertices[i].Pos = p;

		//color vertex based on height
		if (p.y < -10.f)
		{
			//sandy beach
			vertices[i].Color = XMFLOAT4(10.f, 0.96f, 0.62f, 1.0f);
		}
		else if (p.y < 0.5f)
		{
			//light yellow-green
			vertices[i].Color = XMFLOAT4(0.48f, 0.77f, 0.46f, 1.0f);
		}
		else if (p.y < 12.0f)
		{
			//dark yellow-green
			vertices[i].Color = XMFLOAT4(0.1f, 0.48f, 0.19f, 1.0f);
		}
		else if (p.y < 20.f)
		{
			//dark brown
			vertices[i].Color = XMFLOAT4(0.45f, 0.39f, 0.34f, 1.0f);
		}
		else
		{
			//white snow
			vertices[i].Color = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);

		}
	}

	//Create our vertex buffers
	D3D11_BUFFER_DESC vbd, ibd;
	vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vbd.ByteWidth = mGridVertexCount * sizeof(HillsVertex);
	vbd.CPUAccessFlags = vbd.MiscFlags = vbd.StructureByteStride = 0;
	vbd.Usage = D3D11_USAGE_IMMUTABLE;
	
	D3D11_SUBRESOURCE_DATA vinitData, iinitData;
	vinitData.SysMemPitch = vinitData.SysMemSlicePitch = 0;
	vinitData.pSysMem = &vertices[0];

	//Create our index buffer
	ibd.BindFlags = D3D11_BIND_INDEX_BUFFER;
	ibd.Usage = D3D11_USAGE_IMMUTABLE;
	ibd.CPUAccessFlags = ibd.MiscFlags = ibd.StructureByteStride = 0;
	ibd.ByteWidth = mGridIndexCount * sizeof(UINT);

	iinitData.SysMemPitch = iinitData.SysMemSlicePitch = 0;

	//We don't need to change our indices, we can use them directly from grid
	iinitData.pSysMem = &grid.Indices[0];

	//lock mgr
	_dxMgr.LockMgr();

	//create our vertex buffer
	HR(_dxMgr.CurrentDevice()->CreateBuffer(&vbd, &vinitData, &pHillsVertexBuffer));

	//Create our index buffer
	HR(_dxMgr.CurrentDevice()->CreateBuffer(&ibd, &iinitData, &pHillsIndexBuffer));

	//unlock mgr
	_dxMgr.UnlockMgr();
}

void HillsDemo::BuildFX()
{
	std::string targetFile;
#ifdef _DEBUG
	targetFile = "fx/color.cso";
#else
	targetFile = "fx/color.cso";
#endif
	std::ifstream fin(targetFile, std::ios::binary);

	if (fin.fail())
	{
		fin.close();
		DXTrace(__FILE__, (DWORD)__LINE__, 0, 0, true);
		//panic!
		MessageBox(NULL, _T("Couldnt open shader"), _T("error"), MB_ABORTRETRYIGNORE);
	}

	fin.seekg(0, std::ios_base::end);
	int size = (int)fin.tellg();
	fin.seekg(0, std::ios_base::beg);
	std::vector<char> compiledShader(size);

	fin.read(&compiledShader[0], size);
	fin.close();

	_dxMgr.LockMgr();

	HR(D3DX11CreateEffectFromMemory(&compiledShader[0], size, 0, _dxMgr.CurrentDevice(), &pFX));

	_dxMgr.UnlockMgr();

	pTech = pFX->GetTechniqueByName("ColorTech");

	pfxWorldViewProj = pFX->GetVariableByName("gWorldViewProj")->AsMatrix();

}

void HillsDemo::BuildVertexLayout()
{

	//create vertex input layout

	D3D11_INPUT_ELEMENT_DESC zeLayout[] =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 }
	};

	//Create the input layout
	_dxMgr.LockMgr();

	D3DX11_PASS_DESC pDesc;
	pTech->GetPassByIndex(0)->GetDesc(&pDesc);
	HR(_dxMgr.CurrentDevice()->CreateInputLayout(zeLayout, 2, pDesc.pIAInputSignature, pDesc.IAInputSignatureSize, &pInputLayout));

	_dxMgr.UnlockMgr();

}


int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE prevInstance,
	PSTR cmdLine, int showCmd)
{
#ifdef _DEBUG
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

#ifdef _testingMASM


#endif

	HillsDemo theApp(hInstance);

	if (!theApp.InitApp())
	{
		return 0;
	}

	return theApp.Run();

}