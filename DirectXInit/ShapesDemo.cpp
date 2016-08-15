#define DIRECTINPUT_VERSION 0x0800

#include <Windows.h>
#include <d3d11.h>
#include "DxAppBase.h"
#include "d3dUtil.h"
#include "MathHelper.h"
#include <stdint.h>
#include <tchar.h>
#include "d3dx11effect.h"
#include <xnamath.h>
#include "ShapesDemo.h"
#include <iostream>
#include "EPGeometryGenerator.h"
#include "DirectInputHandler.h"

using namespace EPDirectInput;

ShapesDemo::ShapesDemo(HINSTANCE hWnd)
	: DxAppBase(hWnd), pShapesVertexBuffer(NULL), pShapesIndexBuffer(NULL), pFX(NULL), pTech(NULL),
	pfxWorldViewProj(NULL), pInputLayout(NULL), mTheta(1.5f * MathHelper::Pi), mPhi(0.5f * MathHelper::Pi), mRadius(5.0f)
	,mShapesIndexCountSphere(0), mShapesVertexCountSphere(0), mSphereStacks(8), mSphereSlices(8), mSphereRadius(0.2f)

{

	//EPDirectInputMgr mgr(hWnd);
	//mgr.OpenInputMgr();
	//mgr.OpenInputDeviceOfType(EP_INPUT_KBD);

	strMainWindowCaption = _T("Shapes Demo");

	mLastMousePos.x = mLastMousePos.y = 0;

	XMMATRIX I = XMMatrixIdentity();

	//set our world, view, proj matrices to identity
	XMStoreFloat4x4(&mWorld, I);
	XMStoreFloat4x4(&mView, I);
	XMStoreFloat4x4(&mProj, I);
}

ShapesDemo::~ShapesDemo()
{
	ReleaseCOM(pShapesVertexBuffer);
	ReleaseCOM(pShapesIndexBuffer);

	ReleaseCOM(pFX);
	ReleaseCOM(pInputLayout);
}


bool ShapesDemo::InitApp()
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


bool ShapesDemo::OnResizeHandler()
{

	DxAppBase::OnResizeHandler();

	//Window resized, update aspect ratio, recompute the projection matrix

	XMMATRIX p = XMMatrixPerspectiveFovLH(0.25f*MathHelper::Pi, CurAspectRatio(), 1.0f, 1000.0f);
	XMStoreFloat4x4(&mProj, p);
	return true;
}

void ShapesDemo::ProcSceneUpdate(float dt)
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



void ShapesDemo::ProcSceneDraw()
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

	UINT stride = sizeof(ShapesVertex);
	UINT offset = 0;

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

	//draw our sphere, we can use the same input layout, same technique, same constant buffer (but update with sphere world matrix), but need to change vertex/index buffers
	pContext->IASetVertexBuffers(0, 1, &pShapesVertexBuffer, &stride, &offset);
	pContext->IASetIndexBuffer(pShapesIndexBuffer, DXGI_FORMAT_R32_UINT, 0);

	for (UINT i = 0; i < 12; ++i)
	{
		for (UINT p = 0; p < techDesc.Passes; ++p)
		{
			pTech->GetPassByIndex(p)->Apply(0, pContext);

			pContext->DrawIndexed(mShapesIndexCountSphere, 0, 0);
		}
	}

	//flip the back buffer and... front buffer?
	HR(_dxMgr.CurrentSwapChain()->Present(0, 0));

	//unlock mgr
	_dxMgr.UnlockMgr();

}


void ShapesDemo::HandleMouseDown(WPARAM bState, int x, int y)
{

	mLastMousePos.x = x;
	mLastMousePos.y = y;

	SetCapture(this->handleMainWindow);
}

void ShapesDemo::HandleMouseUp(WPARAM bState, int x, int y)
{
	ReleaseCapture();
}


void ShapesDemo::HandleMouseMove(WPARAM bState, int x, int y)
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

float ShapesDemo::GetHeight(float x, float z) const
{
	return 0.3f*(z*sinf(0.1f*x) + x*cosf(0.1f*z));
}

void ShapesDemo::BuildGeometryBuffers()
{

	EPGeometry::EPGeometryGenerator::MeshData sphereTestMesh;

	//quick shapes test
	HRESULT sphereTestResult = EPGeometry::EPGeometryGenerator::CreateGeosphere(5.0f, 5, sphereTestMesh);

	if (sphereTestResult != S_OK)
		_CrtDbgBreak();

	mShapesVertexCountSphere = sphereTestMesh.Vertices.size();
	mShapesIndexCountSphere = sphereTestMesh.Indices.size();

	std::vector<ShapesVertex> sphereTestVertices;

	//end shape test

	//do another pass for the sphere, change color slightly as height changes
	FLOAT colorDelta = 1.0f /( mSphereStacks -2);

	for (UINT i = 0; i < mShapesVertexCountSphere; ++i)
	{
		ShapesVertex v;
		v.Pos = sphereTestMesh.Vertices[i].Position;
		if (v.Pos.y >= 4.5f)
			v.Color = XMFLOAT4(Colors::Green.f);
		else if (v.Pos.y >= 3.0f)
			v.Color = XMFLOAT4(Colors::Yellow.f);
		else if (v.Pos.y >= 1.0f)
			v.Color = XMFLOAT4(Colors::Magenta.f);
		else if (v.Pos.y >= 0.0f)
			v.Color = XMFLOAT4(Colors::Red.f);
		else if (v.Pos.y <= -1.5f && v.Pos.y > -3.25f)
			v.Color = XMFLOAT4(Colors::White.f);
		else
			v.Color = XMFLOAT4(Colors::Green.f);

		sphereTestVertices.push_back(v);
	}

	//end sphere test


	//Create sphere test vertex buffer
	D3D11_BUFFER_DESC stbd;
	stbd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	stbd.ByteWidth = mShapesVertexCountSphere * sizeof(ShapesVertex);
	stbd.CPUAccessFlags = stbd.MiscFlags = stbd.StructureByteStride = 0;
	stbd.Usage = D3D11_USAGE_IMMUTABLE;

	D3D11_SUBRESOURCE_DATA stinitData;
	stinitData.pSysMem = &sphereTestVertices[0];
	stinitData.SysMemPitch = stinitData.SysMemSlicePitch = 0;

	//sphere test index buffer
	D3D11_BUFFER_DESC stid;
	stid.BindFlags = D3D11_BIND_INDEX_BUFFER;
	stid.ByteWidth = mShapesIndexCountSphere * sizeof(UINT);
	stid.CPUAccessFlags = stid.MiscFlags = stid.StructureByteStride = 0;
	stid.Usage = D3D11_USAGE_IMMUTABLE;

	D3D11_SUBRESOURCE_DATA stiinitData;
	stiinitData.pSysMem = &sphereTestMesh.Indices[0];
	stiinitData.SysMemPitch = stiinitData.SysMemSlicePitch = 0;

	//end sphere test stuff

	//lock mgr
	_dxMgr.LockMgr();

	//sphere test start
	
	//create vertex buffer
	HR(_dxMgr.CurrentDevice()->CreateBuffer(&stbd, &stinitData, &pShapesVertexBuffer));

	//index buffer
	HR(_dxMgr.CurrentDevice()->CreateBuffer(&stid, &stiinitData, &pShapesIndexBuffer));


	//sphere test end

	//unlock mgr
	_dxMgr.UnlockMgr();
}

void ShapesDemo::BuildFX()
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

void ShapesDemo::BuildVertexLayout()
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

	ShapesDemo theApp(hInstance);

	if (!theApp.InitApp())
	{
		return 0;
	}

	return theApp.Run();

}
