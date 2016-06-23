//#include "stdafx.h"
//#include "DxAppBase.h"
//#include <Windows.h>
//#include <assert.h>
//#include "d3dUtil.h"
//
//#ifdef _DEBUG
//#include <crtdbg.h>
//#endif
//
////Demo showing the use of the framework in a derived class
//
//class TestDxInit : public DxAppBase
//{
//public:
//	TestDxInit(HINSTANCE hInstance);
//	~TestDxInit();
//
//	//Override four of the virtual methods
//	bool InitApp();
//	bool OnResizeHandler();
//
//	//Required abstract methods
//	void ProcSceneUpdate(float _dt);
//	void ProcSceneDraw();
//
//};
//
//int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE prevInstance,
//	PSTR cmdLine, int showCmd)
//{
//#ifdef _DEBUG
//	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
//#endif
//
//	TestDxInit theApp(hInstance);
//
//	if (!theApp.InitApp())
//	{
//		return 0;
//	}
//
//	return theApp.Run();
//
//}
//
//TestDxInit::TestDxInit(HINSTANCE hInstance)
//	: DxAppBase(hInstance)
//{
//
//}
//
//
//TestDxInit::~TestDxInit()
//{
//
//}
//
//bool TestDxInit::InitApp()
//{
//	if (!DxAppBase::InitApp())
//		return false;
//
//	return true;
//}
//
//bool TestDxInit::OnResizeHandler()
//{
//	return DxAppBase::OnResizeHandler();
//}
//
//void TestDxInit::ProcSceneUpdate(float _dt)
//{
//
//}
//
//void TestDxInit::ProcSceneDraw()
//{
//	_dxMgr.LockMgr();
//	if (_dxMgr.GetCurrentState() != STATE_MGR_VIEWPORT_CREATED)
//	{
//		_dxMgr.UnlockMgr();
//		return;
//	}
//
//	assert(_dxMgr.CurrentDeviceContext());
//	assert(_dxMgr.CurrentSwapChain());
//
//	//Clear back buffer blue.
//
//	_dxMgr.CurrentDeviceContext()->ClearRenderTargetView(_dxMgr.CurrentRenderTargetView(), (const float*)&Colors::Blue);
//
//	//clear depth buffer to 1.0f and stencil buffer to 0.
//	_dxMgr.CurrentDeviceContext()->ClearDepthStencilView(_dxMgr.CurrentDepthStencilView(), D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0f, 0);
//
//	//Present back buffer to screen
//	_dxMgr.CurrentSwapChain()->Present(0, 0);
//	_dxMgr.UnlockMgr();
//	return;
//}