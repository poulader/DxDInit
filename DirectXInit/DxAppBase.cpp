#include "stdafx.h"


/*
Copyright (c) 2016, Eric Pouladian

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


Parts inspired by/adapted from Frank D. Lunas' excellent book "Intro to 3D Game Programming w/ DX11"

*/



#include "DxAppBase.h"
#include <Windows.h>
#include "InitManager.h"
#include <sstream>
#include "ScopeLock.h"
#include <windowsx.h>
#include <assert.h>

using namespace std;

//global instance, more or less singleton model (or half of it anyways)....
DxAppBase *globalDxApp = NULL;

//Alternative window proc

LRESULT CALLBACK AppWndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{

	//We can get some messages before CreateWindow is done, and therefore mainwnd is not valid.
	//So forward hwnd.

	return globalDxApp->WndMsgProc(hwnd, msg, wParam, lParam);
}


DxAppBase::DxAppBase(HINSTANCE wndInstance)
	:
	handleAppInstance(NULL), strMainWindowCaption(_T("DX11 Application")), bEnforce4xMSAA(true),
	handleMainWindow(NULL), bAppPaused(false), bAppMinimized(false), bAppMaximized(false),
	bIsResizing(false), mClientWidth(1920), mClientHeight(1080), bFullScreen(false), maxFramesPerSecond(120)
{

	globalDxApp = this;


	resizeLock = CreateMutex(NULL, false, _T("BASE_LOCK"));

	assert(resizeLock > 0);

}

DxAppBase::~DxAppBase()
{
	//_dxMgr destructor gets called after we go out of scope here
}

HINSTANCE DxAppBase::ProcInstance() const
{
	return handleAppInstance;
}

HWND DxAppBase::ProcWnd() const
{
	return handleMainWindow;
}

float DxAppBase::CurAspectRatio() const
{
	return static_cast<float>(_dxMgr.GetClientWidth()) / _dxMgr.GetClientHeight();
}

int DxAppBase::Run()
{

	MSG curMsg = { NULL };

	//Reset timer...
	_gameTimer.Reset();

	while (curMsg.message != WM_QUIT)
	{

		if (PeekMessage(&curMsg, NULL, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&curMsg);
			DispatchMessage(&curMsg);
		}
		else if (_gameTimer.GetIsValid())
		{

			//Increment timer and get new delta
			_gameTimer.Tick();

			if (!bAppPaused)
			{
				FrameStatUpdate();



				ProcSceneUpdate(_gameTimer.DeltaTime());
				ProcSceneDraw();
			}
			else
			{
				Sleep(100);
			}


		}
		else
		{
			//Something went terribly wrong, game timer is not valid, bail
			//TODO: Add SetLastError of some sort
			return (int)curMsg.wParam;
		}

	}

	return (int)curMsg.wParam;

}

//Initialization code goes here, then overrides can do other stuff
bool DxAppBase::InitApp()
{
	if (!ProcWndInit())
		return FALSE;
	
	if (!D3DInit())
		return FALSE;

	//subclass would call if (!DxAppBase::InitApp()) then do their stuff on success.

	return TRUE;
}

//Create the window
bool DxAppBase::ProcWndInit()
{
	WNDCLASS wc;
	wc.style = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc = AppWndProc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = handleAppInstance;
	wc.hIcon = LoadIcon(0, IDI_APPLICATION);
	wc.hCursor = LoadCursor(0, IDC_ARROW);
	wc.hbrBackground = NULL;
	wc.lpszMenuName = NULL;
	wc.lpszClassName = _T("D3DWndClassName");

	if (!RegisterClass(&wc))
	{
		MessageBox(NULL, _T("RegisterClass Failed."), NULL, MB_OK);
		return false;
	}

	//Get window rectangle dimensions based on requested client area dimensions.
	RECT r = { 0, 0, mClientWidth, mClientHeight };
	if (!AdjustWindowRect(&r, WS_OVERLAPPEDWINDOW, FALSE))
	{
		MessageBox(NULL, _T("AdjustWindowRect failed."), NULL, MB_OK);
		return false;
	}

	int width = r.right - r.left;
	int height = r.bottom - r.top;

	//make-a da window
	handleMainWindow = CreateWindow(_T("D3DWndClassName"), strMainWindowCaption.c_str(), WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, width, height, NULL, NULL, handleAppInstance, NULL);

	if (!handleMainWindow)
	{
		MessageBox(NULL, _T("CreateWindow failed"), NULL, MB_OK);
		return false;
	}

	ShowWindow(handleMainWindow, SW_SHOW);
	UpdateWindow(handleMainWindow);

	return true;
}

bool DxAppBase::D3DInit()
{
	//We need to ignore any resize messages while in this function
	if (WaitForSingleObject(resizeLock, 100))
		return false;

	//If state isn't free, we need to close it first
	if (_dxMgr.GetCurrentState() != STATE_MGR_FREE)
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.CreateDeviceAndContext()))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.Check4xMSAASupport()))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.DescribeSwapChain(bEnforce4xMSAA, bFullScreen, mClientWidth, mClientHeight, handleMainWindow)))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.CreateSwapChain()))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.CreateRenderTargetView()))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.CreateDepthStencilBufferAndView()))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.BindBackBufferAndDepthBufferViewsToOutput()))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (FAILED(_dxMgr.SetDefaultViewport()))
	{
		ReleaseMutex(resizeLock);
		return false;
	}

	if (!ReleaseMutex(resizeLock))
		return false;
	else
		return true;
}



//Window resize handler - have this call something in dxmgr so i dont have to mess around with
//poking at private data members and locking
bool DxAppBase::OnResizeHandler()
{

	//We could potentially get resize messages from the message pump before we have finished intializing everything.
	//However, we will make access to this interface BEFORE the intialization is done, single threaded only.

	//Therefore it will either be free, or completely initialized.
	return _dxMgr.ResizeHandler();

}




//Windows message pump dispatcher/handler, mostly from Frank Luna's code with a few modifications
LRESULT DxAppBase::WndMsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{

	switch (msg)
	{

	//Window activated or deactivated
	case WM_ACTIVATE:
	{
		if (LOWORD(wParam) == WA_INACTIVE)
		{
			bAppPaused = true;
			_gameTimer.Stop();
		}
		else
		{
			bAppPaused = false;
			_gameTimer.Start();
		}
		return 0;

	}

	//Size of window changed
	case WM_SIZE:
	{
		DWORD waitResult = 0;
		if ( (waitResult = WaitForSingleObject(resizeLock, 50)) == WAIT_TIMEOUT)
		{
			return -WAIT_TIMEOUT;
		}
		else if (waitResult != 0)
		{
			return (LONG_PTR)-waitResult;
		}

		mClientWidth = LOWORD(lParam);
		mClientHeight = HIWORD(lParam);

		//At this point we will either have our dxmgr in the FREE state, or it will be completely initalized.
		_dxMgr.LockMgr();

		//Update our dxMgr's height and width
		_dxMgr.SetClientDimensions(mClientHeight, mClientWidth);

		if (wParam == SIZE_MINIMIZED)
		{

			bAppPaused = true;
			bAppMinimized = true;
			bAppMaximized = false;
		}
		else if (wParam == SIZE_MAXIMIZED)
		{ 
			bAppPaused = false;
			bAppMinimized = false;
			bAppMaximized = true;
			_dxMgr.UnlockMgr();
			OnResizeHandler();
		}
		else if (wParam == SIZE_RESTORED)
		{

			if (bAppMinimized)
			{
				//restore from minimized
				bAppPaused = false;
				bAppMinimized = false;
				_dxMgr.UnlockMgr();
				OnResizeHandler();
			}
			else if (bAppMaximized)
			{
				bAppPaused = false;
				bAppMaximized = false;
				_dxMgr.UnlockMgr();
				OnResizeHandler();
			}
			else if (bIsResizing)
			{
				//Still dragging size bar around...



			}
			else
			{
				_dxMgr.UnlockMgr();
				OnResizeHandler();
			}
		}

		return 0;

	}

	//WM_ENTERSIZEMOVE - user grabs resize bar
	case WM_ENTERSIZEMOVE:
	{
		bAppPaused = true;
		bIsResizing = true;
		_gameTimer.Stop();
		return 0;
	}

	//WM_EXITSIZEMOVE - user releases resize bar
	case WM_EXITSIZEMOVE:
	{
		bAppPaused = false;
		bIsResizing = false;
		_gameTimer.Start();
		_dxMgr.UnlockMgr();
		_dxMgr.ResizeHandler();
		return 0;

	}

	case WM_DESTROY:
	{
		PostQuitMessage(0);
		return 0;
	}

	case WM_MENUCHAR:
	{
		return MAKELRESULT(0, MNC_CLOSE);
	}


	//Dont let them make the window too small
	case WM_GETMINMAXINFO:
	{
		((MINMAXINFO*)lParam)->ptMinTrackSize.x = 200;
		((MINMAXINFO*)lParam)->ptMinTrackSize.y = 200;
		return 0;
	}

	case WM_LBUTTONDOWN:
	case WM_MBUTTONDOWN:
	case WM_RBUTTONDOWN:
		HandleMouseDown(wParam, GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
		return 0;

	case WM_LBUTTONUP:
	case WM_MBUTTONUP:
	case WM_RBUTTONUP:
		HandleMouseUp(wParam, GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
		return 0;
	case WM_MOUSEMOVE:
		HandleMouseMove(wParam, GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
		return 0;

	}

	return DefWindowProc(hwnd, msg, wParam, lParam);

}



//Direct from Frank Luna's code w/ some additions to account for unicode/ansi
//Added max FPS enforcement, demo was running at 4k frames per second
void DxAppBase::FrameStatUpdate()
{
	// Code computes the average frames per second, and also the 
	// average time it takes to render one frame.  These stats 
	// are appended to the window caption bar.


	static int nFrames = 0;
	static float timeElapsed = 0;

	//min time between frames in seconds (will be fraction of a second)
	static float frameInterval = (1.0f) / maxFramesPerSecond;
	nFrames++;

	//Check to see if our current timer delta is less than the min time between frames
	if (_gameTimer.DeltaTime() < frameInterval)
	{
		//If so, sleep for the difference
		Sleep((frameInterval - _gameTimer.DeltaTime()) * 1000);
		//Update timer delta
		_gameTimer.Tick();
	}

	// Compute averages over one second period.
	if ((_gameTimer.TotalTime() - timeElapsed) >= 1.0f)
	{
		float fps = (float)nFrames; // fps = frameCnt / 1
		float mspf = 1000.0f / fps;

		std::wostringstream outs;
		outs.precision(6);
		outs << strMainWindowCaption << _T("    ")
			<< _T("FPS: ") << fps << _T("    ")
			<< _T("Frame Time: ") << mspf << _T(" (ms)");
		SetWindowText(handleMainWindow, outs.str().c_str());

		// Reset for next average.
		nFrames = 0;
		timeElapsed += 1.0f;
	}
}