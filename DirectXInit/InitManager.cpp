#pragma once
#include "stdafx.h"
#include "DirectXInit.h"
#include <xnamath.h>
#include <Windows.h>
#include <d3d11.h>
#include "InitManager.h"
#include <tchar.h>
#include <string>
#include <map>
#include <algorithm>
#include "ScopeLock.h"
#include "d3dUtil.h"


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

Inspired, and in some cases copied, from Frank D. Luna's book "Intro to game programming with dx11"

*/




inline bool COMRelease(IUnknown *targetCOM)
{

	if (targetCOM <= 0)
		return false;
	else
		targetCOM->Release();
	targetCOM = NULL;
	return true;

}


//TODO: Make an inline function for if any calls fail, to automatically release any com interfaces we have aqquired depending on mgrState

DirectXManager::DirectXManager() : 
	curDevice(NULL), curDeviceContext(NULL), mutexHandle(NULL), mgrState(STATE_MGR_FREE), lastValidState(STATE_MGR_FREE), 
	curSwapChain(NULL), bbRenderTargetView(NULL)
{
	mutexHandle = CreateMutex(NULL, false, _T("MGR_MUTEX"));
	if (mutexHandle <= 0)
	{
		mgrState = STATE_INIT_ERROR;
	}
}


HRESULT DirectXManager::CreateDeviceAndContext()
{
	if (mutexHandle <= 0 || mgrState == STATE_INIT_ERROR)
		return -1;

	//Lock is released in destructor when going out of context
	ScopeLock lock(mutexHandle);

	D3D_FEATURE_LEVEL checkForDX11[1];
	checkForDX11[0] = D3D_FEATURE_LEVEL_11_0;
	D3D_FEATURE_LEVEL highestFeatureLevel;

	HRESULT retRes = D3D11CreateDevice(NULL, D3D_DRIVER_TYPE_HARDWARE, NULL,

		//Specify debug flag if in debug mode...
#ifdef _DEBUG
		D3D11_CREATE_DEVICE_DEBUG
#else
		NULL
#endif
		, checkForDX11,
		1,
		D3D11_SDK_VERSION,
		&curDevice,
		&highestFeatureLevel,
		&curDeviceContext
		);

	//Make sure we didn't fail and that the device supports D3D 11
	if (FAILED(retRes))
	{
		mgrState = STATE_INIT_ERROR;
		return retRes;
	}
	else if (highestFeatureLevel < D3D_FEATURE_LEVEL_11_0)
	{
		lastValidState = mgrState = STATE_MGR_INIT;
		mgrState = STATE_INIT_ERROR;
		return retRes;
	}

	lastValidState = mgrState = STATE_MGR_INIT;
	return retRes;

}

HRESULT DirectXManager::Check4xMSAASupport()
{
	if (mutexHandle <= 0 || mgrState != STATE_MGR_INIT)
		return -1;

	//Lock is released in destructor when going out of context
	ScopeLock lock(mutexHandle);

	UINT retQuality = 0;
	HRESULT retRes = curDevice->CheckMultisampleQualityLevels(DXGI_FORMAT_R8G8B8A8_UNORM, 4, &retQuality);

	if (FAILED(retRes))
	{
		mgrState = STATE_INIT_ERROR;
		return retRes;
	}
	
	assert(retQuality > 0);

	m4xMsaaQuality = retQuality;

	return retRes;
}

HRESULT DirectXManager::DescribeSwapChain(bool switchMSAA, bool fullScreen, UINT width, UINT height, HWND nCurWnd)
{
	//Add support for fullscreen later, will need to refactor a bit
	if (mutexHandle <= 0 || mgrState != STATE_MGR_INIT)
	{
		return -1;
	}

	//Lock is released in destructor when going out of context
	ScopeLock lock(mutexHandle);

	if (nCurWnd <= 0)
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	//Zero out the swap chain descriptor structure
	ZeroMemory(&curSwapChainDesc, sizeof(DXGI_SWAP_CHAIN_DESC));

	//Set the flag for msaa, we don't just specify it in a param so we can return later whether
	//enabled or not.
	use4XMSAA = switchMSAA;

	wHeight = height;
	wWidth = width;
	wCurWnd = nCurWnd;

	if (fullScreen)
		wWindowed = 0;
	else
		wWindowed = 1;

	//Fill out of the swap chain back buffer
	curSwapChainDesc.BufferDesc.Width = width;
	curSwapChainDesc.BufferDesc.Height = height;

	//TODO: Add non-default refresh rate check
	curSwapChainDesc.BufferDesc.RefreshRate.Numerator = 60;
	curSwapChainDesc.BufferDesc.RefreshRate.Denominator = 1;

	curSwapChainDesc.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;	//Monitor doesn't output the alpha but we can use it for extra effects later
	curSwapChainDesc.BufferDesc.ScanlineOrdering = DXGI_MODE_SCANLINE_ORDER_UNSPECIFIED;	//Leave it up to the adapter
	curSwapChainDesc.BufferDesc.Scaling = DXGI_MODE_SCALING_UNSPECIFIED;	//Leave it up to the adapter

	if (use4XMSAA)
	{
		curSwapChainDesc.SampleDesc.Count = 4;

		//We got quality from CheckMultisampleQualityLevels earlier..
		curSwapChainDesc.SampleDesc.Quality = m4xMsaaQuality - 1;
	}
	//Else no MSAA
	else
	{
		curSwapChainDesc.SampleDesc.Count = 1;
		curSwapChainDesc.SampleDesc.Quality = 0;
	}

	//We will be rendering to the back buffer...
	curSwapChainDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	//Double buffering for now, make triple buffering an option later TODO
	curSwapChainDesc.BufferCount = 1;
	//Set output window to current window...
	curSwapChainDesc.OutputWindow = wCurWnd;
	curSwapChainDesc.Windowed = wWindowed > 0 ? true : false;

	//Let the adapter choose the most efficient presentation method
	curSwapChainDesc.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;
	curSwapChainDesc.Flags = 0;

	lastValidState = mgrState = STATE_MGR_SWAP_CHAIN_DESCR_CREATED;
	return 0;
}

//Create an instance of the swap chain
HRESULT DirectXManager::CreateSwapChain()
{

	if (mgrState != STATE_MGR_SWAP_CHAIN_DESCR_CREATED || mutexHandle <= 0)
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	//We need to get the instanc eof the IDXGIFactory used to create the device...
	//Time for COM queries...

	IDXGIDevice *dxgiDevice = NULL;

	if (FAILED(curDevice->QueryInterface(__uuidof(IDXGIDevice),
		(void**)&dxgiDevice)))
	{
		DWORD errorWord = GetLastError();
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	//Get the adapter
	IDXGIAdapter *dxgiAdapter = NULL;

	if (FAILED(dxgiDevice->GetParent(__uuidof(IDXGIAdapter),
		(void**)&dxgiAdapter)))
	{
		//Release what we have so far
		COMRelease(dxgiDevice);

		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	//Finally get the factory interface
	IDXGIFactory *dxgiFactory = NULL;

	if (FAILED(dxgiAdapter->GetParent(__uuidof(IDXGIFactory),
		(void**)&dxgiFactory)))
	{
		//Free the other interfaces we have succesfully aqquired
		COMRelease(dxgiDevice);
		COMRelease(dxgiAdapter);

		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	//Now create the swap chain
	IDXGISwapChain *mSwapChain = NULL;
	if (FAILED(dxgiFactory->CreateSwapChain(curDevice, &curSwapChainDesc, &mSwapChain)))
	{
		COMRelease(dxgiDevice);
		COMRelease(dxgiAdapter);
		COMRelease(dxgiFactory);
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	COMRelease(dxgiDevice);
	COMRelease(dxgiAdapter);
	COMRelease(dxgiFactory);

	curSwapChain = mSwapChain;
	lastValidState = mgrState = STATE_MGR_SWAP_CHAIN_CREATED;
	
	return 0;
}

//Create a render target view for the back buffer of the swap chain
HRESULT DirectXManager::CreateRenderTargetView()
{

	if (mgrState != STATE_MGR_SWAP_CHAIN_CREATED || mutexHandle <= 0)
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	ScopeLock lock(mutexHandle);

	//Handle to back buffer
	ID3D11Texture2D *backBuffer;

	//Get a pointer to the swap chain back buffer, for now just double buffering so index 0
	if (FAILED(curSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (void**)&backBuffer)))
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	
	//Create a render target view to the back buffer. We can leave this NULL since we specified the type of the resource earlier.
	if (FAILED(curDevice->CreateRenderTargetView(backBuffer, NULL, &bbRenderTargetView)))
	{
		//Release the handle to the back buffer in either case
		COMRelease(backBuffer);
		mgrState = STATE_INIT_ERROR;
		return -1;
	}
	else
		COMRelease(backBuffer);

	lastValidState = mgrState = STATE_MGR_RENDER_TARGET_VIEW_CREATED;
	return 0;
}

//Create the depth/stencil texture and a view which we can bind to

HRESULT DirectXManager::CreateDepthStencilBufferAndView()
{
	if (mgrState != STATE_MGR_RENDER_TARGET_VIEW_CREATED || mutexHandle <= 0)
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}


	//Zero out the depthStencilDesc struct
	ZeroMemory(&depthStencilDesc, sizeof(D3D11_TEXTURE2D_DESC));

	//Set the width and height based on our stored values
	depthStencilDesc.Width = wWidth;
	depthStencilDesc.Height = wHeight;

	//Mip levels and array size are 1 for depth/stencil buffer
	depthStencilDesc.MipLevels = 1;
	depthStencilDesc.ArraySize = 1;

	//24 bits normalized to [0,1] for depth and 8 bits for -128, 127, for stencil.
	depthStencilDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;

	if (use4XMSAA)
	{
		depthStencilDesc.SampleDesc.Count = 4;
		depthStencilDesc.SampleDesc.Quality = m4xMsaaQuality - 1;
	}
	else
	{
		depthStencilDesc.SampleDesc.Count = 1;
		depthStencilDesc.SampleDesc.Quality = 0;
	}

	depthStencilDesc.Usage = D3D11_USAGE_DEFAULT;
	depthStencilDesc.BindFlags = D3D11_BIND_DEPTH_STENCIL;
	depthStencilDesc.CPUAccessFlags = 0;
	depthStencilDesc.MiscFlags = 0;

	//returns pointer to depth/stencil buffer in mDepthStencilBuffer on success
	if (FAILED(curDevice->CreateTexture2D(&depthStencilDesc, NULL, &mDepthStencilBuffer)))
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	//takes pointer to the resource we want to create a view for, returns pointer to view in mDepthStencilView
	if (FAILED(curDevice->CreateDepthStencilView(mDepthStencilBuffer, NULL, &mDepthStencilView)))
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	lastValidState = mgrState = STATE_MGR_DEPTH_STENCIL_BUFFER_CREATED;
	return 0;
}

//Bind the views to the output merger state
HRESULT DirectXManager::BindBackBufferAndDepthBufferViewsToOutput()
{
	if (mgrState != STATE_MGR_DEPTH_STENCIL_BUFFER_CREATED || mutexHandle <= 0)
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	curDeviceContext->OMSetRenderTargets(1, &bbRenderTargetView, mDepthStencilView);

	lastValidState = mgrState = STATE_MGR_VIEWS_BOUND_TO_OUTPUT;
	return 0;
}


//Set the viewport to the back buffer
//Leave these default 0 for now
HRESULT DirectXManager::SetDefaultViewport(float altX, float altY)
{
	if (mgrState != STATE_MGR_VIEWS_BOUND_TO_OUTPUT || mutexHandle <= 0)
	{
		mgrState = STATE_INIT_ERROR;
		return -1;
	}

	//Zero out the viewport struct...
	ZeroMemory(&curViewport, sizeof(D3D11_VIEWPORT));

	curViewport.TopLeftX = altX;
	curViewport.TopLeftY = altY;
	curViewport.Width = wWidth;
	curViewport.Height = wHeight;

	//Keep the min/max depth at the default of [0,1]
	curViewport.MinDepth = 0.0f;
	curViewport.MaxDepth = 1.0f;

	curDeviceContext->RSSetViewports(1, &curViewport);

	lastValidState = mgrState = STATE_MGR_VIEWPORT_CREATED;

	return 0;

}


//Minimum which needs to be done when a resize occurs.
bool DirectXManager::ResizeHandler()
{

	if (mgrState != STATE_MGR_VIEWPORT_CREATED || mutexHandle <= 0)
	{
		return false;
	}

	ScopeLock lock(mutexHandle);

	assert(curDeviceContext);
	assert(curDevice);
	assert(curSwapChain);

	//Any old views which have a reference to buffers we will destroy need to be released.
	//Stencil/depth buffer needs to go too. We can resize the back buffer, but still need to make a new render target view and new depth/stencil view.

	//Release render target view
	COMRelease(bbRenderTargetView);

	//Release depth/stencil view
	COMRelease(mDepthStencilView);

	//Release depth stencil buffer
	COMRelease(mDepthStencilBuffer);

	//Update swap chain buffer, make render target view again

	//Note that even if we fail at some stage in here, we can still use lastValidState as our reference for what to do in the
	//destructor. If a COM interface has already been released, COMRelease just returns.

	//For now we are doing double buffering, just 1 buffer. Change to new width/height.
	if (FAILED(curSwapChain->ResizeBuffers(1, wWidth, wHeight, DXGI_FORMAT_R8G8B8A8_UNORM, 0)))
	{
		mgrState = STATE_INIT_ERROR;
		return false;
	}

	//Get a pointer to the swap chain's back buffer.
	ID3D11Texture2D *backBuf;

	if (FAILED(curSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (void**)(&backBuf))))
	{
		mgrState = STATE_INIT_ERROR;
		return false;
	}

	//Create the render target view, again we have already defined the swap chain back buffer
	//properties through the swap chain descriptor struct, so we can pass in NULL for the second param.
	if (FAILED(curDevice->CreateRenderTargetView(backBuf, NULL, &bbRenderTargetView)))
	{
		mgrState = STATE_INIT_ERROR;
		return false;
	}

	//Release the interface to the buffer, we don't need anymore...
	COMRelease(backBuf);

	//Now we have to create the depth/stencil buffer and view again
	ZeroMemory(&depthStencilDesc, sizeof(D3D11_TEXTURE2D_DESC));

	depthStencilDesc.Width = wWidth;
	depthStencilDesc.Height = wHeight;
	depthStencilDesc.ArraySize = 1;
	depthStencilDesc.MipLevels = 1;
	depthStencilDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;

	if (use4XMSAA)
	{
		depthStencilDesc.SampleDesc.Count = 4;
		depthStencilDesc.SampleDesc.Quality = m4xMsaaQuality - 1;
	}
	else
	{
		depthStencilDesc.SampleDesc.Count = 1;
		depthStencilDesc.SampleDesc.Quality = 0;
	}

	depthStencilDesc.Usage = D3D11_USAGE_DEFAULT;
	depthStencilDesc.BindFlags = D3D11_BIND_DEPTH_STENCIL;
	depthStencilDesc.CPUAccessFlags = 0;
	depthStencilDesc.MiscFlags = 0;

	//create buffer
	if (FAILED(curDevice->CreateTexture2D(&depthStencilDesc, NULL, &mDepthStencilBuffer)))
	{
		mgrState = STATE_INIT_ERROR;
		return false;
	}

	//create view to depth stencil buffer
	if (FAILED(curDevice->CreateDepthStencilView(mDepthStencilBuffer, NULL, &mDepthStencilView)))
	{
		mgrState = STATE_INIT_ERROR;
		return false;
	}

	//bind depth/stencil view and render target view to the pipeline
	curDeviceContext->OMSetRenderTargets(1, &bbRenderTargetView, mDepthStencilView);

	//Set viewport
	curViewport.TopLeftX = 0;
	curViewport.TopLeftY = 0;
	curViewport.Width = (float)wWidth;
	curViewport.Height = (float)wHeight;
	curViewport.MinDepth = 0.0f;
	curViewport.MaxDepth = 1.0f;

	curDeviceContext->RSSetViewports(1, &curViewport);

	return true;
}


//We need to put releases for com interfaces depending on state when destructor is hit...
DirectXManager::~DirectXManager()
{
	//If mutex handle is not valid, we don't have anything to clean up.
	if (mutexHandle > 0)
	{
		Clean();
		CloseHandle(mutexHandle);
	}
}

//Meh figure out a nicer way to do this, since there is a lot of code reptition. Maybe inline functions for each step.
//Test each of these states too.
void DirectXManager::Clean()
{

	//We don't need to check mutex handle or mgrState, since mutex handle was checked before this was called.

	//lastValidState tracks what we need to clean, so if we encounter an error, we can clean the stuff up to that point.

	ScopeLock lock(mutexHandle);

	if (lastValidState >= STATE_MGR_DEPTH_STENCIL_BUFFER_CREATED)
	{
		COMRelease(bbRenderTargetView);
		COMRelease(mDepthStencilView);
		COMRelease(curSwapChain);
		COMRelease(mDepthStencilBuffer);

		if (curDeviceContext)
		{
			try
			{
				curDeviceContext->ClearState();
			}
			catch (...)
			{

			}
		}

		COMRelease(curDeviceContext);
		COMRelease(curDevice);
	}
	else if (lastValidState >= STATE_MGR_RENDER_TARGET_VIEW_CREATED)
	{
		COMRelease(bbRenderTargetView);
		COMRelease(curSwapChain);
		if (curDeviceContext)
		{
			try
			{
				curDeviceContext->ClearState();
			}
			catch (...)
			{

			}
		}
		COMRelease(curDeviceContext);
		COMRelease(curDevice);
	}
	else if (lastValidState >= STATE_MGR_SWAP_CHAIN_CREATED)
	{
		COMRelease(curSwapChain);
		if (curDeviceContext)
		{
			try
			{
				curDeviceContext->ClearState();
			}
			catch (...)
			{

			}
		}
		COMRelease(curDeviceContext);
		COMRelease(curDevice);
	}
	else if (lastValidState >= STATE_MGR_INIT)
	{
		if (curDeviceContext)
		{
			try
			{
				curDeviceContext->ClearState();
			}
			catch (...)
			{

			}
		}
		COMRelease(curDeviceContext);
		COMRelease(curDevice);
	}

}

//Borrowed from Frank D Luna's excellent DX11 book, added mutex and isValid flag

GameTimer::GameTimer()
	: mSecondsPerCount(0.0), mDeltaTime(-1.0), mBaseTime(0),
	mPausedTime(0), mPrevTime(0), mCurrTime(0), mStopped(false), isValid(false)
{

	__int64 ticksPerSec;
	if (!QueryPerformanceFrequency((LARGE_INTEGER*)&ticksPerSec))
	{
		isValid = false;
		return;
	}
	else
	{
		mSecondsPerCount = 1.0 / (double)ticksPerSec;
	}

	//Try to create the mutex...
	timerMutex = CreateMutex(NULL, false, _T("TIMER_MUTEX"));

	if (timerMutex <= 0)
	{
		isValid = false;
		return;
	}

	isValid = true;
	return;
}


void GameTimer::Tick()
{

	if (!isValid)
	{
		return;
	}

	ScopeLock lock(timerMutex);

	if (mStopped)
	{
		mDeltaTime = 0.0;
		return;
	}

	//Get the time this frame.
	__int64 currTime;

	if (!QueryPerformanceCounter((LARGE_INTEGER*)&currTime))
	{
		mDeltaTime = 0.0;
		isValid = false;
		return;
	}

	mCurrTime = currTime;

	//Diff between this frame and prev
	mDeltaTime = (mCurrTime - mPrevTime)*mSecondsPerCount;

	//For next frame
	mPrevTime = mCurrTime;

	//Force nonnegative.
	if (mDeltaTime < 0.0)
	{
		mDeltaTime = 0.0;
	}

}

void GameTimer::Reset()
{
	if (!isValid)
	{
		return;
	}

	ScopeLock lock(timerMutex);

	__int64 currTime;

	if (!QueryPerformanceCounter((LARGE_INTEGER*)&currTime))
	{
		isValid = false;
		return;
	}

	mBaseTime = currTime;
	mPrevTime = currTime;
	mStopTime = 0;
	mStopped = false;
}

void GameTimer::Stop()
{

	if (!isValid)
		return;

	ScopeLock lock(timerMutex);

	//If stopped, return
	if (!mStopped)
	{
		__int64 currTime;
		if (!QueryPerformanceCounter((LARGE_INTEGER*)&currTime))
		{
			isValid = false;
			return;
		}

		mStopTime = currTime;
		mStopped = true;

	}

}

void GameTimer::Start()
{

	if (!isValid)
	{
		return;
	}

	ScopeLock lock(timerMutex);

	__int64 startTime;

	if (!QueryPerformanceCounter((LARGE_INTEGER*)&startTime))
	{
		isValid = false;
		return;
	}

	//If resuming...
	if (mStopped)
	{

		//accumulate paused time
		mPausedTime += (startTime - mStopTime);

		//current prev time not valid, as it was last updated before paused.
		//Reset prev time to current time.

		mPrevTime = startTime;

		//no longer stopped
		mStopTime = 0;
		mStopped = false;

	}

}


//And finally TotalTime, returns time since Reset was called (not counting pause time)
float GameTimer::TotalTime() const
{

	//we care about logical constness rather than bitwise constness..
	ScopeLock lock(static_cast<HANDLE>(timerMutex));

	//If stopped, do not count time passed since stopped.
	//If we already had a pause, mStopTime - mBaseTime includes paused time, so we subtract paused time from mStopTime

	if (mStopped)
	{
		return (float)(((mStopTime - mPausedTime) - mBaseTime) * mSecondsPerCount);
	}
	else
	{
		return (float)(((mCurrTime - mPausedTime) - mBaseTime) * mSecondsPerCount);
	}
}


float GameTimer::DeltaTime() const
{
	return (float)mDeltaTime;
}
