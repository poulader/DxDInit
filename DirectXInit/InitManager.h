#pragma once


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
*/


#include <Windows.h>
#include <d3d11.h>
#include <map>

using namespace std;

enum CurState
{
	STATE_INIT_ERROR = -1,
	STATE_MGR_FREE = 0,
	STATE_MGR_INIT,
	STATE_MGR_SWAP_CHAIN_DESCR_CREATED,
	STATE_MGR_SWAP_CHAIN_CREATED,
	STATE_MGR_RENDER_TARGET_VIEW_CREATED,
	STATE_MGR_DEPTH_STENCIL_BUFFER_CREATED,
	STATE_MGR_VIEWS_BOUND_TO_OUTPUT,
	STATE_MGR_VIEWPORT_CREATED,
};


inline bool COMRelease(IUnknown *targetCOM);


//Provides access to D3D device and devicecontext, initialization methods

class DirectXManager
{
public:
	DirectXManager();
	virtual ~DirectXManager();
	HRESULT CreateDeviceAndContext();
	HRESULT Check4xMSAASupport();
	HRESULT DescribeSwapChain(bool switchMSAA, bool fullScreen, UINT width, UINT height, HWND nCurWnd);
	HRESULT CreateSwapChain();
	HRESULT CreateRenderTargetView();
	HRESULT CreateDepthStencilBufferAndView();
	HRESULT BindBackBufferAndDepthBufferViewsToOutput();
	//Leave these default 0 for now
	HRESULT SetDefaultViewport(float altX = 0, float altY = 0);

	inline UINT    GetClientHeight() const { return wHeight; };
	inline UINT    GetClientWidth()  const { return wWidth; };
	inline CurState GetCurrentState() const { return mgrState; };

	inline void	   SetClientDimensions(UINT height, UINT width) { wHeight = height; wWidth = width; };


	//The minimum which needs to be done when a resize occurs
	bool ResizeHandler();


	//These are used if a data member needs to be directly accessed by another class, the scopelocks are used in member functions
	inline bool LockMgr() {
		if (mgrState >= STATE_MGR_FREE) { if (WaitForSingleObject(mutexHandle, 2000)) { return false; } else { isLocked = true; return true; } }
		else { return false; }
	}

	inline bool UnlockMgr() {
		if (isLocked && mutexHandle > 0)
		{
			ReleaseMutex(mutexHandle);
			isLocked = false;
			return true;
		}
		else
			return false;
	}

	//Lock should be obtained before calling any of these, and released after.
	inline ID3D11Device *CurrentDevice() { return curDevice; };
	inline ID3D11DeviceContext *CurrentDeviceContext() const { return curDeviceContext; };
	inline IDXGISwapChain *CurrentSwapChain() { return curSwapChain; };
	inline ID3D11RenderTargetView *CurrentRenderTargetView() const { return bbRenderTargetView; };
	inline DXGI_SWAP_CHAIN_DESC& CurrentSwapChainDesc() { return curSwapChainDesc; };
	inline D3D11_TEXTURE2D_DESC& CurrentDepthStencilDesc() { return depthStencilDesc; };
	inline ID3D11Texture2D* CurrentDepthStencilBuffer() { return mDepthStencilBuffer; };
	inline ID3D11DepthStencilView* CurrentDepthStencilView() { return mDepthStencilView; };
	inline D3D11_VIEWPORT& GetCurrentViewPort() { return curViewport; };

private:

	//Called in destructor, lastValidState tracks where we are as far as COM interface reference counts we need to decrement
	void Clean();


	//Current device
	ID3D11Device *curDevice;

	//Current device context
	ID3D11DeviceContext *curDeviceContext;

	//Current swap chain descriptor
	DXGI_SWAP_CHAIN_DESC curSwapChainDesc;

	//Pointer to swap chain interface
	IDXGISwapChain *curSwapChain;

	//Handle to render target view for swap chain back buffer
	ID3D11RenderTargetView *bbRenderTargetView;

	//depth stencil texture descriptor
	D3D11_TEXTURE2D_DESC depthStencilDesc;

	//depth/stencil texture
	ID3D11Texture2D *mDepthStencilBuffer;

	//depth/stencil view to bind to output pipeline
	ID3D11DepthStencilView *mDepthStencilView;

	//Keep a copy of the viewport, by default we will set it to fill the entire backbuffer.
	//We will just use one viewport for now.
	D3D11_VIEWPORT curViewport;

	//State (error, free, init, disposing)
	CurState mgrState;

	//Last valid state (same as mgrState if no error)
	CurState lastValidState;

	//switch for 4x msaa
	bool use4XMSAA;

	//Height, width, windowed
	UINT wHeight, wWidth, wWindowed;

	//Quality returned from CheckmultisampleQualityLevels, if we are using 4x msaa
	UINT m4xMsaaQuality;

	//Output window handle (mostly for swap chain descriptor)
	HWND wCurWnd;

	//Mutex associated with this instance
	HANDLE mutexHandle;

	//for when an owner class needs to directly lock and unlock,
	//scopelock is used for member functions
	bool isLocked;

};


typedef UINT TimerHandle;

//More or less Frank D. Lunas' "Intro to game programming with Dx11" game timer class, added mutex and isValid

class GameTimer
{

public:
	GameTimer();
	GameTimer(bool setThreadAffinity);


	float TotalTime() const;		//In seconds
	float DeltaTime() const;	//In seconds

	void Reset();	// Call before message loop
	void Start();	// Call when unpaused
	void Stop();	// Call when paused
	void Tick();	// Call every frame

	inline bool GetIsValid() const { return isValid; };

private:

	double mSecondsPerCount;
	double mDeltaTime;

	__int64 mBaseTime;
	__int64 mPausedTime;
	__int64 mStopTime;
	__int64 mPrevTime;
	__int64 mCurrTime;

	bool mStopped;
	bool isValid;
	HANDLE timerMutex;

};