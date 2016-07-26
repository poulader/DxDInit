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


Parts inspired by/adapted from Frank D. Lunas' excellent book "Intro to 3D Game Programming w/ DX11"

*/

#include "InitManager.h"
#include <tchar.h>
#include <string>
#include <Windows.h>


//provides abstract base class with window and d3d init stuff taken care of.

class DxAppBase
{
public:

	DxAppBase(HINSTANCE wndInstance);
	virtual ~DxAppBase();

	HINSTANCE ProcInstance() const;
	HWND	  ProcWnd()		 const;
	float	  CurAspectRatio() const;

	int		  Run();


	virtual bool InitApp();
	virtual bool OnResizeHandler();
	virtual void ProcSceneUpdate(float _dt) = 0;
	virtual void ProcSceneDraw() = 0;

	virtual LRESULT WndMsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);

	//Overrides for mouse input
	virtual void HandleMouseDown(WPARAM bState, int x, int y) { }
	virtual void HandleMouseUp(WPARAM bState, int x, int y) { }
	virtual void HandleMouseMove(WPARAM bState, int x, int y) { }

	inline void	   SetMaxFramesPerSecond(float nMaxFPS) { if (nMaxFPS <= 30 || nMaxFPS > 200) return; maxFramesPerSecond = nMaxFPS; };


protected:

	bool ProcWndInit();
	bool D3DInit();
	void FrameStatUpdate();

protected:


	HINSTANCE handleAppInstance;
	HWND	  handleMainWindow;
	bool	  bAppPaused;
	bool	  bAppMinimized;
	bool	  bAppMaximized;
	bool      bIsResizing;
	bool      bEnforce4xMSAA;
	bool	  bFullScreen;

	HANDLE    resizeLock;


	DirectXManager _dxMgr;
	GameTimer	   _gameTimer;

	int mClientWidth;
	int mClientHeight;

	//default to 120
	float maxFramesPerSecond;

#ifdef UNICODE
	std::wstring strMainWindowCaption;
#else
	std::string  strMainWindowCaption;
#endif

};