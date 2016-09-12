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

#define _testingMASM


int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE prevInstance,
	PSTR cmdLine, int showCmd)
{
#ifdef _DEBUG
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

	ShapesDemo theApp(hInstance);

	if (!theApp.InitApp())
	{
		return 0;
	}

	return theApp.Run();

}