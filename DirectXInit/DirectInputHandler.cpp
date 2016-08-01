#include <Windows.h>
#include <D3D11.h>
#include <dinput.h>
#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <tchar.h>
#include "ScopeLock.h"
#include "DirectInputHandler.h"
#include <assert.h>

using namespace EPDirectInput;

EXTERN_C IMAGE_DOS_HEADER __ImageBase;
#define HINST_THISCOMPONENT ((HINSTANCE)&__ImageBase)

#pragma comment(lib, "dinput8.lib")

EPDirectInputMgr::EPDirectInputMgr(HINSTANCE nhWnd)
	: pDirectInputObject(NULL), mgrState(IMGR_FREE), pdIkbd(NULL), pdImouse(NULL), _lastErrorHR(0), _lastErrorStr(_T("No error")), _isKbdCreated(false), _isMouseCreated(false), hWnd(nhWnd)
{

};

HRESULT EPDirectInputMgr::OpenInputMgr()
{
	return CreateDirectInputInstance(hWnd);
}

 inline HRESULT EPDirectInputMgr::CreateDirectInputInstance(HINSTANCE hWnd)
{
	HRESULT res = DirectInput8Create(hWnd, DIRECTINPUT_VERSION, IID_IDirectInput8, (void**)&pDirectInputObject, NULL);
	if (res != S_OK)
	{
		mgrState = IMGR_DID_INIT_ERROR;
		LogLastError(res, _T("DirectInput8Create."));
	}
	else
	{
		LogLastError(res, std::wstring(_T("DirectInput8Create OK.")));
		mgrState = IMGR_ACTIVE;
	}

	return res;
}

void EPDirectInputMgr::LogLastError(HRESULT code, const std::wstring& str)
{
	_lastErrorHR = code;
	_lastErrorStr = (str);
	_logBuffer.push_back(std::make_pair(code, str));
}

int EPDirectInputMgr::GetMessageLog(std::vector<std::pair<HRESULT, std::wstring>>& buf, uint32_t count, uint32_t index) const
{

	if (index > _logBuffer.size() - 1)
		return -1;
	else if (count + index > _logBuffer.size() - 1)
		return -2;
	
	std::vector<std::pair<HRESULT, std::wstring>>::const_iterator it = _logBuffer.begin() + index;
	buf.resize(count);
	int copied = 0;
	for (; it != _logBuffer.end() && copied < count; ++it, ++copied)
	{
		buf[copied] = *it;
	}
	return copied;
}

HRESULT EPDirectInputMgr::OpenInputDeviceOfType(INPUT_CLASS deviceClass)
{
	return ImplOpenInputDeviceOfType(deviceClass);
}

HRESULT EPDirectInputMgr::ImplOpenInputDeviceOfType(INPUT_CLASS deviceClass)
{
	if (GetMgrState() != IMGR_ACTIVE)
	{
		LogLastError(-1, _T("DirectInputMgr is not active."));
		return -1;
	}
	else if (deviceClass != INPUT_CLASS::EP_INPUT_KBD && deviceClass != INPUT_CLASS::EP_INPUT_MOUSE)
	{
		LogLastError(-2, _T("Invalid INPUT_CLASS argument."));
		return -2;
	}
	else if (deviceClass == EP_INPUT_KBD && _isKbdCreated)
	{
		LogLastError(0, _T("Keyboard is already created."));
		return 0;
	}
	else if (deviceClass == EP_INPUT_MOUSE && _isMouseCreated)
	{
		LogLastError(0, _T("Mouse is already created."));
		return 0;
	}

	HRESULT createResult = 0;

	if (deviceClass == EP_INPUT_KBD)
	{
		pdIkbd = new EPDirectInputDevice(EP_INPUT_KBD);
		assert(pdIkbd != NULL);
		_isKbdCreated = true;

		if (pdIkbd->GetDeviceState() != DID_FREE)
		{
			delete  pdIkbd;
			_isKbdCreated = false;
			LogLastError(-3, _T("Error creating the device."));
			return -3;
		}

		//try to create the device
		createResult = pDirectInputObject->CreateDevice(GUID_SysKeyboard, &pdIkbd->pDIDevice, NULL);

		if (createResult != S_OK)
		{
			delete pdIkbd;
			pdIkbd = NULL;
			_isKbdCreated = false;
			LogLastError(createResult, _T("CreateDevice returned an error."));
		}
		return createResult;
	}
	else
	{
		pdImouse = new EPDirectInputDevice(EP_INPUT_KBD);
		assert(pdImouse != NULL);
		_isMouseCreated = true;

		if (pdImouse->GetDeviceState() != DID_FREE)
		{
			delete  pdImouse;
			pdImouse = NULL;
			_isMouseCreated = false;
			LogLastError(-3, _T("Error creating the device."));
			return -3;
		}

		//try to create the device
		createResult = pDirectInputObject->CreateDevice(GUID_SysMouse, &pdImouse->pDIDevice, NULL);

		if (createResult != S_OK)
		{
			delete pdImouse;
			pdImouse = NULL;
			_isMouseCreated = false;
			LogLastError(createResult, _T("CreateDevice returned an error."));
		}
		return createResult;
	}
}

EPDirectInputDevice::EPDirectInputDevice(INPUT_CLASS deviceClass)
	: dIDGUID(deviceClass == EP_INPUT_KBD ? GUID_SysKeyboard : (deviceClass == EP_INPUT_MOUSE ? GUID_SysMouse : GUID_NULL)), pDIDevice(NULL), dIDState(DID_FREE),
	dInputClass(deviceClass == EP_INPUT_KBD ? EP_INPUT_KBD : (deviceClass == EP_INPUT_MOUSE ? EP_INPUT_MOUSE : EP_INPUT_UNKNOWN))
{
	if (dInputClass != EP_INPUT_KBD && dInputClass != EP_INPUT_MOUSE)
	{
		dIDState = DID_ERROR;
	}
}

EPDirectInputDevice::EPDirectInputDevice(REFGUID deviceGUID)
	: dIDGUID(deviceGUID == GUID_SysKeyboard ? GUID_SysKeyboard : (deviceGUID == GUID_SysMouse ? GUID_SysMouse : GUID_NULL)), pDIDevice(NULL), dIDState(DID_FREE),
	dInputClass(deviceGUID == GUID_SysKeyboard ? EP_INPUT_KBD : (deviceGUID == GUID_SysMouse ? EP_INPUT_MOUSE : EP_INPUT_UNKNOWN))
{
	if (dInputClass != EP_INPUT_KBD && dInputClass != EP_INPUT_MOUSE)
	{
		dIDState = DID_ERROR;
	}
}




EPDirectInputMgr::~EPDirectInputMgr()
{
	//TODO
}

HRESULT EPDirectInputMgr::CloseInputDeviceOfType(INPUT_CLASS deviceClass)
{
	throw std::exception("Not implemented");
}

HRESULT EPDirectInputMgr::CloseInputMgr()
{
	throw std::exception("Not implemented");
}

EPDirectInputDevice::~EPDirectInputDevice()
{
	//TODO

}