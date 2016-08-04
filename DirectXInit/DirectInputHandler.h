#pragma once

#include <Windows.h>
#include <D3D11.h>
#include <dinput.h>
#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <tchar.h>
#include "ScopeLock.h"
#include <string>

namespace EPDirectInput
{

	using namespace std;

	//Typedefs
	typedef enum INPUT_CLASS
	{
		EP_INPUT_UNKNOWN = -1,
		EP_INPUT_MOUSE = 0,
		EP_INPUT_KBD

	} INPUT_CLASS;

	typedef enum INPUT_MGR_STATE : int32_t
	{
		IMGR_DID_INIT_ERROR = -2,
		IMGR_DI_INIT_ERROR = -1,
		IMGR_FREE = 0,
		IMGR_ACTIVE
	} INPUT_MGR_STATE;

	typedef enum DIRECTID_STATE
	{
		DID_ERROR = -1,
		DID_FREE = 0,
		DID_ACTIVE
	} DIRECTID_STATE;


	//I don't want these to be copyable... need to come up with a safe way to store them
	//And I don't want a default constructor, once this thing is initialized, its type must be set.
	class EPDirectInputDevice : private Uncopyable
	{
	public:

		//Inputmgr can access the private fields of an input device ofc
		friend class EPDirectInputMgr;

		EPDirectInputDevice(INPUT_CLASS deviceClass);
		EPDirectInputDevice(REFGUID deviceGUID);

		virtual ~EPDirectInputDevice();

		DIRECTID_STATE GetDeviceState() const { return dIDState; };
		const GUID &    GetDeviceGUID() const { return dIDGUID; };
		INPUT_CLASS    GetDeviceInputClass() const { return dInputClass; };

	private:

		//GUID for the device, since we cannot copy these objects, and the type must be declared
		//at initialization, we can make this const
		const GUID dIDGUID;

		//Handle to DirectInputDevice
		IDirectInputDevice8 *pDIDevice;

		const INPUT_CLASS dInputClass;

		//state
		DIRECTID_STATE dIDState;
	};

	//Anything that holds a pointer to a DirectInput interface must not be copyable
	//I would like to abstract this further so I can use a factory method to return either a directinput mgr
	//or a message pump input mgr which both implement the same itnerface, but this is a learning project and that is getting way off track
	class EPDirectInputMgr : private Uncopyable
	{
	public:
		//constructors + destructors
		EPDirectInputMgr(HINSTANCE nhWnd);
		virtual ~EPDirectInputMgr();

		//Give user a chance to close before destructor in case of exceptions
		virtual HRESULT CloseInputMgr();

		//provide a function to open input mgr to avoid any excpetions in constructor
		virtual HRESULT OpenInputMgr();

		//get state
		inline INPUT_MGR_STATE GetMgrState() const { return mgrState; };

		//logging
		inline void   GetLastError(HRESULT& code, std::wstring& str) const { code = _lastErrorHR; str = _lastErrorStr; };
		int GetMessageLog(std::vector<std::pair<HRESULT, std::wstring>>& buf, uint32_t count, uint32_t index) const;
		inline size_t GetMessageLogSize() const { return _logBuffer.size(); };

		//todo: Add something that returns an ID per instance of EPDirectInputDevice so we can have more than 1 instance of a device class

		//init input device of type INPUT_CLASS specified
		virtual HRESULT OpenInputDeviceOfType(INPUT_CLASS);

		//give user a chance to close before destructor in case an exception occurs, they will at least have the
		//opportunity to check for it. If they don't call this, be it on their own head.
		virtual HRESULT CloseInputDeviceOfType(INPUT_CLASS);


	private:

		HRESULT ImplOpenInputDeviceOfType(INPUT_CLASS);

		inline HRESULT CreateDirectInputInstance(HINSTANCE hWnd);

		void		LogLastError(HRESULT code, const std::wstring& str);

		//Handle to DirectInput 
		IDirectInput8 *pDirectInputObject;

		//Handle to top level window of application
		const HINSTANCE hWnd;

		//current state
		INPUT_MGR_STATE mgrState;

		//screw it, this is a learning project, for now just have two instead of going full generic
		EPDirectInputDevice *pdIkbd, *pdImouse;
		bool _isKbdCreated, _isMouseCreated;

		//last error, if any
		std::wstring _lastErrorStr;
		HRESULT		 _lastErrorHR;

		//buffer for errors
		std::vector<std::pair<HRESULT, std::wstring>> _logBuffer;

	};

}
