#pragma once

#include "stdafx.h"

#include <tchar.h>
#include <string>
#include <Windows.h>
#include "ScopeLock.h"

/*

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

ScopeLock::ScopeLock(HANDLE &m) :
	mutexHandle(m)
{
	if (mutexHandle <= 0)
#ifdef _DEBUG
		_CrtDbgBreak();
#else
		return;
#endif
	else if (WaitForSingleObject(mutexHandle, 2000))
	{
#ifdef _DEBUG
		_CrtDbgBreak();
#else
		return;
#endif
	}
}

//TODO: Add some optional logging
ScopeLock::~ScopeLock()
{
	if (mutexHandle <= 0)
		return;
	else
		ReleaseMutex(mutexHandle);
}