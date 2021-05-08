#ifdef _WIN32
#include "Windows.h"
#endif

#include "Utf8Example1CLP.h"

#include <iostream>
#include <fstream>

#ifdef _WIN32

// Helper function to get Windows version
typedef LONG NTSTATUS, * PNTSTATUS;
#define STATUS_SUCCESS (0x00000000)
typedef NTSTATUS(WINAPI* RtlGetVersionPtr)(PRTL_OSVERSIONINFOW);
RTL_OSVERSIONINFOW GetRealOSVersion()
{
  HMODULE hMod = ::GetModuleHandleW(L"ntdll.dll");
  if (hMod)
  {
    RtlGetVersionPtr fxPtr = (RtlGetVersionPtr)::GetProcAddress(hMod, "RtlGetVersion");
    if (fxPtr != nullptr)
    {
      RTL_OSVERSIONINFOW rovi = { 0 };
      rovi.dwOSVersionInfoSize = sizeof(rovi);
      if (STATUS_SUCCESS == fxPtr(&rovi))
      {
        return rovi;
      }
    }
  }
  RTL_OSVERSIONINFOW rovi = { 0 };
  return rovi;
}

#endif

int main (int argc, char *argv[])
{
  PARSE_ARGS;

  std::cout << "inputStringEnum: " << inputStringEnum << std::endl;
  std::cout << "inputFile: " << inputFile << std::endl;
  std::cout << "outputFile: " << outputFile << std::endl;

  // Test file writing
  std::ofstream outputFileStream;
  outputFileStream.open(outputFile);
  if (!outputFileStream.is_open())
  {
    std::cerr << "Failed to open file for writing: " << outputFile << std::endl;
    return 1;
  }
  outputFileStream << inputStringEnum << std::endl;
  outputFileStream.close();

  // Test file reading
  std::ifstream inputFileStream;
  inputFileStream.open(outputFile);
  if (!inputFileStream.is_open())
  {
    std::cerr << "Failed to open file for reading: " << inputFile << std::endl;
  }
  std::string writtenReadStringEnum;
  getline(inputFileStream, writtenReadStringEnum);
  inputFileStream.close();
  std::cout << "written and read inputStringEnum: " << writtenReadStringEnum << std::endl;

  // Test output string
  // Write out the return parameters in "name = value" form
  std::ofstream rts;
  rts.open(returnParameterFile.c_str());
  rts << "outputString = " << inputStringEnum << std::endl;
  rts.close();

#ifdef _WIN32
  
  // Check current windows version before proceeding
  RTL_OSVERSIONINFOW rovi = GetRealOSVersion();
  std::cout << "Windows version: " << rovi.dwMajorVersion << "." << rovi.dwMinorVersion << " build " << rovi.dwBuildNumber << std::endl;
  if (rovi.dwBuildNumber < 18362)
  {
    std::cout << "This Windows version does not support UTF-8 as active code page in application (minimum build number is 18362). Further testing is skipped." << std::endl;
    return 0;
  }

  // Check that active code page is UTF-8 on Windows
  UINT activeCodePage = GetACP();
  std::cout << "Active code page: " << activeCodePage << std::endl;
  if (activeCodePage != CP_UTF8)
  {
    std::cerr << "Error: active code page is " << activeCodePage << ", expected " << CP_UTF8 << " (UTF-8)" << std::endl;
    return 1;
  }
#endif

  // Check that we can create a file wit utf8 filename using standard file API
  std::string filenameUtf8 = u8"alpha(\u03b1).txt";
  outputFileStream.open(filenameUtf8);
  if (!outputFileStream.is_open())
  {
    std::cerr << "Failed to open file for writing: " << outputFile << std::endl;
    return 1;
  }
  outputFileStream << inputStringEnum << std::endl;
  outputFileStream.close();

  // Test if file is created successfully
  // We use special wide character API on Windows to ensure we properly check for the existence of the file
  // regardless of active code page.
#if _WIN32
  std::wstring filenameW;
  filenameW += L"alpha(";
  filenameW.push_back((wchar_t)(0x03B1));
  filenameW += L").txt";
  FILE* tmp = _wfopen(filenameW.c_str(), L"r");
#else
  FILE* tmp = fopen(filenameUtf8.c_str(), "r");
#endif
  if (!tmp)
  {
    std::cerr << "Expected file does not exist: " << filenameUtf8 << std::endl;
    return 1;
  }
  fclose(tmp);

  // Delete the temporary file
#if _WIN32
  _wunlink(filenameW.c_str());
#else
  unlink(utf8_str.c_str());
#endif

  return 0;
}
