
#include "CLPTestMultipleCLP.h"

#include <vector>

template<class T>
bool Compare(T* expectedValues,
             unsigned int expectedNumberOfValues,
             std::vector<T>& values
             )
{
  bool success = values.size() == expectedNumberOfValues;
  if (!success)
    {
    std::cerr << "Wrong number of values, got: " << values.size()
      << " expected: " << expectedNumberOfValues << std::endl;
    return false;
    }

  for (size_t i = 0; i < expectedNumberOfValues; ++i)
    {
    success &= values[i] == expectedValues[i];
    if (!success)
      {
      std::cerr << "Wrong values, got: " << values[i]
        << " expected: " << expectedValues[i] << std::endl;
      }
    }
  return success;
}

//-----------------------------------------------------------------------------
int main(int argc, char * argv[])
{
  PARSE_ARGS;

  const unsigned int numberOfExpectedValues = 6;
  int expectedValues[numberOfExpectedValues] = {1, 2, 3, 4, 5, 6};
  int negativeExpectedValues[numberOfExpectedValues] = {-1, -2, -3, -4, -5, -6};
  std::string stringExpectedValues[numberOfExpectedValues] = {"1", "2", "3", "4", "5", "6"};
  std::string negativeStringExpectedValues[numberOfExpectedValues] = {"-1", "-2", "-3", "-4", "-5", "-6"};

  int floatVectorSizes[3] = {1, 3, 2};
  float floatExpectedValues[3][3] = {
      {1, 0, 0},
      {2, 3, 4},
      {5, 6, 0}
    };
  int negativeFloatVectorSizes[3] = {2, 1, 3};
  float negativeFloatExpectedValues[3][3] = {
      {-1, -2, 0},
      {-3, 0, 0},
      {-4, -5, -6}
    };

  bool success = false;
  if (IntList.size() > 0)
    {
    success = Compare<int>(
      TestCommandLineOverWrite ? negativeExpectedValues : expectedValues,
      numberOfExpectedValues,
      IntList);
    }
  else if (FileList.size() > 0)
    {
    success = Compare<std::string>(
      TestCommandLineOverWrite ? negativeStringExpectedValues : stringExpectedValues,
      numberOfExpectedValues,
      FileList);
    }
  else if (DirList.size() > 0)
    {
    success = Compare<std::string>(
      TestCommandLineOverWrite ? negativeStringExpectedValues : stringExpectedValues,
      numberOfExpectedValues,
      DirList);
    }
  else if (ImageList.size() > 0)
    {
    success = Compare<std::string>(
      TestCommandLineOverWrite ? negativeStringExpectedValues : stringExpectedValues,
      numberOfExpectedValues,
      ImageList);
    }
  else if (GeometryList.size() > 0)
    {
    success = Compare<std::string>(
      TestCommandLineOverWrite ? negativeStringExpectedValues : stringExpectedValues,
      numberOfExpectedValues,
      GeometryList);
    }
  else if (PointList.size() > 0)
    {
    success = PointList.size() == 3;
    if (!success)
      {
      std::cerr << "Wrong number of values, got: " << PointList.size()
        << " expected: " << 3 << std::endl;
      return EXIT_FAILURE;
      }
    for (size_t i = 0; i < 3; ++i)
      {
      success &=
        Compare<float>(
          TestCommandLineOverWrite ? negativeFloatExpectedValues[i] : floatExpectedValues[i],
          TestCommandLineOverWrite ? negativeFloatVectorSizes[i] : floatVectorSizes[i],
          PointList[i]);
      }
    }
  else if (PointFileList.size() > 0)
    {
    success = Compare<std::string>(
      TestCommandLineOverWrite ? negativeStringExpectedValues : stringExpectedValues,
      numberOfExpectedValues,
      PointFileList);
    }
  else if (RegionList.size() > 0)
    {
    success = RegionList.size() == 3;
    if (!success)
      {
      std::cerr << "Wrong number of values, got: " << RegionList.size()
        << " expected: " << 3 << std::endl;
      return EXIT_FAILURE;
      }
    for (size_t i = 0; i < 3; ++i)
      {
      success &=
        Compare<float>(
          TestCommandLineOverWrite ? negativeFloatExpectedValues[i] : floatExpectedValues[i],
          TestCommandLineOverWrite ? negativeFloatVectorSizes[i] : floatVectorSizes[i],
          RegionList[i]);
      }
    }

  if (success)
    {
    return EXIT_SUCCESS;
    }
  std::cerr << "Fail" << std::endl;
  return EXIT_FAILURE;
}
