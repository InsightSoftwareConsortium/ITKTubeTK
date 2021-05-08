#include "CLPExample1CLP.h"

#include <string>
#include <vector>

//----------------------------------------------------------------------------
template<class T>
bool HasExpectedValue(
  T value, std::vector<T>& expectedValues, bool printError=true)
{
  bool success = false;
  for (size_t i = 0; i < expectedValues.size(); ++i)
    {
    success |= (expectedValues[i] == value);
    }
  if (!success && printError)
    {
    std::cerr << "Did not find value: " << value
      << " inside the expected values [";
    for (size_t i = 0; i < expectedValues.size(); ++i)
      {
      std::cerr << expectedValues[i] << ", ";
      }
    std::cerr << "]" << std::endl;
    }
  return success;
}

//----------------------------------------------------------------------------
template<class T>
bool HasExpectedValue(
  T value, T expectedValue1, bool printError=true)
{
  std::vector<T> expectedValues;
  expectedValues.push_back(expectedValue1);
  return HasExpectedValue(value, expectedValues, printError);
}

template<class T>
bool HasExpectedValue(
  T value, T expectedValue1,  T expectedValue2, bool printError=true)
{
  std::vector<T> expectedValues;
  expectedValues.push_back(expectedValue1);
  expectedValues.push_back(expectedValue2);
  return HasExpectedValue(value, expectedValues, printError);
}

//----------------------------------------------------------------------------
template<class T>
bool VectorHasExpectedValue(
  std::vector<T> values, std::vector<T> expectedValues, bool printError=true)
{
  if (values.size() != expectedValues.size())
    {
    if (printError)
      {
      std::cerr << "Vector not the same size" << std::endl;
      }
    return false;
    }

  bool success = true;
  for (size_t i = 0; i < expectedValues.size(); ++i)
    {
    success &= HasExpectedValue(values[i], expectedValues[i], printError);
    }
  return success;
}

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
int main (int argc, char *argv[])
{
  PARSE_ARGS;

  bool success = true;

  // interpolationType
  success &= HasExpectedValue<std::string>(interpolationType, "sinc");

  // selection
  success &= HasExpectedValue<int>(selection, 4);

  // HistogramBins
  success &= HasExpectedValue<int>(HistogramBins, 30);

  // SpatialSamples
  success &= HasExpectedValue<int>(SpatialSamples, 3324, 2112);

  // InitializeTransform
  success &= HasExpectedValue<bool>(InitializeTransform, false);

  // Iterations
  std::vector<int> expectedIterations;
  expectedIterations.push_back(100);
  expectedIterations.push_back(100);
  expectedIterations.push_back(100);
  expectedIterations.push_back(200);
  success &= VectorHasExpectedValue<int>(Iterations, expectedIterations);

  // LearningRate
  std::vector<float> expectedLearningRate1;
  expectedLearningRate1.push_back(0.002);
  expectedLearningRate1.push_back(0.001);
  expectedLearningRate1.push_back(0.0007);
  expectedLearningRate1.push_back(0.0002);
  bool learningRateSuccess = false;
  learningRateSuccess |= VectorHasExpectedValue<float>(LearningRate, expectedLearningRate1, false);

  std::vector<float> expectedLearningRate2;
  expectedLearningRate2.push_back(0.001);
  expectedLearningRate2.push_back(0.001);
  expectedLearningRate2.push_back(0.0005);
  expectedLearningRate2.push_back(0.0003);
  learningRateSuccess |= VectorHasExpectedValue<float>(LearningRate, expectedLearningRate2, false);
  if (!learningRateSuccess)
    {
    std::cerr << "Error, unexpected LearningRate value, got: [";
    for (size_t i = 0; i < LearningRate.size(); ++i)
      {
      std::cerr << LearningRate[i] << ", ";
      }
    std::cerr << "]" << std::endl;
    }
  success &= learningRateSuccess;

  // TranslationScale
  success &= HasExpectedValue<double>(TranslationScale, 20);

  // DownsamplingFactor
  success &= HasExpectedValue<int>(DownsamplingFactor, 2, 5);

  // MultipleValues
  if (DefaultMultipleValues.size() != 0)
    {
    std::cerr << "Error, unexpected DefaultMultipleValues size, got: "
      << DefaultMultipleValues.size() << " expected 0." << std::endl;
    return EXIT_FAILURE;
    }

  // FixedImage
  success &= HasExpectedValue<std::string>(FixedImage, std::string("Head.mha"), std::string("OtherHead.mha"));

  // MovingImage
  success &= HasExpectedValue<std::string>(MovingImage, std::string("ProgrammingHead.mha"), std::string("OtherProgrammingHead.mha"));

  // OutputImage
  success &= HasExpectedValue<std::string>(OutputImage, std::string("ShrunkHead.mha"), std::string("OtherShrunkHead.mha"));

  if (!success)
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

