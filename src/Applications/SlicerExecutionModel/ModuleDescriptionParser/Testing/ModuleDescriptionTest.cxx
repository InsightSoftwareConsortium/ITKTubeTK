
// ModuleDescriptionParser includes
#include "ModuleDescription.h"
#include "ModuleDescriptionTestingMacros.h"

// STD includes
#include <cstdlib>
#include <string>

//---------------------------------------------------------------------------
int TestValues();
int TestReadParameterFileWithMissingValue();
int TestParameterFileWithPointFile();
int TestTargetCallback();

//---------------------------------------------------------------------------
namespace
{
  std::string INPUT_DIR;
}

//---------------------------------------------------------------------------
int ModuleDescriptionTest(int argc, char * argv[])
{
  if (argc < 2)
    {
    std::cout << "Usage: " << argv[0] << " /path/to/inputs" << std::endl;
    return EXIT_FAILURE;
    }

  INPUT_DIR = std::string(argv[1]);

  CHECK_EXIT_SUCCESS(TestValues());
  CHECK_EXIT_SUCCESS(TestReadParameterFileWithMissingValue());
  CHECK_EXIT_SUCCESS(TestParameterFileWithPointFile());
  CHECK_EXIT_SUCCESS(TestTargetCallback());

  return EXIT_SUCCESS;
}

//---------------------------------------------------------------------------
int TestValues()
{
  ModuleLogo logo;
  CHECK_INT(logo.GetWidth(),        0);
  CHECK_INT(logo.GetHeight(),       0);
  CHECK_INT(logo.GetPixelSize(),    0);
  CHECK_INT(logo.GetBufferLength(), 0);
  CHECK_INT(logo.GetOptions(),      0);
  CHECK_STRING(logo.GetLogo(),      "");

  ModuleDescription desc;
  CHECK_STD_STRING(desc.GetCategory(),         "Unspecified");
  CHECK_STD_STRING(desc.GetIndex(),            "65535");
  CHECK_STD_STRING(desc.GetTitle(),            "Unknown");
  CHECK_STD_STRING(desc.GetDescription(),      "No description provided");
  CHECK_STD_STRING(desc.GetVersion(),          "Unspecified");
  CHECK_STD_STRING(desc.GetDocumentationURL(), "");
  CHECK_STD_STRING(desc.GetLicense(),          "");
  CHECK_STD_STRING(desc.GetAcknowledgements(), "Thank you everyone.");
  CHECK_STD_STRING(desc.GetContributor(),      "Anonymous");
  CHECK_STD_STRING(desc.GetType(),             "Unknown");
  CHECK_STD_STRING(desc.GetTarget(),           "");
  CHECK_POINTER(desc.GetLibraryLoader(),       0);
  CHECK_STD_STRING(desc.GetLocation(),         "");
  CHECK_BOOL(desc.HasParameter(""),            false);
  CHECK_BOOL(desc.HasReturnParameters(),       false);
  CHECK_INT(desc.GetParameterGroups().size(),  0);
  CHECK_BOOL(desc.ReadParameterFile(""),       false);
  CHECK_BOOL(desc.WriteParameterFile(""),      false);

  CHECK_POINTER_DIFFERENT(desc.GetProcessInformation(), 0);

  return EXIT_SUCCESS;
}

//---------------------------------------------------------------------------
int TestReadParameterFileWithMissingValue()
{
  std::string input = INPUT_DIR
      + "/parameter-file-with-missing-value-slicer-issue2712.params";

  ModuleParameterGroup group;

  {
    ModuleParameter parameter;
    parameter.SetName("OutputLabel");
    group.AddParameter(parameter);
  }

  {
    ModuleParameter parameter;
    parameter.SetName("SUVMean");
    group.AddParameter(parameter);
  }

  ModuleDescription desc;
  desc.AddParameterGroup(group);

  if (!desc.HasParameter("OutputLabel") || !desc.HasParameter("SUVMean"))
    {
    std::cerr << "Line " << __LINE__
              << " - Parameters are expected."
              << std::endl;
    return EXIT_FAILURE;
    }

  if (!desc.ReadParameterFile(input))
    {
    std::cerr << "Line " << __LINE__
              << " - 'SUVMean' set to a new value. Modification are expected."
              << std::endl;
    return EXIT_FAILURE;
    }

  if (!desc.HasParameter("OutputLabel") || !desc.HasParameter("SUVMean"))
    {
    std::cerr << "Line " << __LINE__
              << " - Problem reading parameters - Parameters are expected."
              << std::endl;
    return EXIT_FAILURE;
    }

  if (desc.GetParameterValue("OutputLabel") != "")
    {
    std::cerr << "Line " << __LINE__
              << " - Problem reading parameters - Value is expected to be empty."
              << std::endl;
    return EXIT_FAILURE;
    }

  if (desc.GetParameterValue("SUVMean") != "2")
    {
    std::cerr << "Line " << __LINE__
              << " - Problem reading parameters - Value is expected to be '2'."
              << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

//---------------------------------------------------------------------------
int TestParameterFileWithPointFile()
{
  std::string input = INPUT_DIR
      + "/parameter-file-with-pointfile-slicer-issue2979.params";

  ModuleParameterGroup group;

  {
    ModuleParameter parameter;
    parameter.SetName("Input Fiducial File");
    parameter.SetValue("input.fcsv");
    parameter.SetTag("pointfile");
    parameter.SetMultiple("false");
    parameter.SetFileExtensionsAsString(".fcsv");
    parameter.SetCoordinateSystem("ras");
    parameter.SetChannel("input");
    group.AddParameter(parameter);
  }

  {
    ModuleParameter parameter;
    parameter.SetName("Output Fiducial File");
    parameter.SetValue("output.fcsv");
    parameter.SetTag("pointfile");
    parameter.SetMultiple("false");
    parameter.SetFileExtensionsAsString(".fcsv");
    parameter.SetCoordinateSystem("lps");
    parameter.SetChannel("output");
    group.AddParameter(parameter);
  }

  ModuleDescription desc;
  desc.AddParameterGroup(group);

  if (!desc.HasParameter("Input Fiducial File") || !desc.HasParameter("Output Fiducial File"))
    {
    std::cerr << "Line " << __LINE__
              << " - Parameters are expected."
              << std::endl;
    return EXIT_FAILURE;
    }
  if (!desc.WriteParameterFile(input, true))
    {
    std::cerr << "Line " << __LINE__
              << " - Unable to write parameter file "
              << input
              << std::endl;
    return EXIT_FAILURE;
    }

  ModuleDescription readDesc;
  if (readDesc.ReadParameterFile(input))
    {
    std::cerr << "Line " << __LINE__
              << " - Unable to read parameter file, something changed"
              << ", but it was reading into an empty description "
              << input
              << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

//---------------------------------------------------------------------------
int TestTargetCallback()
{
  struct Loader
  {
    Loader():Invoked(0){}
    static void loadAndResolve(void* libraryLoader, ModuleDescription& desc)
    {
      Loader* loader = reinterpret_cast<Loader*>(libraryLoader);
      loader->Invoked++;
      desc.SetTarget("loaded-from-callback");
    }
    int Invoked;
  };

  Loader loader;
  ModuleDescription desc;

  // Without lazy loading
  desc.SetTarget("loaded");
  CHECK_STD_STRING(desc.GetTarget(), "loaded");
  CHECK_INT(loader.Invoked, 0);

  // Invalid library loader
  desc.SetTargetCallback(0, Loader::loadAndResolve);
  desc.SetTarget("");
  CHECK_STD_STRING(desc.GetTarget(), "");
  CHECK_INT(loader.Invoked, 0);

  // Invalid callback
  desc.SetTargetCallback(&loader, 0);
  desc.SetTarget("");
  CHECK_STD_STRING(desc.GetTarget(), "");
  CHECK_INT(loader.Invoked, 0);

  // Lazy loading should work
  desc.SetTargetCallback(&loader, Loader::loadAndResolve);
  desc.SetTarget("");
  CHECK_STD_STRING(desc.GetTarget(), "loaded-from-callback");
  CHECK_STD_STRING(desc.GetTarget(), "loaded-from-callback");
  CHECK_INT(loader.Invoked, 1);


  return EXIT_SUCCESS;
}
