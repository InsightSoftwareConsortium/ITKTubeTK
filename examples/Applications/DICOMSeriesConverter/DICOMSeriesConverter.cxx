/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include <fstream>
#include <iostream>

#include "tubeMessage.h"

// #include <gdcmCommon.h>
// #include <gdcmDictEntry.h>
#include <gdcmFile.h>
#include <gdcmFileHelper.h>
#include <gdcmSerieHelper.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageSeriesReader.h>
#include <itkMetaDataObject.h>
#include <itksys/SystemTools.hxx>

int
main(int argc, char * argv[])
{

  if (argc < 3)
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " DicomDirectory  outputFileName" << std::endl;
    return EXIT_FAILURE;
  }

  typedef signed short PixelType;
  const unsigned int   Dimension = 3;

  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef itk::ImageFileReader<ImageType>   Reader2DType;
  typedef itk::ImageSeriesReader<ImageType> ReaderType;
  typedef itk::GDCMImageIO                  ImageIOType;

  typedef std::vector<std::string>        SeriesIdContainer;
  typedef std::vector<std::string>        FileNamesContainer;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;

  ReaderType::Pointer  reader;
  ImageIOType::Pointer dicomIO;

  Reader2DType::Pointer reader2D;

  WriterType::Pointer writer;

  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->SetDirectory(argv[1]);

  unsigned int groupID[100];
  unsigned int elementID[100];
  unsigned int count = 0;
  //  InstitutionCodeSequence
  groupID[count] = 0x0008;
  elementID[count++] = 0x0082;
  //  ReferringPhysiciansName
  groupID[count] = 0x0008;
  elementID[count++] = 0x0090;
  //  ReferringPhysiciansAddress
  groupID[count] = 0x0008;
  elementID[count++] = 0x0092;
  //  ReferringPhysiciansTelephoneNumbers
  groupID[count] = 0x0008;
  elementID[count++] = 0x0094;
  //  ReferringPhysiciansIdentificationSequence
  groupID[count] = 0x0008;
  elementID[count++] = 0x0096;
  //  ResponsibleOrganization
  groupID[count] = 0x0008;
  elementID[count++] = 0x0116;
  //  StationName
  groupID[count] = 0x0008;
  elementID[count++] = 0x1010;
  //  InstitutionalDepartmentName
  groupID[count] = 0x0008;
  elementID[count++] = 0x1040;
  //  PhysiciansOfRecord
  groupID[count] = 0x0008;
  elementID[count++] = 0x1048;
  //  PhysiciansOfRecordIdentificationSequence
  groupID[count] = 0x0008;
  elementID[count++] = 0x1049;
  //  PerformingPhysiciansName
  groupID[count] = 0x0008;
  elementID[count++] = 0x1050;
  //  PerformingPhysicianIdentificationSequence
  groupID[count] = 0x0008;
  elementID[count++] = 0x1052;
  //  NameOfPhysiciansReadingStudy
  groupID[count] = 0x0008;
  elementID[count++] = 0x1060;
  //  NameOfPhysiciansReadingStudyIdentificationSequence
  groupID[count] = 0x0008;
  elementID[count++] = 0x1062;
  //  OperatorsName
  groupID[count] = 0x0008;
  elementID[count++] = 0x1070;
  //  OperatorIdentificationSequence
  groupID[count] = 0x0008;
  elementID[count++] = 0x1072;

  //  IssuerOfPatientID
  groupID[count] = 0x0010;
  elementID[count++] = 0x0021;
  //  PatientsInsurancePlanCodeSequence
  groupID[count] = 0x0010;
  elementID[count++] = 0x0050;
  //  OtherPatientsIDs
  groupID[count] = 0x0010;
  elementID[count++] = 0x1000;
  //  OtherPatientNames
  groupID[count] = 0x0010;
  elementID[count++] = 0x1001;
  //  PatientsBirthName
  groupID[count] = 0x0010;
  elementID[count++] = 0x1005;
  //  PatientsAddress
  groupID[count] = 0x0010;
  elementID[count++] = 0x1040;
  //  PatientsMothersBirthName
  groupID[count] = 0x0010;
  elementID[count++] = 0x1060;
  //  MilitaryRank
  groupID[count] = 0x0010;
  elementID[count++] = 0x1080;
  //  MedicalRecordLocator
  groupID[count] = 0x0010;
  elementID[count++] = 0x1090;
  //  PatientsTelephoneNumbers
  groupID[count] = 0x0010;
  elementID[count++] = 0x2154;

  //  DeviceSerialNumber
  groupID[count] = 0x0018;
  elementID[count++] = 0x1000;

  // RequestingPhysicianIdentificationSequence
  groupID[count] = 0x0032;
  elementID[count++] = 0x1031;
  // RequestingPhysician
  groupID[count] = 0x0032;
  elementID[count++] = 0x1032;

  // AdmissionID
  groupID[count] = 0x0038;
  elementID[count++] = 0x0010;
  // IssuerOfAdmissionID
  groupID[count] = 0x0038;
  elementID[count++] = 0x0011;
  // PatientsInstitutionResidence
  groupID[count] = 0x0038;
  elementID[count++] = 0x0400;

  // ScheduledPerformingPhysiciansName
  groupID[count] = 0x0040;
  elementID[count++] = 0x0006;
  // ScheduledPerformingPhysiciansIdentificationSequence
  groupID[count] = 0x0040;
  elementID[count++] = 0x000B;
  // PerformedLocation
  groupID[count] = 0x0040;
  elementID[count++] = 0x0243;
  // NamesOfIntendedRecipientsOfResults
  groupID[count] = 0x0040;
  elementID[count++] = 0x1010;
  // IntendedRecipientsOfResultsIdentificationSequence
  groupID[count] = 0x0040;
  elementID[count++] = 0x1011;
  // PersonIdentificationCodeSequence
  groupID[count] = 0x0040;
  elementID[count++] = 0x1101;
  // PersonAddress
  groupID[count] = 0x0040;
  elementID[count++] = 0x1102;
  // PersonTelephoneNumbers
  groupID[count] = 0x0040;
  elementID[count++] = 0x1103;
  // OrderEnteredBy
  groupID[count] = 0x0040;
  elementID[count++] = 0x2008;
  // OrderEnterersLocation
  groupID[count] = 0x0040;
  elementID[count++] = 0x2009;
  // OrderCallbackPhoneNumber
  groupID[count] = 0x0040;
  elementID[count++] = 0x2010;
  // HumanPerformersOrganization
  groupID[count] = 0x0040;
  elementID[count++] = 0x4036;
  // HumanPerformersName
  groupID[count] = 0x0040;
  elementID[count++] = 0x4037;
  // VerifyingObserverName
  groupID[count] = 0x0040;
  elementID[count++] = 0xA075;
  // PersonName
  groupID[count] = 0x0040;
  elementID[count++] = 0xA123;

  // PhysicianApprovingInterpretation
  groupID[count] = 0x4008;
  elementID[count++] = 0x0114;

  int numID = count;

  std::string privateDir = argv[2];
  privateDir += "-PrivateDicom";
  itksys::SystemTools::MakeDirectory(privateDir.c_str());

  std::string publicDir = argv[2];
  publicDir += "-PublicDicom";
  itksys::SystemTools::MakeDirectory(publicDir.c_str());

  std::string metaDir = argv[2];
  metaDir += "-MetaImage";
  itksys::SystemTools::MakeDirectory(metaDir.c_str());

  try
  {
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    FileNamesContainer        fileNames;

    std::string                       seriesIdentifier;
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while (seriesItr != seriesEnd)
    {
      seriesIdentifier = *seriesItr;

      dicomIO = ImageIOType::New();

      bool imageIs3D;
      fileNames = nameGenerator->GetFileNames(seriesIdentifier);
      if (fileNames.size() < 2)
      {
        imageIs3D = false;

        reader2D = Reader2DType::New();
        reader2D->SetFileName(fileNames.begin()->c_str());

        try
        {
          reader2D->Update();
        }
        catch (itk::ExceptionObject & ex)
        {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }
      else
      {
        imageIs3D = true;

        reader = ReaderType::New();
        reader->SetImageIO(dicomIO);

        reader->SetFileNames(fileNames);

        try
        {
          reader->Update();
        }
        catch (itk::ExceptionObject & ex)
        {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }
      }

      gdcm::File * file = (*(nameGenerator->GetSeriesHelper()->GetCoherentFileList(seriesIdentifier)))[0];

      std::string modality = file->GetEntryValue(0x0008, 0x0060);

      std::string seriesNum = file->GetEntryValue(0x0020, 0x0011);

      std::string sequenceName;
      if (modality.c_str() == "CT")
      {
        sequenceName = file->GetEntryValue(0x0018, 0x0050);
        if (sequenceName == gdcm::GDCM_UNFOUND)
        {
          sequenceName = "";
        }
        else
        {
          double sliceThickness = atof(sequenceName.c_str());
          char   st[80];
          sprintf(st, "%0.2f", sliceThickness);
          sequenceName = st;
        }
      }
      else
      {
        sequenceName = file->GetEntryValue(0x0018, 0x0024);
        if (sequenceName == gdcm::GDCM_UNFOUND)
        {
          sequenceName = "";
        }
      }

      std::string protocolName = file->GetEntryValue(0x0018, 0x1030);
      if (protocolName == gdcm::GDCM_UNFOUND)
      {
        protocolName = "";
      }

      char coord[80];
      if (imageIs3D)
      {
        std::string spacing = file->GetEntryValue(0x0028, 0x0030);
        int         split = spacing.find_first_of("\\");
        int         len = spacing.size() - split - 1;
        double      xSpacing = atof(spacing.substr(0, split).c_str());
        double      ySpacing = atof(spacing.substr(split + 1, len).c_str());
        std::string pos = file->GetEntryValue(0x0020, 0x0032);
        int         splitX = pos.find_first_of("\\");
        int         splitY = pos.find_first_of("\\", splitX + 1);
        int         lenY = splitY - splitX - 1;
        int         lenZ = pos.size() - splitY - 1;
        double      xPos = atof(pos.substr(0, splitY).c_str());
        double      yPos = atof(pos.substr(splitX + 1, lenY).c_str());
        double      zPos = atof(pos.substr(splitY + 1, lenZ).c_str());
        file = (*(nameGenerator->GetSeriesHelper()->GetCoherentFileList(seriesIdentifier)))[1];
        pos = file->GetEntryValue(0x0020, 0x0032);
        splitX = pos.find_first_of("\\");
        splitY = pos.find_first_of("\\", splitX + 1);
        lenY = splitY - splitX - 1;
        lenZ = pos.size() - splitY - 1;
        xPos = (xPos - atof(pos.substr(0, splitY).c_str()));
        yPos = (yPos - atof(pos.substr(splitX + 1, lenY).c_str()));
        zPos = (zPos - atof(pos.substr(splitY + 1, lenZ).c_str()));
        double zSpacing = sqrt(xPos * xPos + yPos * yPos + zPos * zPos);
        sprintf(coord, "%0.2fx%0.2fx%0.2f", xSpacing, ySpacing, zSpacing);
      }
      else
      {
        sprintf(coord, "");
      }

      for (int i = 0; i < seriesNum.size(); i++)
      {
        while (i < seriesNum.size() &&
               !((seriesNum[i] >= 'a' && seriesNum[i] <= 'z') || (seriesNum[i] >= '0' && seriesNum[i] <= '9') ||
                 (seriesNum[i] >= 'A' && seriesNum[i] <= 'Z')))
        {
          seriesNum = seriesNum.erase(i, 1).c_str();
        }
      }

      for (int i = 0; i < sequenceName.size(); i++)
      {
        while (i < sequenceName.size() && !((sequenceName[i] >= 'a' && sequenceName[i] <= 'z') ||
                                            (sequenceName[i] >= '0' && sequenceName[i] <= '9') ||
                                            (sequenceName[i] >= 'A' && sequenceName[i] <= 'Z')))
        {
          sequenceName = sequenceName.erase(i, 1).c_str();
        }
      }

      for (int i = 0; i < protocolName.size(); i++)
      {
        while (i < protocolName.size() && !((protocolName[i] >= 'a' && protocolName[i] <= 'z') ||
                                            (protocolName[i] >= '0' && protocolName[i] <= '9') ||
                                            (protocolName[i] >= 'A' && protocolName[i] <= 'Z')))
        {
          protocolName = protocolName.erase(i, 1).c_str();
        }
      }

      std::string metaFilename;
      metaFilename = metaDir;
      metaFilename += "/";
      metaFilename += argv[2];
      metaFilename += "_";
      metaFilename += modality.c_str();
      metaFilename += "_";
      metaFilename += seriesNum.c_str();
      if (sequenceName.size() > 1)
      {
        metaFilename += "_";
        metaFilename += sequenceName.c_str();
      }
      if (protocolName.size() > 1)
      {
        metaFilename += "_";
        metaFilename += protocolName.c_str();
      }
      metaFilename += "_";
      if (imageIs3D)
      {
        metaFilename += coord;
      }
      else
      {
        metaFilename += "2D";
      }
      metaFilename += ".mha";

      writer = WriterType::New();
      writer->SetFileName(metaFilename.c_str());
      if (imageIs3D)
      {
        writer->SetInput(reader->GetOutput());
      }
      else
      {
        writer->SetInput(reader2D->GetOutput());
      }
      writer->SetUseCompression(true);

      try
      {
        writer->Update();
      }
      catch (itk::ExceptionObject & ex)
      {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }

      gdcm::FileHelper * fileReader;

      // Now Anonymize the DICOM object
      for (int fileNum = 0; fileNum < fileNames.size(); fileNum++)
      {
        gdcm::File * file = (*(nameGenerator->GetSeriesHelper()->GetCoherentFileList(seriesIdentifier)))[fileNum];

        std::string oldFilenameWithPath = file->GetFileName().c_str();

        std::string filename = oldFilenameWithPath.c_str();
        itksys::SystemTools::ConvertToUnixSlashes(filename);
        int len = filename.size();
        int start = len - 20;
        int split = filename.find_last_of("/");
        if (split == std::string::npos)
        {
          split = 0;
        }
        else
        {
          split++;
        }
        if (start < split)
        {
          start = split;
        }

        // Preserve a copy of the unchanged file
        std::string privateFilename = privateDir.c_str();
        privateFilename += "/";
        privateFilename += filename.substr(split, 255).c_str();

        itksys::SystemTools::CopyFileAlways(oldFilenameWithPath.c_str(), privateFilename.c_str());

        // Open the file and start changing it
        // file->SetLoadMode( gdcm::LD_ALL );
        // file->Load();

        fileReader = new gdcm::FileHelper(file);

        uint8_t * imageData = fileReader->GetImageData();

        std::string publicFilename = publicDir;
        publicFilename += "/";
        publicFilename += argv[2];
        publicFilename += "_";
        publicFilename += modality.c_str();
        publicFilename += "_";
        publicFilename += seriesNum.c_str();
        if (sequenceName.size() > 1)
        {
          publicFilename += "_";
          publicFilename += sequenceName.c_str();
        }
        if (protocolName.size() > 1)
        {
          publicFilename += "_";
          publicFilename += protocolName.c_str();
        }
        std::string  pos = file->GetEntryValue(0x0020, 0x0032);
        unsigned int splitX = pos.find_first_of("\\");
        unsigned int splitY = pos.find_first_of("\\", splitX + 1);
        unsigned int lenY = splitY - splitX - 1;
        unsigned int lenZ = pos.size() - splitY - 1;
        double       xPos = atof(pos.substr(0, splitY).c_str());
        double       yPos = atof(pos.substr(splitX + 1, lenY).c_str());
        double       zPos = atof(pos.substr(splitY + 1, lenZ).c_str());
        char         coord[80];
        sprintf(coord, "%0.2fx%0.2fx%0.2f", xPos, yPos, zPos);
        publicFilename += "_";
        publicFilename += coord;
        publicFilename += "_";
        publicFilename += filename.substr(start, 20).c_str();
        publicFilename += ".dcm";

        // StudyDate : 0x0008, 0x0020);
        std::string studyDate = file->GetEntryValue(0x0008, 0x0020).c_str();
        std::string newStudyDate = studyDate.c_str();
        newStudyDate[6] = '0';
        newStudyDate[7] = '1';

        // BirthDate : 0x0010, 0x0030);
        std::string birthDate = file->GetEntryValue(0x0010, 0x0030).c_str();
        std::string newBirthDate = birthDate.c_str();
        newBirthDate[4] = '0';
        newBirthDate[5] = '1';
        newBirthDate[6] = '0';
        newBirthDate[7] = '1';

        file->ClearAnonymizeList();
        // InstitutionName -> UNC-CH: MR Research Center: CADDLab
        // file->SetValEntry( "UNC-CH: MR Research Center: CADDLab",
        // 0x0008, 0x0080 );
        file->AddAnonymizeElement(0x0008, 0x0080, "UNC-CH: MR Research Center: CADDLab");

        // InstitutionAddress -> http://caddlab.rad.unc.edu
        // file->SetValEntry( "http://caddlab.rad.unc.edu",
        // 0x0008, 0x0081 );
        file->AddAnonymizeElement(0x0008, 0x0081, "http://caddlab.rad.unc.edu");

        // PatientsName -> argv[2]
        // file->SetValEntry( argv[2], 0x0010, 0x0010);
        file->AddAnonymizeElement(0x0010, 0x0010, argv[2]);

        // PatientsID -> argv[2]
        // file->SetValEntry( argv[2], 0x0010, 0x0020);
        file->AddAnonymizeElement(0x0010, 0x0020, argv[2]);

        gdcm::DictEntry * it = file->GetFirstEntry();
        while (it != NULL)
        {
          unsigned int group = it->GetGroup();
          unsigned int element = it->GetElement();
          bool         found = false;
          for (int i = 0; i < numID; i++)
          {
            if (group == groupID[i] && element == elementID[i])
            {
              gdcm::DictEntry * tmpIt = it;
              it = file->GetNextEntry();
              // file->RemoveEntry(tmpIt);
              file->AddAnonymizeElement(group, element, "");
              found = true;
              break;
            }
            else
            {
              // groupID is an ordered sequence - abort if past group
              if (group < groupID[i])
              {
                break;
              }
            }
          }
          if (!found)
          {
            // check for a date
            int pos = file->GetEntryValue(group, element).find(studyDate);
            if (pos != std::string::npos)
            {
              std::string newV = file->GetEntryValue(group, element).c_str();
              newV.replace(pos, studyDate.size(), newStudyDate);
              // file->SetValEntry(newV, group, element);
              file->AddAnonymizeElement(group, element, newV);
            }
            // check for a date
            pos = file->GetEntryValue(group, element).find(birthDate);
            if (pos != std::string::npos)
            {
              std::string newV = file->GetEntryValue(group, element).c_str();
              newV.replace(pos, birthDate.size(), newBirthDate);
              // file->SetValEntry(newV, group, element);
              file->AddAnonymizeElement(group, element, newV);
            }

            it = file->GetNextEntry();
          }
        }

        file->AnonymizeNoLoad();

        // fileReader->WriteDcmExplVR( publicFilename );

        file->ClearAnonymizeList();

        file->CloseFile();

        itksys::SystemTools::CopyFileAlways(oldFilenameWithPath.c_str(), publicFilename.c_str());

        itksys::SystemTools::RemoveFile(oldFilenameWithPath.c_str());
      }
      seriesItr++;
    }
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
