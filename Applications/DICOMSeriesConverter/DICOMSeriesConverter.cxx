/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#include <iostream>
#include <fstream>
#include <MetaIO/metaCommand.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkImageFileReader.h>
#include <itkImageSeriesReader.h>
#include <itkImageFileWriter.h>
#include <itksys/SystemTools.hxx>
#include <gdcm/src/gdcmFile.h>
#include <gdcm/src/gdcmFileHelper.h>
#include <gdcm/src/gdcmCommon.h>
#include <gdcm/src/gdcmDocEntry.h>
#include <gdcm/src/gdcmSerieHelper.h>

#include <stdio.h>

const int PatientNameGroup = 0x0010;
const int PatientNameElement = 0x0010;

const int StudyDateGroup = 0x0008;
const int StudyDateElement = 0x0020;

const int ModalityGroup = 0x0008;
const int ModalityElement = 0x0060;

const int SeriesNumberGroup = 0x0020;
const int SeriesNumberElement = 0x0011;

const int SequenceNameGroup = 0x0018;
const int SequenceNameElement = 0x0024;

const int SliceThicknessGroup = 0x0018;
const int SliceThicknessElement = 0x0050;

const int ProtocolNameGroup = 0x0018;
const int ProtocolNameElement = 0x1030;

const int SpacingGroup = 0x0028;
const int SpacingElement = 0x0030;

const int PositionGroup = 0x0020;
const int PositionElement = 0x0032;

const int StudyIDGroup = 0x0020;
const int StudyIDElement = 0x0010;

const std::string CleanString(std::string inString, int strLen)
{
  if(inString == gdcm::GDCM_UNFOUND)
    {
    inString = "";
    }
  else
    {
    for(int i=0; i<(int)inString.size(); i++)
      {
      while(i<(int)inString.size()
            && !(   (inString[i] >= 'a' && inString[i] <= 'z')
                 || (inString[i] >= '0' && inString[i] <= '9')
                 || (inString[i] >= 'A' && inString[i] <= 'Z')
                 || inString[i] == '_'
                 || inString[i] == ' '))
        {
        inString = inString.erase(i,1).c_str();
        }
      if(i<(int)inString.size() && inString[i] == ' ')
        {
        inString[i] = '_';
        }
      }

    for(int i=0; i<(int)inString.size(); i++)
      {
      while(i<(int)inString.size()
             && !(   (inString[i] >= 'a' && inString[i] <= 'z')
                  || (inString[i] >= '0' && inString[i] <= '9')
                  || (inString[i] >= 'A' && inString[i] <= 'Z')))
         {
         inString = inString.erase(i,1).c_str();
         }
      }

    if(strLen>0 && (int)inString.size()>strLen)
      {
      inString.erase(strLen, inString.size()-strLen);
      }
    }

  if(strLen>0 && (int)inString.size()<strLen)
    {
    int i = inString.size();
    while(i<strLen)
      {
      inString[i] = 'x';
      }
    }

  return inString;
}


int main( int argc, char* argv[] )
{

  MetaCommand command;

  command.SetOption("Split", "s", false,
                    "If <SequenceName>, split into <Number> slices / volume.");
  command.AddOptionField("Split", "SequenceName", MetaCommand::STRING, true);
  command.AddOptionField("Split", "Number", MetaCommand::STRING, true);

  command.SetOption("Recursive", "r", false,
                    "Process subdirectories recursively.");

  command.SetOption("IgnoreSeries", "i", false,
                    "Ignore series numbers; try to merge to a single volume.");

  command.SetOption("Anonymize", "a", false,
                    "Create anonymized dicom files using <newName>.");
  command.AddOptionField("Anonymize", "Name", MetaCommand::STRING, true);
  command.AddOptionField("Anonymize", "ID", MetaCommand::STRING, true);

  command.SetOption("CompressionOff", "c", false,
                    "Turn off metaImage pixel data compression.");

  command.SetOption("WriteRaw", "w", false,
                    "Split metaImage into header and raw files per image.");

  command.SetOption("OutputDirectory", "o", false,
                    "Write output files to the specified directory.");
  command.AddOptionField("OutputDirectory", "Directory",
                    MetaCommand::STRING, true);

  command.SetOption("OutputBasename", "b", false,
                    "Pre-append <basename> to each output filename.");
  command.AddOptionField("OutputBasename", "Basename",
                    MetaCommand::STRING, true);

  command.AddField("InputDirectory", "Directory to be processed",
                    MetaCommand::STRING, true);

  if(!command.Parse(argc, argv))
    {
    std::cerr << "Command line error" << std::endl;
    return 1;
    }

  bool optSplit = false;
  std::string splitSequenceName;
  int splitNumber = 0;
  if(command.GetOptionWasSet("Split"))
    {
    optSplit = true;
    splitSequenceName = command.GetValueAsString("Split", "SequenceName");
    splitNumber = command.GetValueAsInt("Split", "Number");
    }

  bool optAnonymize = false;
  std::string anonymizeName = "";
  std::string anonymizeID = "";
  if(command.GetOptionWasSet("Anonymize"))
    {
    optAnonymize = true;
    anonymizeName = command.GetValueAsString("Anonymize", "Name");
    anonymizeID = command.GetValueAsString("Anonymize", "ID");
    }

  bool optCompressionOff = false;
  if(command.GetOptionWasSet("CompressionOff"))
    {
    optCompressionOff = true;
    }
  bool optWriteRaw = false;
  if(command.GetOptionWasSet("WriteRaw"))
    {
    optWriteRaw = true;
    }
  bool optRecursive = false;
  if(command.GetOptionWasSet("Recursive"))
    {
    optRecursive = true;
    }
  bool optBasename = false;
  std::string basename;
  if(command.GetOptionWasSet("OutputBasename"))
    {
    optBasename = true;
    basename = command.GetValueAsString("OutputBasename", "Basename").c_str();
    }
  bool optOutputDirectory = false;
  std::string outputDirectory;
  if(command.GetOptionWasSet("OutputDirectory"))
    {
    optOutputDirectory = true;
    outputDirectory = command.GetValueAsString("OutputDirectory",
                                               "Directory").c_str();
    }

  std::string inputDirectory = command.GetValueAsString("InputDirectory")
                                      .c_str();

  //
  // Parameters to be used as template arguments
  //
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;

  //
  // Typedefs
  //
  typedef itk::Image< PixelType, Dimension >         ImageType;

  typedef itk::ImageFileReader< ImageType >          Reader2DType;
  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  typedef itk::GDCMImageIO                           ImageIOType;

  typedef std::vector< std::string >                 SeriesIdContainer;
  typedef std::vector< std::string >                 FileNamesContainer;
  typedef FileNamesContainer::iterator               FileNamesIterator;
  typedef itk::ImageFileWriter< ImageType >          WriterType;

  typedef itk::GDCMSeriesFileNames                   NamesGeneratorType;

  //
  // IO variables
  //
  ReaderType::Pointer     reader;
  Reader2DType::Pointer   reader2D;

  ImageIOType::Pointer    dicomIO;

  WriterType::Pointer     writer;

  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->SetRecursive( optRecursive );
  if( !command.GetOptionWasSet("IgnoreSeries") )
    {
    nameGenerator->AddSeriesRestriction( "0020|0012" );
    }
  try
    {
    nameGenerator->SetInputDirectory( inputDirectory );
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << "ERROR: 2D Read: " << ex << std::endl;
    return EXIT_FAILURE;
    }

  //
  // TAGS TO ANONYMIZE
  //
  unsigned int groupID[100];
  unsigned int elementID[100];
  unsigned int count = 0;
  //  InstitutionCodeSequence
  groupID[count] = 0x0008; elementID[count++] = 0x0082;
  //  ReferringPhysiciansName
  groupID[count] = 0x0008; elementID[count++] = 0x0090;
  //  ReferringPhysiciansAddress
  groupID[count] = 0x0008; elementID[count++] = 0x0092;
  //  ReferringPhysiciansTelephoneNumbers
  groupID[count] = 0x0008; elementID[count++] = 0x0094;
  //  ReferringPhysiciansIdentificationSequence
  groupID[count] = 0x0008; elementID[count++] = 0x0096;
  //  ResponsibleOrganization
  groupID[count] = 0x0008; elementID[count++] = 0x0116;
  //  StationName
  groupID[count] = 0x0008; elementID[count++] = 0x1010;
  //  InstitutionalDepartmentName
  groupID[count] = 0x0008; elementID[count++] = 0x1040;
  //  PhysiciansOfRecord
  groupID[count] = 0x0008; elementID[count++] = 0x1048;
  //  PhysiciansOfRecordIdentificationSequence
  groupID[count] = 0x0008; elementID[count++] = 0x1049;
  //  PerformingPhysiciansName
  groupID[count] = 0x0008; elementID[count++] = 0x1050;
  //  PerformingPhysicianIdentificationSequence
  groupID[count] = 0x0008; elementID[count++] = 0x1052;
  //  NameOfPhysiciansReadingStudy
  groupID[count] = 0x0008; elementID[count++] = 0x1060;
  //  NameOfPhysiciansReadingStudyIdentificationSequence
  groupID[count] = 0x0008; elementID[count++] = 0x1062;
  //  OperatorsName
  groupID[count] = 0x0008; elementID[count++] = 0x1070;
  //  OperatorIdentificationSequence
  groupID[count] = 0x0008; elementID[count++] = 0x1072;
 
  //  IssuerOfPatientID
  groupID[count] = 0x0010; elementID[count++] = 0x0021;
  //  PatientsInsurancePlanCodeSequence
  groupID[count] = 0x0010; elementID[count++] = 0x0050;
  //  OtherPatientsIDs
  groupID[count] = 0x0010; elementID[count++] = 0x1000;
  //  OtherPatientNames
  groupID[count] = 0x0010; elementID[count++] = 0x1001;
  //  PatientsBirthName
  groupID[count] = 0x0010; elementID[count++] = 0x1005;
  //  PatientsAddress
  groupID[count] = 0x0010; elementID[count++] = 0x1040;
  //  PatientsMothersBirthName
  groupID[count] = 0x0010; elementID[count++] = 0x1060;
  //  MilitaryRank
  groupID[count] = 0x0010; elementID[count++] = 0x1080;
  //  MedicalRecordLocator
  groupID[count] = 0x0010; elementID[count++] = 0x1090;
  //  PatientsTelephoneNumbers
  groupID[count] = 0x0010; elementID[count++] = 0x2154;
 
  //  DeviceSerialNumber
  groupID[count] = 0x0018; elementID[count++] = 0x1000;
 
  // RequestingPhysicianIdentificationSequence
  groupID[count] = 0x0032; elementID[count++] = 0x1031;
  // RequestingPhysician
  groupID[count] = 0x0032; elementID[count++] = 0x1032;
 
  // AdmissionID
  groupID[count] = 0x0038; elementID[count++] = 0x0010;
  // IssuerOfAdmissionID
  groupID[count] = 0x0038; elementID[count++] = 0x0011;
  // PatientsInstitutionResidence
  groupID[count] = 0x0038; elementID[count++] = 0x0400;
 
  // ScheduledPerformingPhysiciansName
  groupID[count] = 0x0040; elementID[count++] = 0x0006;
  // ScheduledPerformingPhysiciansIdentificationSequence
  groupID[count] = 0x0040; elementID[count++] = 0x000B;
  // PerformedLocation
  groupID[count] = 0x0040; elementID[count++] = 0x0243;
  // NamesOfIntendedRecipientsOfResults
  groupID[count] = 0x0040; elementID[count++] = 0x1010;
  // IntendedRecipientsOfResultsIdentificationSequence
  groupID[count] = 0x0040; elementID[count++] = 0x1011;
  // PersonIdentificationCodeSequence
  groupID[count] = 0x0040; elementID[count++] = 0x1101;
  // PersonAddress
  groupID[count] = 0x0040; elementID[count++] = 0x1102;
  // PersonTelephoneNumbers
  groupID[count] = 0x0040; elementID[count++] = 0x1103;
  // OrderEnteredBy
  groupID[count] = 0x0040; elementID[count++] = 0x2008;
  // OrderEnterersLocation
  groupID[count] = 0x0040; elementID[count++] = 0x2009;
  // OrderCallbackPhoneNumber
  groupID[count] = 0x0040; elementID[count++] = 0x2010;
  // HumanPerformersOrganization
  groupID[count] = 0x0040; elementID[count++] = 0x4036;
  // HumanPerformersName
  groupID[count] = 0x0040; elementID[count++] = 0x4037;
  // VerifyingObserverName
  groupID[count] = 0x0040; elementID[count++] = 0xA075;
  // PersonName
  groupID[count] = 0x0040; elementID[count++] = 0xA123;
 
  // PhysicianApprovingInterpretation
  groupID[count] = 0x4008; elementID[count++] = 0x0114;

  // Number of tags to anonymize
  int numID = count;

  //
  // Begin forming the output filenames
  //
  std::string prependFilename = "";
  if(optOutputDirectory)
    {
    prependFilename = outputDirectory.c_str();
    itksys::SystemTools::MakeDirectory(outputDirectory.c_str());
    prependFilename += "/";
    }
  if(optBasename)
    {
    prependFilename += basename.c_str();
    }

  try
    {
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    FileNamesContainer fileNames;

    std::string seriesIdentifier;
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

    int  unsplitNumberOfFileNames;
    bool processingSplit = false;
    int  splitVolumeStart = 0;
    while( seriesItr != seriesEnd )
      {
      seriesIdentifier = *seriesItr;

      dicomIO = ImageIOType::New();

      bool imageIs3D;
      fileNames = nameGenerator->GetFileNames( seriesIdentifier );
      unsplitNumberOfFileNames = fileNames.size();

      if(optSplit)
        {
        if(processingSplit)
          {
          fileNames = nameGenerator->GetFileNames( seriesIdentifier );
          splitVolumeStart += splitNumber;
          FileNamesIterator  fileNamesStart;
          FileNamesIterator  fileNamesEnd;
          fileNamesStart = fileNames.begin();
          fileNamesEnd = fileNames.begin() + splitVolumeStart;
          fileNames.erase(fileNamesStart, fileNamesEnd);
          if(splitVolumeStart+splitNumber < unsplitNumberOfFileNames)
            {
            fileNamesStart = fileNames.begin() + splitNumber;
            fileNamesEnd = fileNames.end();
            fileNames.erase(fileNamesStart, fileNamesEnd);
            }
          }
        else
          {
          gdcm::File *file = (*(nameGenerator->GetSeriesHelper()
                                             ->GetSingleSerieUIDFileSet(
                                                    seriesIdentifier )))[0];
          std::string sequenceName = file->GetEntryValue( SequenceNameGroup,
                                                          SequenceNameElement );
          if( sequenceName == splitSequenceName
              && (int)fileNames.size() > splitNumber)
            {
            processingSplit = true;
            splitVolumeStart = 0;
            FileNamesIterator  fileNamesStart;
            FileNamesIterator  fileNamesEnd;
            fileNamesStart = fileNames.begin() + splitNumber;
            fileNamesEnd = fileNames.end();
            fileNames.erase(fileNamesStart, fileNamesEnd);
            }
          else
            {
            processingSplit = false;
            splitVolumeStart = 0;
            }
          }
        }
      if(fileNames.size() < 2)
        {
        imageIs3D = false;

        reader2D = Reader2DType::New();
        reader2D->SetFileName( fileNames.begin()->c_str() );

        try
          {
          reader2D->Update();
          }
        catch (itk::ExceptionObject &ex)
          {
          std::cout << "ERROR: 2D Read: " << ex << std::endl;
          return EXIT_FAILURE;
          }
        }
      else
        {
        imageIs3D = true;

        reader = ReaderType::New();
        reader->SetImageIO( dicomIO );
  
        reader->SetFileNames( fileNames );
    
        try
          {
          reader->Update();
          }
        catch (itk::ExceptionObject &ex)
          {
          std::cout << "ERROR: 3D Read: " << ex << std::endl;
          return EXIT_FAILURE;
          }
        }
  
      gdcm::File *file = (*(nameGenerator->GetSeriesHelper()
                                         ->GetSingleSerieUIDFileSet(
                                             seriesIdentifier )))[0];

      std::string patientName;
      if(!optAnonymize)
        {
        patientName = file->GetEntryValue( PatientNameGroup,
                                           PatientNameElement );
        }
      else
        {
        patientName = anonymizeName;
        }
      patientName = CleanString(patientName, 0);

      std::string studyDate = file->GetEntryValue( StudyDateGroup,
                                                   StudyDateElement );
      studyDate = CleanString(studyDate, 0);

      std::string modality = file->GetEntryValue( ModalityGroup,
                                                  ModalityElement );
      modality = CleanString(modality, 0);

      std::string studyID = file->GetEntryValue( StudyIDGroup,
                                                 StudyIDElement );
      studyID = CleanString(studyID, 0);

      std::string seriesNumber = file->GetEntryValue( SeriesNumberGroup,
                                                      SeriesNumberElement );
      seriesNumber = CleanString(seriesNumber, 0);

      std::string sequenceName;
      if( modality == "CT")
        {
        // No real meaning to sequenceName in CT, so use sliceThickness
        sequenceName = file->GetEntryValue( SliceThicknessGroup,
                                            SliceThicknessElement );
        if( sequenceName == gdcm::GDCM_UNFOUND )
          {
          sequenceName = "";
          }
        else
          {
          double sliceThickness = atof(sequenceName.c_str());
          char st[80];
          sprintf(st, "%0.2f", sliceThickness);
          sequenceName = st;
          }
        }
      else
        {
        sequenceName = file->GetEntryValue( SequenceNameGroup,
                                            SequenceNameElement );
        sequenceName = CleanString(sequenceName, 0);
        }

      std::string protocolName = file->GetEntryValue( ProtocolNameGroup,
                                                      ProtocolNameElement );
      protocolName = CleanString(protocolName, 0);

      char coord[80];
      if( imageIs3D )
        {
        std::string spacing = file->GetEntryValue( SpacingGroup,
                                                   SpacingElement );
        int split = spacing.find_first_of("\\");
        // int len = spacing.size()-split-1; UNUSED
        double xSpacing = atof(spacing.substr(0,split).c_str());
        // double ySpacing = atof(spacing.substr(split+1,len).c_str()); UNUSED
        std::string pos = file->GetEntryValue( PositionGroup, PositionElement );
        int splitX = pos.find_first_of("\\");
        int splitY = pos.find_first_of("\\", splitX+1);
        int lenY = splitY-splitX-1;
        int lenZ = pos.size()-splitY-1;
        double xPos = atof(pos.substr(0, splitY).c_str());
        double yPos = atof(pos.substr(splitX+1, lenY).c_str());
        double zPos = atof(pos.substr(splitY+1, lenZ).c_str());
        file = (*(nameGenerator->GetSeriesHelper()
                               ->GetSingleSerieUIDFileSet(
                                        seriesIdentifier )))[1];
        pos = file->GetEntryValue( PositionGroup, PositionElement );
        splitX = pos.find_first_of("\\");
        splitY = pos.find_first_of("\\", splitX+1);
        lenY = splitY-splitX-1;
        lenZ = pos.size()-splitY-1;
        xPos = (xPos - atof(pos.substr(0, splitY).c_str()));
        yPos = (yPos - atof(pos.substr(splitX+1, lenY).c_str()));
        zPos = (zPos - atof(pos.substr(splitY+1, lenZ).c_str()));
        double zSpacing = sqrt(xPos*xPos + yPos*yPos + zPos*zPos);
        if(processingSplit)
          {
          sprintf(coord, "%0.2f_%0.2f_%d_%03d", xSpacing, zSpacing,
                                                (int)(fileNames.size()),
                                                (int)(splitVolumeStart
                                                     /(double)splitNumber));
          }
        else
          {
          sprintf(coord, "%0.2f_%0.2f_%d", xSpacing, zSpacing,
                                           (int)(fileNames.size()));
          }
        }
      else
        {
        coord[0] = '\0';
        }
  
      if(protocolName == sequenceName)
        {
        // if protocol name and sequence name are the same - use study
        // ID instead of sequence name
        sequenceName = file->GetEntryValue( StudyIDGroup, StudyIDElement );
        sequenceName = CleanString(sequenceName, 0);
        }

      std::string metaFilename;
      metaFilename = prependFilename.c_str();
      metaFilename += patientName.c_str();
      metaFilename += "_";
      metaFilename += studyDate.c_str();
      metaFilename += "_";
      metaFilename += sequenceName.c_str();
      metaFilename += "_";
      metaFilename += seriesNumber.c_str();
      metaFilename += "_";
      metaFilename += protocolName.c_str();
      metaFilename += "_";
      if( imageIs3D )
        {
        metaFilename += coord;
        }
      else
        {
        metaFilename += "2D";
        }

      std::string filename = metaFilename.c_str();

      if(optWriteRaw)
        {
        metaFilename += ".mhd";
        }
      else
        {
        metaFilename += ".mha";
        }

      writer = WriterType::New();
      writer->SetFileName( metaFilename.c_str() );
      if( imageIs3D )
        {
        writer->SetInput( reader->GetOutput() );
        }
      else
        {
        writer->SetInput( reader2D->GetOutput() );
        }
      if(!optCompressionOff)
        {
        writer->SetUseCompression(true);
        }
  
      try
        {
        writer->Update();
        }
      catch (itk::ExceptionObject &ex)
        {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
        }

      std::string mappingFilename = filename.c_str();
      mappingFilename += "_DCM";
      itksys::SystemTools::MakeDirectory(mappingFilename.c_str());
      mappingFilename += "/DCM_To_MHA.txt";
      std::ofstream mappingFilenameStream(mappingFilename.c_str(),
                                          std::ios::out);

      //gdcm::FileHelper * fileReader;

      // Now Anonymize the DICOM object
      for(int fileNum=splitVolumeStart;
          fileNum<(int)(splitVolumeStart+fileNames.size());
          fileNum++)
        {
        char fileNumString[100];
        sprintf(fileNumString, "%04d", fileNum-splitVolumeStart);

        gdcm::File *newfile = (*(nameGenerator->GetSeriesHelper()
                                           ->GetSingleSerieUIDFileSet(
                                               seriesIdentifier )))[fileNum];

        std::string oldFilenameWithPath = newfile->GetFileName().c_str();
        std::string newFilename = oldFilenameWithPath.c_str();

        itksys::SystemTools::ConvertToUnixSlashes( newFilename );
        int len = newFilename.size();
        int start = len - 20;
        int split = newFilename.find_last_of("/");
        if(split == (int)std::string::npos)
          {
          split = 0;
          }
        else
          {
          split++;
          }
        if(start < split)
          {
          start = split;
          }
        std::string oldFilename = oldFilenameWithPath.substr(start,
                                      oldFilenameWithPath.size()-start).c_str();
        std::string oldFilenamePath = oldFilenameWithPath.substr(0,
                                          start).c_str();

        std::string newDicomFilename = filename;
        newDicomFilename += "_DCM/";
        newDicomFilename += oldFilename;

        newfile->CloseFile();

        itksys::SystemTools::CopyFileAlways(oldFilenameWithPath.c_str(),
                                            newDicomFilename.c_str());

        mappingFilenameStream << oldFilenameWithPath << " "
                              << newDicomFilename << "  "
                              << metaFilename << "  "
                              << fileNumString << std::endl;

        if(optAnonymize)
          {
          // Open the newfile and start changing it
          //newfile->SetLoadMode( gdcm::LD_ALL );
          //newfile->Load();

          newfile = new gdcm::File;
          newfile->SetFileName(newDicomFilename);
          newfile->Load();

          // fileReader = new gdcm::FileHelper( newfile );

          // uint8_t * imageData = fileReader->GetImageData();
    
          //StudyDate : 0x0008, 0x0020);
          //std::string studyDate = newfile->GetEntryValue(0x0008, 0x0020).c_str();
          //std::string newStudyDate = studyDate.c_str();
          //newStudyDate[6] = '0';
          //newStudyDate[7] = '1';
  
          //BirthDate : 0x0010, 0x0030);
          std::string birthDate = newfile->GetEntryValue( 0x0010, 0x0030 ).c_str();
          std::string newBirthDate = birthDate.c_str();
          newBirthDate[4] = '0';
          newBirthDate[5] = '1';
          newBirthDate[6] = '0';
          newBirthDate[7] = '1';
      
          newfile->ClearAnonymizeList();
          // InstitutionName -> Kitware Anonymization Software
          //newfile->SetValEntry( "Kitware Anonymization Software",
                             // 0x0008, 0x0080 );
          newfile->AddAnonymizeElement( 0x0008, 0x0080,
                                     "Kitware Anonymization Software" );
                            
          // InstitutionAddress -> http://kitware.com
          //newfile->SetValEntry( "http://kitware.com",
                             // 0x0008, 0x0081 );
          newfile->AddAnonymizeElement( 0x0008, 0x0081,
                                     "http://kitware.com" );
      
          // PatientsName -> argv[2]
          //newfile->SetValEntry( argv[2], 0x0010, 0x0010);
          newfile->AddAnonymizeElement( 0x0010, 0x0010, anonymizeName );
      
          // PatientsID -> argv[2]
          //newfile->SetValEntry( argv[2], 0x0010, 0x0020);
          newfile->AddAnonymizeElement( 0x0010, 0x0020, anonymizeID );
      
          gdcm::DocEntry * it = newfile->GetFirstEntry();
          while(it != NULL)
            {
            unsigned int group = it->GetGroup();
            unsigned int element = it->GetElement();
            bool found = false;
            for(int i = 0; i<numID; i++)
              {
              if(group == groupID[i]
              && element == elementID[i])
                {
                //gdcm::DocEntry * tmpIt = it;
                it = newfile->GetNextEntry();
                //newfile->RemoveEntry(tmpIt);
                newfile->AddAnonymizeElement( group, element, "");
                found = true;
                break;
                }
              else
                {
                // groupID is an ordered sequence - abort if past group
                if(group < groupID[i])
                  {
                  break;
                  }
                }
              }
            if(!found)
              {
              // check for a date
              //int pos = newfile->GetEntryValue(group, element).find( studyDate );
              //if(pos != std::string::npos)
              //  {
              //  std::string newV = newfile->GetEntryValue(group,element).c_str();
              //  newV.replace(pos, studyDate.size(), newStudyDate);
              //  newfile->SetValEntry(newV, group, element);
              //  newfile->AddAnonymizeElement( group, element, newV );
              //  }
              //
              // check for a date
              int pos = newfile->GetEntryValue(group, element).find( birthDate );
              if(pos != (int)std::string::npos)
                {
                std::string newV = newfile->GetEntryValue(group, element).c_str();
                newV.replace(pos, birthDate.size(), newBirthDate);
                //newfile->SetValEntry(newV, group, element);
                newfile->AddAnonymizeElement( group, element, newV );
                }
  
              it = newfile->GetNextEntry();
              }
            }
  
          newfile->AnonymizeNoLoad();
  
          //fileReader->WriteDcmExplVR( publicFilename );
          
          newfile->ClearAnonymizeList();
  
          newfile->CloseFile();
          }
        }
      if(!processingSplit)
        {
        seriesItr++;
        }
      else
        {
        if(splitVolumeStart+splitNumber >= unsplitNumberOfFileNames)
          {
          processingSplit = false;
          splitVolumeStart = 0;
          seriesItr++;
          }
        }
      }
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << "ERROR : " << ex << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
