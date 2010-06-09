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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// The following three should be used in every CLI application
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"

// Includes specific to this CLI application
#include "igtlPositionMessage.h"
#include "igtlServerSocket.h"
#include "igtlOSUtil.h"

// Must include CLP before including Helper functions
#include "ReceivePositionMessageUsingOpenIGTLinkCLP.h"

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  //Establish connection
  igtl::ClientSocket::Pointer clientSocket;
  clientSocket = igtl::ClientSocket::New();
  int connectionStatus = clientSocket->ConnectToServer(hostname.c_str(),port);

  if (connectionStatus < 0)
    {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(1);
    }
  else
    {
    std::cout << "Connected to the server :\t" << hostname << "\t" << port << std::endl; 
    }

  //------------------------------------------------------------
  // Create a message buffer to receive header
  igtl::MessageHeader::Pointer headerMessage;
  headerMessage = igtl::MessageHeader::New();
  
  while (1)
    {
    //------------------------------------------------------------
    // loop
    for (int i = 0; i < 100; i ++)
      {
      // Initialize receive buffer
      headerMessage->InitPack();
      
      // Receive generic header from the socket
      int receiveReturnValue = clientSocket->Receive(headerMessage->GetPackPointer(), headerMessage->GetPackSize());
      if (receiveReturnValue == 0)
        {
        clientSocket->CloseSocket();
        exit(0);
        }
      if (receiveReturnValue != headerMessage->GetPackSize())
        {
        continue;
        }
      
      // Deserialize the header
      headerMessage->Unpack();
      
      // Check for position message type and receive data body
      if (strcmp(headerMessage->GetDeviceType(), "POSITION") == 0)
        {
        std::cerr << "Receiving POSITION data type." << std::endl;
        
        // Create a message buffer to receive transform data
        igtl::PositionMessage::Pointer positionMessage;
        positionMessage = igtl::PositionMessage::New();
        positionMessage->SetMessageHeader(headerMessage);
        positionMessage->AllocatePack();
        
        // Receive position position data from the socket
        clientSocket->Receive(positionMessage->GetPackBodyPointer(), positionMessage->GetPackBodySize());
        
        // Deserialize the transform data
        // If you want to skip CRC check, call Unpack() without argument.
        int c = positionMessage->Unpack(1);
        
        if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
          {
          // Retrive the transform data
          float position[3];
          float quaternion[4];

          positionMessage->GetPosition(position);
          positionMessage->GetQuaternion(quaternion);

          std::cerr << "position   = (" << position[0] << ", " << position[1] << ", " << position[2] << ")" << std::endl;
          std::cerr << "quaternion = (" << quaternion[0] << ", " << quaternion[1] << ", "
                    << quaternion[2] << ", " << quaternion[3] << ")" << std::endl << std::endl;

          }
        }
      }
    }

  //------------------------------------------------------------
  // Close connection 
  
  clientSocket->CloseSocket();

  return EXIT_SUCCESS;
}


